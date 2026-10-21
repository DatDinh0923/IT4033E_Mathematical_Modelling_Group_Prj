# =============================================================================
#  run_experiments.jl — numerical experiments on SetA
#
#  usage:  julia --project=. run_experiments.jl <mode> <timelimit> <runs> [instances...]
#     mode = p0   : Step 6, period 0 only (no budget), e.g. timelimit 1800
#     mode = p01  : Step 9, periods 0 and 1 with budget, e.g. timelimit 3600
#  Results are appended to results/<mode>-<rmethod>.csv (one line per run) and
#  the solution of run 1 is written to results/<inst>-<mode>-<rmethod>-srpaths.json
# =============================================================================

include(joinpath(@__DIR__, "src", "data.jl"))
include(joinpath(@__DIR__, "src", "model.jl"))
using Printf, Statistics

const DATADIR = get(ENV, "DATADIR",
    joinpath(@__DIR__, "..", "SOICT_mathematical_modelling", "project", "setA"))
const OUTDIR  = joinpath(@__DIR__, "results")
mkpath(OUTDIR)

mode      = length(ARGS) >= 1 ? ARGS[1] : "p0"
timelimit = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1800.0
runs      = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
insts     = length(ARGS) >= 4 ? ARGS[4:end] : ["setA-$(lpad(k,2,'0'))" for k in 1:12]
rmethod   = Symbol(get(ENV, "RMETHOD", "pathcount"))
periods   = mode == "p0" ? [0] : [0, 1]

csv = joinpath(OUTDIR, "$mode-$rmethod.csv")
isfile(csv) || open(csv, "w") do io
    println(io, "instance,run,V,A,D,maxseg,budget,vars,cons,t_data,t_build,status,obj,lb,gap,nodes,cpu,mlu_check,dist_check,maxrss_MB")
end

for name in insts
    for run in 1:runs
        GC.gc()
        t_data = @elapsed inst = read_instance(DATADIR, name; periods = periods, rmethod = rmethod)
        t_build = @elapsed (m, x, λ, S) = build_model(inst; periods = periods,
                                                       budget = mode != "p0",
                                                       timelimit = timelimit)
        nv_, nc_ = num_variables(m), num_constraints(m; count_variable_in_set_constraints = false)
        res = solve!(m)

        # independent verification of the reported objective
        chk = NaN; dchk = -1
        if !isnan(res.obj)
            Wt = Dict((d, t) => waypoints(inst, x, S, d, t)
                      for d in eachindex(inst.demands), t in periods)
            chk = maximum(mlu_from_waypoints(inst, Dict(d => Wt[(d, t)] for d in eachindex(inst.demands)), t)
                          for t in periods)
            if mode != "p0"
                dchk = scheme_distance(inst, Dict(d => Wt[(d, 0)] for d in eachindex(inst.demands)),
                                             Dict(d => Wt[(d, 1)] for d in eachindex(inst.demands)))
            end
            run == 1 && open(io -> write_srpaths(io, inst, Wt),
                             joinpath(OUTDIR, "$name-$mode-$rmethod-srpaths.json"), "w")
        end
        bud = mode == "p0" ? "" : string(get(inst.budget, 1, ""))
        line = @sprintf("%s,%d,%d,%d,%d,%d,%s,%d,%d,%.2f,%.2f,%s,%.6f,%.6f,%.6f,%d,%.2f,%.6f,%d,%d",
                        name, run, inst.n, length(inst.arcs), length(inst.demands), inst.maxseg, bud,
                        nv_, nc_, t_data, t_build, res.status, res.obj, res.lb, res.gap,
                        res.nodes, res.cpu, chk, dchk, Sys.maxrss() ÷ 2^20)
        println(line); flush(stdout)
        open(io -> println(io, line), csv, "a")
    end
end