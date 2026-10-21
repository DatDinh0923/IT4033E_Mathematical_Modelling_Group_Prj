# =============================================================================
#  two_stage.jl — two-stage method for the two-period problem (readme Step 9)
#
#  Stage 1 = the period-0 run of Step 6: routing P0, MLU λ0, lower bound LB0.
#  Stage 2 = period-1 model with P0 FIXED. The budget constraint (4) becomes
#            linear in x1 alone (no z variables):
#              Σ_d [ Σ_{e∈P0_d} (1 - x1[d,e]) + Σ_{e∉P0_d} x1[d,e] ] ≤ κ(1)
#            and λ ≥ λ0 makes the objective max(λ0, λ1), the same objective
#            as the joint model. Stage 2 has the size of the period-0 model.
#
#  Validity: (P0, P1) is feasible for the joint model, so UB = max(λ0, λ1) is
#  an upper bound; LB0 (period-0 bound) is a lower bound of the joint problem.
#  gap_joint = 1 - LB0 / UB; gap_joint = 0 proves joint optimality.
#
#  usage: julia --project=. two_stage.jl <timelimit_s> <runs> [instances...]
#  needs results/p0-<rmethod>.csv and results/<inst>-p0-<rmethod>-srpaths.json
# =============================================================================

include(joinpath(@__DIR__, "src", "data.jl"))
include(joinpath(@__DIR__, "src", "model.jl"))
using Printf

const DATADIR = get(ENV, "DATADIR",
    joinpath(@__DIR__, "..", "SOICT_mathematical_modelling", "project", "setA"))
const OUTDIR  = joinpath(@__DIR__, "results")

timelimit = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3600.0
runs      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
insts     = length(ARGS) >= 3 ? ARGS[3:end] : ["setA-$(lpad(k,2,'0'))" for k in 1:9]
rmethod   = Symbol(get(ENV, "RMETHOD", "pathcount"))

# ---- period-0 results (first CSV line of each instance) ----
function p0_row(name)
    f = joinpath(OUTDIR, "p0-$rmethod.csv")
    isfile(f) || error("missing $f: run Step 6 first")
    L = readlines(f); hdr = split(L[1], ",")
    for l in L[2:end]
        r = split(l, ","); r[1] == name || continue
        get(c) = parse(Float64, r[findfirst(==(c), hdr)])
        return (obj = get("obj"), lb = get("lb"), cpu = get("cpu"), status = r[findfirst(==("status"), hdr)])
    end
    error("no period-0 result for $name in $f")
end

function p0_waypoints(inst, name)
    f = joinpath(OUTDIR, "$name-p0-$rmethod-srpaths.json")
    isfile(f) || error("missing $f: run Step 6 first")
    W = Dict(d => Int[] for d in eachindex(inst.demands))
    for e in JSON3.read(read(f, String))["srpaths"]
        Int(e["t"]) == 0 || continue
        W[Int(e["d"]) + 1] = [inst.id_to_vertex[Int(w)] for w in e["w"]]
    end
    return W
end

segset(dem, w) = (q = [dem.s; w; dem.t]; Set((q[k], q[k+1]) for k in 1:length(q)-1))

csv = joinpath(OUTDIR, "p01ts-$rmethod.csv")
isfile(csv) || open(csv, "w") do io
    println(io, "instance,run,V,A,D,budget,vars,cons,t_data,t_build,status2,obj2,lb2,gap2,nodes,cpu2,",
                "p0_status,p0_cpu,lambda0,lambda1,lb_p0,ub_joint,gap_joint,proven,dist_check,maxrss_MB")
end

for name in insts, run in 1:runs
    GC.gc()
    t_data = @elapsed inst = read_instance(DATADIR, name; periods = [0, 1], rmethod = rmethod)
    p0 = p0_row(name)
    W0 = p0_waypoints(inst, name)
    λ0 = mlu_from_waypoints(inst, W0, 0)
    κ  = inst.budget[1]
    seg0 = [segset(dem, W0[d]) for (d, dem) in enumerate(inst.demands)]

    t_build = @elapsed begin
        # period-1 model: constraints (1)-(3) for t = 1, presolve on (no z rows)
        m, x, λ, S = build_model(inst; periods = [1], budget = false,
                                 timelimit = timelimit, warmstart = false)
        D = eachindex(inst.demands); S1 = S[1]; S1set = Set(S1)
        @constraint(m, λ >= λ0)
        # budget (4) with P0 fixed: changes = const - Σ_{e∈P0} x + Σ_{e∉P0} x
        nconst = sum(length(seg0[d]) for d in D)
        @constraint(m, nconst
                       - sum(x[1][d, e] for d in D for e in S1 if e in seg0[d])
                       + sum(x[1][d, e] for d in D for e in S1 if !(e in seg0[d])) <= κ)
        # MIP start: keep P0 in period 1 (0 changes) when all its segments exist in G(1)
        if all(e in S1set for d in D for e in seg0[d])
            for d in D, e in S1
                set_start_value(x[1][d, e], e in seg0[d] ? 1.0 : 0.0)
            end
            set_start_value(λ, max(λ0, mlu_from_waypoints(inst, W0, 1)))
        end
    end
    nv_, nc_ = num_variables(m), num_constraints(m; count_variable_in_set_constraints = false)
    res = solve!(m)

    λ1 = NaN; dchk = -1; ub = NaN; gapj = NaN; proven = "no"
    if !isnan(res.obj)
        W1 = Dict(d => waypoints(inst, x, S, d, 1) for d in D)
        λ1 = mlu_from_waypoints(inst, W1, 1)
        dchk = scheme_distance(inst, W0, W1)
        ub = max(λ0, λ1)
        gapj = ub > 0 ? max(0.0, 1.0 - p0.lb / ub) : NaN
        proven = gapj <= 1e-6 ? "yes" : "no"
        run == 1 && open(joinpath(OUTDIR, "$name-p01ts-$rmethod-srpaths.json"), "w") do io
            write_srpaths(io, inst, merge(Dict((d, 0) => W0[d] for d in D), Dict((d, 1) => W1[d] for d in D)))
        end
    end
    line = @sprintf("%s,%d,%d,%d,%d,%d,%d,%d,%.2f,%.2f,%s,%.6f,%.6f,%.6f,%d,%.2f,%s,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%d,%d",
                    name, run, inst.n, length(inst.arcs), length(inst.demands), κ, nv_, nc_, t_data, t_build,
                    res.status, res.obj, res.lb, res.gap, res.nodes, res.cpu,
                    p0.status, p0.cpu, λ0, λ1, p0.lb, ub, gapj, proven, dchk, Sys.maxrss() ÷ 2^20)
    println(line); flush(stdout)
    open(io -> println(io, line), csv, "a")
end
