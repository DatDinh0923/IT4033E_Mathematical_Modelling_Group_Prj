# Steps 6 and 9 — numerical experiments.
#
#   julia --project=. run.jl --model t0  [options] 1 2 3 4     # Step 6: period 0
#   julia --project=. run.jl --model t01 [options] 1 2 3 4     # Step 9: periods 0 and 1
#
# options: --time-limit S (default 1800 for t0, 3600 for t01)   --seed K
#          --mip-lp-solver choose|simplex|ipm  (LP algorithm inside branch-and-bound;
#                          "ipm" solves the root LP of setA-03/04 where dual simplex stalls)
#          --no-warm-start   --formulation aggregated|direct (t0 only)
#          --data-dir DIR (default ../setA)   --results-dir DIR (default results)
#          --out FILE (default <results-dir>/<model>_results.csv)
#
# One CSV row per instance is appended to the output file; the routing found
# is written to <results-dir>/solutions/<instance>-<model>-srpaths.json.

using Dates
using JuMP

include(joinpath(@__DIR__, "src", "OrangeSR.jl"))
using .OrangeSR

function parse_args(args)
    opts = Dict{Symbol, Any}(:model => "t0", :time_limit => nothing, :seed => nothing,
                             :warm_start => true, :formulation => :aggregated,
                             :mip_lp_solver => "choose", :out => nothing,
                             :data_dir => joinpath(@__DIR__, "..", "setA"),
                             :results_dir => joinpath(@__DIR__, "results"))
    instances = Int[]
    k = 1
    while k <= length(args)
        a = args[k]
        if a == "--model"
            opts[:model] = args[k += 1]
        elseif a == "--time-limit"
            opts[:time_limit] = parse(Float64, args[k += 1])
        elseif a == "--seed"
            opts[:seed] = parse(Int, args[k += 1])
        elseif a == "--no-warm-start"
            opts[:warm_start] = false
        elseif a == "--formulation"
            opts[:formulation] = Symbol(args[k += 1])
        elseif a == "--mip-lp-solver"
            opts[:mip_lp_solver] = args[k += 1]
        elseif a == "--data-dir"
            opts[:data_dir] = args[k += 1]
        elseif a == "--results-dir"
            opts[:results_dir] = args[k += 1]
        elseif a == "--out"
            opts[:out] = args[k += 1]
        else
            push!(instances, parse(Int, a))
        end
        k += 1
    end
    opts[:model] in ("t0", "t01") || error("--model must be t0 or t01")
    opts[:time_limit] === nothing &&
        (opts[:time_limit] = opts[:model] == "t0" ? 1800.0 : 3600.0)
    opts[:out] === nothing &&
        (opts[:out] = joinpath(opts[:results_dir], "$(opts[:model])_results.csv"))
    return instances, opts
end

const COLUMNS = Dict(
    "t0" => [:instance, :n_vertices, :n_links, :n_demands, :n_demands_model, :max_seg,
             :n_binary_vars, :n_constraints, :baseline_mlu, :objective, :best_bound, :gap,
             :status, :cpu_time_s, :highs_time_s, :bb_nodes, :mlu_check, :dropped_segments,
             :data_time_s, :build_time_s, :peak_rss_gb, :time_limit_s, :seed, :formulation,
             :mip_lp_solver, :warm_start, :threads, :date],
    "t01" => [:instance, :n_vertices, :n_links, :n_links_down, :n_demands, :n_demands_model,
              :max_seg, :budget, :n_binary_vars, :n_constraints, :baseline_mlu, :objective,
              :best_bound, :gap, :status, :cpu_time_s, :highs_time_s, :bb_nodes, :mlu0_check,
              :mlu1_check, :budget_used, :data_time_s, :build_time_s, :peak_rss_gb,
              :time_limit_s, :seed, :mip_lp_solver, :warm_start, :threads, :date])

function append_row(file, columns, row)
    mkpath(dirname(file))
    new = !isfile(file)
    open(file, "a") do io
        new && println(io, join(columns, ","))
        println(io, join((replace(string(row[c]), "," => ";") for c in columns), ","))
    end
end

# Fields shared by both result tables
run_info(opts, res) = (status = res.status, objective = res.objective, best_bound = res.bound,
                       gap = res.gap, cpu_time_s = res.cpu_time, highs_time_s = res.solve_time,
                       bb_nodes = res.nodes, peak_rss_gb = Sys.maxrss() / 2^30,
                       time_limit_s = opts[:time_limit], seed = something(opts[:seed], "default"),
                       mip_lp_solver = opts[:mip_lp_solver], warm_start = opts[:warm_start],
                       threads = Threads.nthreads(), date = now())

function prepare!(model, opts, silent)
    silent && set_silent(model)
    set_attribute(model, "mip_lp_solver", opts[:mip_lp_solver])
end

solution_file(opts, name) =
    (mkpath(joinpath(opts[:results_dir], "solutions"));
     joinpath(opts[:results_dir], "solutions", "$name-$(opts[:model])-srpaths.json"))

function run_instance_t0(dir, name, opts; silent = false)
    data_time = @elapsed inst = load_instance_t0(dir, name)
    baseline = mlu(inst, direct_paths(inst))

    build_time = @elapsed m = build_model_t0(inst; formulation = opts[:formulation],
                                             warm_start = opts[:warm_start])
    prepare!(m.model, opts, silent)
    println("[$name] data $(round(data_time; digits = 2)) s, build $(round(build_time; digits = 2)) s, ",
            "$(length(m.x)) binary variables"); flush(stdout)

    res = solve_t0!(m; time_limit = opts[:time_limit], seed = opts[:seed])

    check, dropped = NaN, -1
    if res.status in ("Optimal", "Not Optimal")
        paths, dropped = extract_paths(m)
        check = mlu(inst, paths)
        write_srpaths(solution_file(opts, name), inst, paths)
    end

    return (; instance = name, n_vertices = num_vertices(inst), n_links = length(inst.arcs),
            n_demands = length(inst.demands), n_demands_model = length(m.D),
            max_seg = inst.maxSeg, n_binary_vars = length(m.x),
            n_constraints = num_constraints(m.model; count_variable_in_set_constraints = false),
            baseline_mlu = baseline, mlu_check = check, dropped_segments = dropped,
            data_time_s = data_time, build_time_s = build_time,
            formulation = opts[:formulation], run_info(opts, res)...)
end

function run_instance_t01(dir, name, opts; silent = false)
    data_time = @elapsed begin
        inst = load_instance_t0(dir, name)
        p1 = load_period1(inst, dir, name)
    end
    baseline = max(mlu(inst, direct_paths(inst)), mlu(inst, p1, direct_paths(inst)))

    build_time = @elapsed m = build_model_t01(inst, p1; warm_start = opts[:warm_start])
    prepare!(m.model, opts, silent)
    println("[$name] data $(round(data_time; digits = 2)) s, build $(round(build_time; digits = 2)) s, ",
            "$(length(m.x0) + length(m.x1)) binary variables"); flush(stdout)

    res = solve_t01!(m; time_limit = opts[:time_limit], seed = opts[:seed])

    mlu0, mlu1, used = NaN, NaN, -1
    if res.status in ("Optimal", "Not Optimal")
        paths0, paths1 = extract_paths(m)
        mlu0, mlu1 = mlu(inst, paths0), mlu(inst, p1, paths1)
        used = routing_distance(paths0, paths1)
        write_srpaths(solution_file(opts, name), inst, [paths0, paths1])
    end

    return (; instance = name, n_vertices = num_vertices(inst), n_links = length(inst.arcs),
            n_links_down = length(p1.down), n_demands = length(inst.demands),
            n_demands_model = length(m.D), max_seg = inst.maxSeg, budget = p1.budget,
            n_binary_vars = length(m.x0) + length(m.x1),
            n_constraints = num_constraints(m.model; count_variable_in_set_constraints = false),
            baseline_mlu = baseline, mlu0_check = mlu0, mlu1_check = mlu1, budget_used = used,
            data_time_s = data_time, build_time_s = build_time, run_info(opts, res)...)
end

function main(args)
    instances, opts = parse_args(args)
    run_instance = opts[:model] == "t0" ? run_instance_t0 : run_instance_t01

    # Compile everything on the toy instance first, so that the timings
    # below do not include Julia's just-in-time compilation.
    toy_opts = merge(opts, Dict(:time_limit => 10.0, :results_dir => mktempdir()))
    run_instance(joinpath(@__DIR__, "test", "toy"), "toy", toy_opts; silent = true)

    for k in instances
        name = "setA-$(lpad(k, 2, '0'))"
        row = run_instance(opts[:data_dir], name, opts)
        append_row(opts[:out], COLUMNS[opts[:model]], row)
        println("[$name] $(row.status)  MLU = $(row.objective)  (no waypoints: $(row.baseline_mlu))  ",
                "gap = $(row.gap)  optimize! = $(round(row.cpu_time_s; digits = 2)) s  ",
                "nodes = $(row.bb_nodes)  seed = $(row.seed)  lp = $(row.mip_lp_solver)")
        flush(stdout)
    end
end

main(ARGS)
