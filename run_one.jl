# Một lần chạy = (instance, step, seed):
#   julia --project=. run_one.jl <setA-XX> <6|9> <seed> <exact|restricted> <time_limit> [K]
#   step 6 : chỉ t = 0, không budget         (readme Step 5–6)
#   step 9 : t ∈ {0,1}, có budget (4)         (readme Step 7–9)
# Kết quả: results/raw/step<6|9>/<inst>_s<seed>.json  (+ log HiGHS, + srpaths cho seed 1)
include(joinpath(@__DIR__, "src", "orange.jl"))

const DATA = joinpath(@__DIR__, "..", "project", "setA")
const RES = joinpath(@__DIR__, "results")

periods_of(step) = step == 6 ? [0] : [0, 1]

function run_case(prefix, step, seed, mode, time_limit, K; log_file = "")
    periods = periods_of(step)
    inst = load_instance(prefix)
    t_data = @elapsed begin
        pd = [period_data(inst, t) for t in periods]
        alives = [p[1] for p in pd]
        sds = [p[2] for p in pd]
    end
    P0 = baseline_paths(inst)
    base = maximum(mlu(arc_loads(inst, sds[k], t, P0), alives[k]) for (k, t) in enumerate(periods))
    # Mô hình exact: pha 1 giải nhanh mô hình rút gọn (K = 20 waypoint ứng viên, ≤ 10% time
    # limit). Nghiệm của nó khả thi cho mô hình đầy đủ -> dùng làm warm start thay cho baseline.
    # Lý do: với mô hình đầy đủ, heuristic của HiGHS hầu như không tìm được nghiệm tốt hơn
    # baseline trong time limit (setA-03 step 9: 1 giờ, UB vẫn = baseline).
    warmP, warmU, warm_cpu = [P0 for _ in periods], base, 0.0
    if mode == :exact
        srw = build_model(inst, sds, alives, periods, candidate_segments(inst, sds, :restricted; K = 20))
        set_warm_start!(srw, inst, sds, base)
        rw = solve!(srw; time_limit = 0.1 * time_limit, seed = seed)   # mô hình rút gọn: simplex
        warm_cpu = rw.cpu_time
        if rw.has_sol
            warmP = decode_paths(srw, inst)
            warmU = verify_solution(inst, sds, alives, periods, warmP).mlu
        end
        srw = nothing; GC.gc()
    end
    t_build = @elapsed begin
        cand = candidate_segments(inst, sds, mode; K = K)
        srm = build_model(inst, sds, alives, periods, cand)
        set_warm_start!(srm, inst, sds, warmU; paths = warmP)
    end
    nvars = num_variables(srm.model)
    ncons = num_constraints(srm.model; count_variable_in_set_constraints = false)
    GC.gc()
    # Mô hình exact: LP gốc giải bằng interior point. Dual simplex mặc định không giải xong LP gốc
    # của setA-03 step 9 sau 1 giờ (752 nghìn vòng lặp); IPM giải xong trong ~36 s.
    # Mô hình rút gọn: giữ dual simplex (mặc định). Nghiệm đỉnh của simplex thường đã nguyên
    # (phần lớn instance Optimal ngay ở nút gốc); IPM cho nghiệm trong miền -> setA-20 step 9
    # không tìm được nghiệm tốt hơn baseline sau 20 phút, trong khi simplex: Optimal sau 62 s.
    lp = get(ENV, "ORANGE_LP", mode == :exact ? "ipm" : "simplex")   # ORANGE_LP: thử nghiệm bổ sung
    res = solve!(srm; time_limit = time_limit, seed = seed, log_file = log_file,
                 lp_solver = lp == "ipm" ? "ipm" : "")

    # Nghiệm tốt nhất: từ solver nếu có, nếu không thì chính warm start
    if res.has_sol
        P = decode_paths(srm, inst)
    else
        P = warmP
    end
    chk = verify_solution(inst, sds, alives, periods, P)
    best = res.has_sol ? res.obj : warmU
    abs(chk.mlu - best) <= 1e-6 * max(1, best) || chk.mlu < best ||
        error("MLU tính lại ($(chk.mlu)) > objective ($best)")
    return inst, periods, P, (;
        instance = inst.name, step, seed, mode = string(mode), K = mode == :restricted ? K : 0,
        n = inst.net.n, arcs = narcs(inst), demands = ndemands(inst),
        demands_in_model = length(srm.dems), maxseg = inst.maxseg,
        budget = step == 9 ? inst.budget[1] : -1,
        nvars, ncons, time_data = t_data, time_build = t_build, time_limit,
        res.status, res.termination, res.has_sol,
        # ±Inf (chưa có cận dưới khi hết giờ) không ghi được ra JSON -> null
        warm_obj = warmU, warm_cpu_time = warm_cpu, lp_solver = lp,
        best_obj = best, best_bound = isfinite(res.bound) ? res.bound : nothing, gap = rel_gap(res.bound, best),
        nodes = res.nodes, cpu_time = res.cpu_time, wall_time = res.wall_time,
        baseline_mlu = base, improvement = 1 - best / base,
        mlu_check = chk.mlu, mlu_per_period = chk.mlu_per_period,
        changes = chk.changes, waypoints = chk.waypoints,
        maxrss_gb = Sys.maxrss() / 2^30)
end

function main(args)
    name, step, seed, mode, tl = args[1], parse(Int, args[2]), parse(Int, args[3]),
                                 Symbol(args[4]), parse(Float64, args[5])
    K = length(args) >= 6 ? parse(Int, args[6]) : 10
    # Khởi động JIT trên toy để thời gian compile không bị tính vào optimize!
    run_case(joinpath(@__DIR__, "toy", "toy"), step, 1, :exact, 10.0, K)
    run_case(joinpath(@__DIR__, "toy", "toy"), step, 1, :restricted, 10.0, 3)

    sub = "step$step" * get(ENV, "ORANGE_TAG", "")      # vd. ORANGE_TAG=_ipm -> results/raw/step6_ipm
    dir = joinpath(RES, "raw", sub); mkpath(dir)
    logdir = joinpath(RES, "highs_logs", sub); mkpath(logdir)
    tag = "$(name)_s$(seed)"
    inst, periods, P, out = run_case(joinpath(DATA, name), step, seed, mode, tl, K;
                                     log_file = joinpath(logdir, tag * ".log"))
    if seed == 1
        soldir = joinpath(RES, "solutions", sub); mkpath(soldir)
        write_srpaths(joinpath(soldir, name * "-srpaths.json"), inst, periods, P)
    end
    tmp = joinpath(dir, tag * ".json.tmp")
    open(io -> JSON3.pretty(io, Dict(pairs(out))), tmp, "w")
    mv(tmp, joinpath(dir, tag * ".json"); force = true)     # ghi nguyên tử: file .json có nghĩa là run đã xong
    @printf("%s step%d seed%d %s: %s obj=%.6f bound=%.6f gap=%.4f nodes=%d cpu=%.1fs base=%.4f rss=%.1fGB\n",
            name, step, seed, out.mode, out.status, out.best_obj, something(out.best_bound, NaN), out.gap,
            out.nodes, out.cpu_time, out.baseline_mlu, out.maxrss_gb)
end

main(ARGS)
