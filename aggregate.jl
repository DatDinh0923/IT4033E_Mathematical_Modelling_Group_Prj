# Tổng hợp các lần chạy -> bảng kết quả (trung bình trên các seed):
#   julia --project=. aggregate.jl
# Ghi results/table_step6.{csv,md}, results/table_step9.{csv,md}, results/all_runs_step*.csv
using JSON3, Printf, Statistics

const RES = joinpath(@__DIR__, "results")

fmt(x; d = 4) = isnan(x) ? "–" : @sprintf("%.*f", d, x)
mstd(v; d = 1) = length(v) <= 1 ? fmt(mean(v); d) : fmt(mean(v); d) * " ± " * fmt(std(v); d)

function aggregate(step)
    dir = joinpath(RES, "raw", "step$step")
    isdir(dir) || return
    runs = [JSON3.read(read(joinpath(dir, f), String)) for f in readdir(dir) if endswith(f, ".json")]
    isempty(runs) && return
    cols = ["instance", "seed", "mode", "K", "n", "arcs", "demands", "nvars", "ncons", "status",
            "best_obj", "best_bound", "gap", "nodes", "cpu_time", "wall_time", "time_build",
            "baseline_mlu", "improvement", "waypoints", "changes", "maxrss_gb"]
    open(joinpath(RES, "all_runs_step$step.csv"), "w") do io
        println(io, join(cols, ","))
        for r in sort(runs; by = r -> (r["instance"], r["seed"]))
            println(io, join([c == "changes" ? join(r[c], ";") : string(r[c]) for c in cols], ","))
        end
    end

    insts = sort(unique(r["instance"] for r in runs))
    info(nm) = JSON3.read(read(joinpath(RES, "instances", nm * ".json"), String))
    rows = []
    for nm in insts
        rs = [r for r in runs if r["instance"] == nm]
        I = info(nm)
        lb = something(step == 6 ? I["mcf_lb_t0"] : I["mcf_lb_t01"], NaN)
        objs = [r["best_obj"] for r in rs]
        st = Dict{String,Int}()
        foreach(r -> st[r["status"]] = get(st, r["status"], 0) + 1, rs)
        push!(rows, (; instance = nm, n = I["n"], arcs = I["arcs"], demands = I["demands"],
            mode = rs[1]["mode"] == "exact" ? "exact" : "restricted(K=$(rs[1]["K"]))",
            runs = length(rs), nvars = rs[1]["nvars"], ncons = rs[1]["ncons"],
            status = join(["$k×$v" for (k, v) in sort(collect(st))], ", "),
            obj_mean = mean(objs), obj_best = minimum(objs),
            bound_mean = mean(isnothing(r["best_bound"]) ? 0.0 : max(r["best_bound"], 0.0) for r in rs),
            gap_mean = mean(r["gap"] for r in rs),
            cpu = [r["cpu_time"] for r in rs], nodes = [Float64(r["nodes"]) for r in rs],
            baseline = rs[1]["baseline_mlu"], mcf_lb = lb,
            gap_mcf = 1 - lb / minimum(objs),
            improvement = mean(r["improvement"] for r in rs),
            changes = step == 9 ? mean(r["changes"][1] for r in rs) : NaN,
            budget = step == 9 ? rs[1]["budget"] : -1,
            rss = maximum(r["maxrss_gb"] for r in rs),
            warm = rs[1]["mode"] == "exact" ? mean(r["warm_obj"] for r in rs) : NaN))
    end

    open(joinpath(RES, "table_step$step.csv"), "w") do io
        println(io, "instance,n,arcs,demands,mode,runs,nvars,ncons,status,best_obj_mean,best_obj_min,best_bound_mean,gap_mean,cpu_mean,cpu_std,nodes_mean,nodes_std,baseline_mlu,improvement_mean,mcf_lb,gap_vs_mcf_lb,changes_mean,budget,max_rss_gb")
        for r in rows
            println(io, join(Any[r.instance, r.n, r.arcs, r.demands, r.mode, r.runs, r.nvars, r.ncons,
                "\"$(r.status)\"", r.obj_mean, r.obj_best, r.bound_mean, r.gap_mean, mean(r.cpu),
                length(r.cpu) > 1 ? std(r.cpu) : 0.0, mean(r.nodes), length(r.nodes) > 1 ? std(r.nodes) : 0.0,
                r.baseline, r.improvement, r.mcf_lb, r.gap_mcf, r.changes, r.budget, r.rss], ","))
        end
    end

    title = step == 6 ? "Step 6 — chỉ t = 0, time limit $(Int(runs[1]["time_limit"])) s" :
                        "Step 9 — t ∈ {0,1} + budget, time limit $(Int(runs[1]["time_limit"])) s"
    open(joinpath(RES, "table_step$step.md"), "w") do io
        println(io, "## $title\n")
        println(io, "Giá trị là trung bình trên `runs` lần chạy (seed HiGHS khác nhau); CPU time = thời gian của `optimize!`.\n")
        print(io, "| Instance | \\|V\\| | \\|A\\| | \\|D\\| | Mô hình | runs | Status | Best obj (TB) | Best bound (TB) | Gap (TB) | CPU time (s) | Nodes | MLU baseline | Warm start | Giảm MLU | Cận dưới MCF | Gap theo MCF |")
        step == 9 && print(io, " Thay đổi / κ |")
        println(io)
        println(io, "|" * repeat("---|", step == 9 ? 18 : 17))
        for r in rows
            print(io, "| $(r.instance) | $(r.n) | $(r.arcs) | $(r.demands) | $(r.mode) | $(r.runs) | $(r.status) | ",
                  "$(fmt(r.obj_mean)) | $(fmt(r.bound_mean)) | $(fmt(100r.gap_mean; d=2))% | $(mstd(r.cpu)) | ",
                  "$(mstd(r.nodes; d=0)) | $(fmt(r.baseline)) | $(fmt(r.warm)) | $(fmt(100r.improvement; d=1))% | $(fmt(r.mcf_lb)) | $(fmt(100r.gap_mcf; d=1))% |")
            step == 9 && print(io, " $(fmt(r.changes; d=1)) / $(r.budget) |")
            println(io)
        end
    end
    println("step $step: $(length(runs)) runs, $(length(rows)) instances -> results/table_step$step.md")
end

aggregate(6)
aggregate(9)
