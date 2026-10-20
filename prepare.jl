# Phân tích một instance trước khi giải (Step 2–4 + số liệu cho quyết định chạy):
#   julia --project=. prepare.jl setA-01 [setA-02 ...]
# Ghi results/instances/<name>.json
include(joinpath(@__DIR__, "src", "orange.jl"))

const DATA = joinpath(@__DIR__, "..", "project", "setA")
const OUT = joinpath(@__DIR__, "results", "instances")

nan2null(x) = isnan(x) ? nothing : x

function prepare(name)
    inst = load_instance(joinpath(DATA, name))
    n, m, D = inst.net.n, narcs(inst), ndemands(inst)
    t_r = @elapsed begin
        al0, sd0 = period_data(inst, 0)
        al1, sd1 = period_data(inst, 1)
    end
    err0 = check_r_conservation(inst.net, sd0)
    err1 = check_r_conservation(inst.net, sd1)
    # sai số nếu chỉ dùng Graphs.jl (so sánh float tuyệt đối) — để ghi nhận
    g0, W0, _ = graph_at(inst.net)
    err_graphsjl = check_r_conservation(inst.net, compute_r(inst.net, g0, W0, al0))
    npairs0 = count(p -> reachable(sd0, p), 1:n*n)
    npairs1 = count(p -> reachable(sd1, p), 1:n*n)
    fg0 = length(sd0.arc); fg1 = length(sd1.arc)
    P = baseline_paths(inst)
    base0 = mlu(arc_loads(inst, sd0, 0, P), al0)
    base1 = mlu(arc_loads(inst, sd1, 1, P), al1)
    # Không có budget nên bài MCF hai chu kỳ tách rời: min max(U0,U1) = max(LB(t=0), LB(t=1)).
    # Giải từng chu kỳ riêng -> LP nhỏ bằng một nửa (quan trọng với giới hạn RAM 15 GB).
    t_lb = @elapsed begin
        lb0 = mcf_lower_bound(inst, [al0], [0])
        lb1 = mcf_lower_bound(inst, [al1], [1])
        lb01 = max(lb0, lb1)          # NaN nếu một trong hai không giải xong
    end
    Dpos0 = count(>(0), inst.vol[:, 1])
    Dpos01 = count(d -> inst.vol[d, 1] > 0 || inst.vol[d, 2] > 0, 1:D)
    info = Dict(
        "instance" => name, "n" => n, "arcs" => m, "demands" => D,
        "zero_volume_t0" => D - Dpos0, "maxseg" => inst.maxseg,
        "budget_t1" => inst.budget[1], "down_links_t1" => inst.net.link_id[inst.down[1]],
        "time_r_s" => t_r, "r_conservation_err_t0" => err0, "r_conservation_err_t1" => err1,
        "r_err_graphsjl_t0" => err_graphsjl, "tolerant_dijkstra" => err_graphsjl > 1e-6,
        "pairs_t0" => npairs0, "pairs_t1" => npairs1, "sum_fg_t0" => fg0, "sum_fg_t1" => fg1,
        # kích thước mô hình đầy đủ (exact)
        "x_vars_step6" => Dpos0 * npairs0,
        "nnz_load_direct_t0" => Dpos0 * fg0,
        "nnz_load_aggregated_t0" => Dpos0 * npairs0 + fg0,
        "x_vars_step9" => 2 * Dpos01 * max(npairs0, npairs1),
        "baseline_mlu_t0" => base0, "baseline_mlu_t1" => base1,
        "baseline_mlu_t01" => max(base0, base1),
        # NaN (LP không giải xong trong time limit) -> null trong JSON
        "mcf_lb_t0" => nan2null(lb0), "mcf_lb_t1" => nan2null(lb1), "mcf_lb_t01" => nan2null(lb01),
        "time_mcf_s" => t_lb, "maxrss_gb" => Sys.maxrss() / 2^30,
    )
    mkpath(OUT)
    open(io -> JSON3.pretty(io, info), joinpath(OUT, name * ".json"), "w")
    @printf("%-9s n=%3d m=%4d D=%4d |P|=%6d ΣFG=%8d  r-err=%.1e/%.1e (Graphs.jl %.1e)  base=%.4f/%.4f  LB=%.4f/%.4f  t_r=%.1fs t_lb=%.1fs rss=%.1fGB\n",
            name, n, m, D, npairs0, fg0, err0, err1, err_graphsjl, base0, max(base0, base1), lb0, lb01, t_r, t_lb, Sys.maxrss() / 2^30)
    flush(stdout)
end

for name in ARGS
    prepare(name)
end
