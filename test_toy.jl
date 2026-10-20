# Kiểm tra toàn bộ pipeline trên instance 'toy' của đề (mục 4–5 của Subject).
#   julia --project=. test_toy.jl
include(joinpath(@__DIR__, "src", "orange.jl"))
using Test

const TOY = joinpath(@__DIR__, "toy", "toy")

# ECMP chia đều tại từng nút (Fig. 1 / Table 2 của đề) — CHỈ dùng để đối chiếu
# hàm tính tải với các con số chính thức của đề; mô hình dùng r của calculating_r.pdf.
function ecmp_split(net::Network, sd::SplitData)
    n = net.n
    r = similar(sd.r)
    for i in 1:n, j in 1:n
        reachable(sd, i, j) || continue
        rg = fg_range(sd, pid(n, i, j))
        outs = Dict{Int,Vector{Int}}()
        for e in rg
            push!(get!(outs, net.arc_from[sd.arc[e]], Int[]), e)
        end
        inflow = zeros(n); inflow[i] = 1.0
        for u in sortperm(sd.dist[i, :])          # đi theo thứ tự khoảng cách từ i
            es = get(outs, u, Int[])
            for e in es
                r[e] = inflow[u] / length(es)
                inflow[net.arc_to[sd.arc[e]]] += r[e]
            end
        end
    end
    return SplitData(sd.n, sd.dist, sd.sigma, sd.ptr, sd.arc, r)
end

rval(inst, sd, i, j, from, to) = begin
    n = inst.net.n
    for e in fg_range(sd, pid(n, i + 1, j + 1))
        a = sd.arc[e]
        (inst.net.arc_from[a], inst.net.arc_to[a]) == (from + 1, to + 1) && return sd.r[e]
    end
    0.0
end

@testset "toy" begin
    inst = load_instance(TOY)
    @test inst.net.n == 7 && narcs(inst) == 22 && ndemands(inst) == 2
    @test inst.maxseg == 4 && inst.budget[1] == 20
    @test sort(inst.net.link_id[inst.down[1]]) == [2, 3]

    al0, sd0 = period_data(inst, 0)
    al1, sd1 = period_data(inst, 1)
    @test check_r_conservation(inst.net, sd0) < 1e-12
    @test check_r_conservation(inst.net, sd1) < 1e-12

    # r theo calculating_r.pdf (đếm số đường ngắn nhất): demand (v0,v5) có 3 đường
    @test rval(inst, sd0, 0, 5, 0, 1) ≈ 1
    @test rval(inst, sd0, 0, 5, 1, 3) ≈ 1 / 3
    @test rval(inst, sd0, 0, 5, 1, 4) ≈ 2 / 3
    @test rval(inst, sd0, 0, 5, 3, 5) ≈ 2 / 3
    # demand (v2,v5): trùng với Table 2 của đề
    for (a, b, val) in [(2, 4, 1.0), (3, 5, 0.5), (4, 3, 0.5), (4, 6, 0.5), (6, 5, 0.5)]
        @test rval(inst, sd0, 2, 5, a, b) ≈ val
    end

    # Hàm tính tải: với r kiểu ECMP phải ra đúng Table 2 và Fig. 4 của đề
    e0, e1 = ecmp_split(inst.net, sd0), ecmp_split(inst.net, sd1)
    @test rval(inst, e0, 0, 5, 1, 3) ≈ 0.5 && rval(inst, e0, 0, 5, 3, 5) ≈ 0.75
    P = baseline_paths(inst)
    λ0 = arc_loads(inst, e0, 0, P)
    @test λ0[findfirst(k -> (inst.net.arc_from[k], inst.net.arc_to[k]) == (5, 7), 1:22)] ≈ 0.3125
    @test max(mlu(λ0, al0), mlu(arc_loads(inst, e1, 1, P), al1)) ≈ 0.4375
    Pw = copy(P); Pw[1] = [1, 5, 6]                       # ⟨v0, v4, v5⟩ tại t = 0
    @test max(mlu(arc_loads(inst, e0, 0, Pw), al0), mlu(arc_loads(inst, e1, 1, P), al1)) ≈ 0.375

    # Mô hình t = 0 (Step 5) và t ∈ {0,1} (Step 8), cả hai cách viết ràng buộc tải
    for periods in ([0], [0, 1]), form in (:aggregated, :direct), bf in (:one_row, :two_rows)
        (periods == [0] && bf == :two_rows) && continue
        sds = periods == [0] ? [sd0] : [sd0, sd1]
        als = periods == [0] ? [al0] : [al0, al1]
        base = maximum(mlu(arc_loads(inst, sds[k], t, P), als[k]) for (k, t) in enumerate(periods))
        cand = candidate_segments(inst, sds, :exact)
        srm = build_model(inst, sds, als, periods, cand; load_form = form, budget_form = bf)
        set_warm_start!(srm, inst, sds, base)
        res = solve!(srm; time_limit = 60)
        @test res.status == "Optimal"
        @test res.obj <= base + 1e-9
        Psol = decode_paths(srm, inst)
        chk = verify_solution(inst, sds, als, periods, Psol)
        @test chk.mlu ≈ res.obj atol = 1e-6
        println("periods=$periods form=$form budget=$bf  baseline=$(round(base, digits=4))  opt=$(round(res.obj, digits=4))  ",
                "paths=$(Psol)  changes=$(chk.changes)")
    end

    # Budget κ(1) = 0: routing ở t = 1 phải trùng t = 0
    inst0 = Instance(inst.name, inst.net, inst.src, inst.dst, inst.vol, inst.maxseg, Dict(1 => 0), inst.down)
    srm = build_model(inst0, [sd0, sd1], [al0, al1], [0, 1], candidate_segments(inst0, [sd0, sd1], :exact))
    res = solve!(srm; time_limit = 60)
    Psol = decode_paths(srm, inst0)
    @test res.status == "Optimal" && Psol[1] == Psol[2]
    @test verify_solution(inst0, [sd0, sd1], [al0, al1], [0, 1], Psol).changes == [0]

    # Cận dưới MCF ≤ tối ưu
    lb = mcf_lower_bound(inst, [al0], [0])
    @test lb <= 0.375 + 1e-9
    println("MCF lower bound t=0: ", lb)
end
