# julia --project=. test/runtests.jl
using Test
using JuMP

include(joinpath(@__DIR__, "..", "src", "OrangeSR.jl"))
using .OrangeSR

const SETA = joinpath(@__DIR__, "..", "..", "setA")

# r(s,t,·) as a Dict JSON-link-id => r, for readable comparisons with the subject
function r_by_link(inst, s_id, t_id)
    s, t = inst.network.id_to_vertex[s_id], inst.network.id_to_vertex[t_id]
    return Dict(inst.link_ids[a] => r for (a, r) in inst.Rseg[(s, t)])
end

@testset "toy instance (subject, Section 4)" begin
    toy = load_instance_t0(joinpath(@__DIR__, "toy"), "toy")
    @test num_vertices(toy) == 7 && length(toy.arcs) == 22 && toy.maxSeg == 4
    @test toy.v == [50.0, 100.0]

    # Demand (v0,v5): shortest paths 0-1-3-5, 0-1-4-3-5, 0-1-4-6-5 (σ = 3).
    # Path counting (calculating_r.pdf) gives 1/3 on a13; the subject's Table 2
    # uses per-node ECMP and gives 1/2. Links: 0=a01 2=a13 3=a14 6=a35 18=a43 9=a46 21=a65
    r05 = r_by_link(toy, 0, 5)
    @test r05[0] ≈ 1
    @test r05[2] ≈ 1 / 3 && r05[3] ≈ 2 / 3
    @test r05[6] ≈ 2 / 3 && r05[18] ≈ 1 / 3
    @test r05[9] ≈ 1 / 3 && r05[21] ≈ 1 / 3
    @test length(r05) == 7

    # Demand (v2,v5): paths 2-4-3-5 and 2-4-6-5; identical to Table 2.
    r25 = r_by_link(toy, 2, 5)
    @test r25 == Dict(4 => 1.0, 18 => 0.5, 6 => 0.5, 9 => 0.5, 21 => 0.5)

    @test check_split_coefficients(7, toy.arcs, toy.R, toy.sigma) < 1e-12

    # λ(a46,0) without waypoints: (1/3·50 + 1/2·100)/200 (subject, ECMP: 0.3125)
    a46 = findfirst(==(9), toy.link_ids)
    @test arc_loads(toy, direct_paths(toy))[a46] / toy.c[a46] ≈ (50 / 3 + 50) / 200

    baseline = mlu(toy, direct_paths(toy))
    m = build_model_t0(toy)
    set_silent(m.model)
    res = solve_t0!(m; time_limit = 60)
    @test res.status == "Optimal"
    @test res.objective <= baseline + 1e-9
    paths, _ = extract_paths(m)
    @test mlu(toy, paths) ≈ res.objective atol = 1e-6
    @test all(length(p) - 1 <= toy.maxSeg for p in paths)
end

@testset "SetA: data and r coefficients" begin
    for k in 1:20
        name = "setA-$(lpad(k, 2, '0'))"
        inst = load_instance_t0(SETA, name)
        n = num_vertices(inst)
        @test all(>(0), inst.sigma[i, j] for i in 1:n, j in 1:n)  # strongly connected
        @test check_split_coefficients(n, inst.arcs, inst.R, inst.sigma) < 1e-9
        @test all(inst.v .>= 0) && inst.maxSeg == 6
        if all(isinteger, inst.omega)  # tolerant count == Graphs.jl path counts
            _, sigma_dijkstra = all_pairs_shortest_data(inst.network.graph, metric_matrix(inst.network))
            @test inst.sigma == sigma_dijkstra
        end
    end
end

@testset "toy instance, period 1 (subject, Section 4.3)" begin
    toy = load_instance_t0(joinpath(@__DIR__, "toy"), "toy")
    p1 = load_period1(toy, joinpath(@__DIR__, "toy"), "toy")
    @test sort(toy.link_ids[p1.down]) == [2, 3]          # a13 and a14 are down
    @test p1.v == [100.0, 50.0] && p1.budget == 20

    # Without a13, a14 the flow of (v0,v5) goes 0-2-4 then splits 4-3-5 / 4-6-5
    # (Figure 3, right). Links: 1=a02 4=a24 18=a43 6=a35 9=a46 21=a65
    s, t = toy.network.id_to_vertex[0], toy.network.id_to_vertex[5]
    r05 = Dict(toy.link_ids[a] => r for (a, r) in p1.Rseg[(s, t)])
    @test r05 == Dict(1 => 1.0, 4 => 1.0, 18 => 0.5, 6 => 0.5, 9 => 0.5, 21 => 0.5)
    @test check_split_coefficients(7, toy.arcs, p1.R, p1.sigma) < 1e-12

    @test routing_distance([[1, 3, 5]], [[1, 5]]) == 3                 # Table 1 of the subject
    @test routing_distance([[1, 2, 3, 6]], [[1, 3, 2, 6]]) == 6
    @test routing_distance([[1, 2, 6], [4, 3, 6]], [[1, 7, 6], [4, 3, 6]]) == 4

    m = build_model_t01(toy, p1)
    set_silent(m.model)
    res = solve_t01!(m; time_limit = 60)
    @test res.status == "Optimal"
    paths0, paths1 = extract_paths(m)
    @test max(mlu(toy, paths0), mlu(toy, p1, paths1)) ≈ res.objective atol = 1e-6
    @test routing_distance(paths0, paths1) <= p1.budget
end

@testset "setA-01, periods 0 and 1" begin
    inst = load_instance_t0(SETA, "setA-01")
    p1 = load_period1(inst, SETA, "setA-01")
    @test length(p1.down) == 1 && p1.budget == 51
    @test all(>(0), p1.sigma)                                 # G_1 strongly connected
    @test check_split_coefficients(num_vertices(inst), inst.arcs, p1.R, p1.sigma) < 1e-9
    @test all(isempty, p1.R[p1.down])

    m = build_model_t01(inst, p1)
    set_silent(m.model)
    res = solve_t01!(m; time_limit = 60)     # not solved to optimality in 600 s either
    @test res.status in ("Optimal", "Not Optimal")
    paths0, paths1 = extract_paths(m)
    @test max(mlu(inst, paths0), mlu(inst, p1, paths1)) <= res.objective + 1e-6
    @test routing_distance(paths0, paths1) <= p1.budget
    # period 0 alone can only do better than both periods together
    @test res.objective >= 0.929383899821109 - 1e-6
end

@testset "setA-01: aggregated (3) == literal (3)" begin
    inst = load_instance_t0(SETA, "setA-01")
    obj = map((:aggregated, :direct)) do formulation
        m = build_model_t0(inst; formulation)
        set_silent(m.model)
        res = solve_t0!(m; time_limit = 600)
        @test res.status == "Optimal"
        paths, _ = extract_paths(m)
        @test mlu(inst, paths) ≈ res.objective atol = 1e-6
        res.objective
    end
    @test obj[1] ≈ obj[2] rtol = 1e-4
end
