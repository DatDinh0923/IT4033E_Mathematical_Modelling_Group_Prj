# Step 7 — graph and parameters at time period 1.
#
#   G_1 = (V, A \ q(1)): the links of the intervention at t = 1 are removed,
#   r(i,j,a,1) is computed on G_1 exactly as in Step 3,
#   v(d,1) comes from the tm file and the budget κ(1) from the scenario file.
#
# Arcs keep their period-0 numbering a = 1..|A|; an arc a ∈ q(1) simply has
# no coefficient r(i,j,a,1) and no load constraint at t = 1.

struct PeriodT1
    down::Vector{Int}                # q(1): arcs out of service at t = 1
    up::Vector{Int}                  # A \ q(1)
    graph::SimpleDiGraph{Int}        # G_1
    v::Vector{Float64}               # v(d,1)
    budget::Int                      # κ(1)
    dist::Matrix{Float64}            # d_ij(1)
    sigma::Matrix{Float64}           # σ_ij(1)
    R::Vector{Vector{Tuple{Int, Int, Float64}}}           # R[a] = [(i, j, r(i,j,a,1)) : r > 0]
    Rseg::Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}
end

"""
    load_period1(inst, dir, name)

Period-1 parameters of the instance whose period-0 data is `inst`.
"""
function load_period1(inst::InstanceT0, dir::AbstractString, name::AbstractString)
    t = 1
    tm = JSON3.read(read(joinpath(dir, "$name-tm.json"), String))
    scenario = JSON3.read(read(joinpath(dir, "$name-scenario.json"), String))
    Int(tm["num_time_slots"]) > t || error("$name has no time period $t")

    # q(1), given as JSON link ids
    arc_of_link = Dict(id => a for (a, id) in enumerate(inst.link_ids))
    down = sort!(unique(arc_of_link[Int(l)] for it in scenario["interventions"]
                        if Int(it["t"]) == t for l in it["links"]))
    up = setdiff(eachindex(inst.arcs), down)

    # G_1 and r(i,j,a,1)
    n = num_vertices(inst)
    g1 = copy(inst.network.graph)
    for a in down
        rem_edge!(g1, inst.arcs[a]...)
    end
    W1 = sparse(first.(inst.arcs[up]), last.(inst.arcs[up]), inst.omega[up], n, n)
    dist, _ = all_pairs_shortest_data(g1, W1)
    sigma = count_shortest_paths(g1, W1, dist)

    R = [Tuple{Int, Int, Float64}[] for _ in inst.arcs]
    R[up] = split_coefficients(inst.arcs[up], inst.omega[up], dist, sigma)

    v = [Float64(dm["v"][t+1]) for dm in tm["demands"]]
    budget = only(Int(b["value"]) for b in scenario["budget"] if Int(b["t"]) == t)

    return PeriodT1(down, up, g1, v, budget, dist, sigma, R, index_by_segment(R))
end

arc_loads(inst::InstanceT0, p1::PeriodT1, paths) =
    arc_loads(p1.Rseg, p1.v, length(inst.arcs), paths)

"Maximum link utilization at t = 1, over the arcs in service."
mlu(inst::InstanceT0, p1::PeriodT1, paths) =
    maximum(arc_loads(inst, p1, paths)[p1.up] ./ inst.c[p1.up])

"""
    routing_distance(paths0, paths1)

dist(P, P') of the subject (Section 2, Budget): number of segments (i,j)
used by a demand in one routing and not in the other, summed over demands.
"""
function routing_distance(paths0, paths1)
    segments(path) = Set(zip(path[1:end-1], path[2:end]))
    return sum(length(symdiff(segments(p0), segments(p1))) for (p0, p1) in zip(paths0, paths1))
end
