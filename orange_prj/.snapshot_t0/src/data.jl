# Step 4 — Data section: model parameters for time period 0 from the 3 JSON files.
#
#   net file      → graph G (Step 2), c(a), and r(i,j,a,0) (Step 3)
#   tm file       → demands d = (s_d, t_d) and volumes v(d,0)
#   scenario file → maxSeg
#
# Vertices are the Julia numbers 1:n; JSON ids are kept for output.

struct InstanceT0
    name::String
    network::NetworkData

    arcs::Vector{Tuple{Int, Int}}    # arc a = (u,v), ordered by JSON link id
    link_ids::Vector{Int}            # JSON id of arc a
    omega::Vector{Float64}           # metric ω(a) (arc length for shortest paths)
    c::Vector{Float64}               # capacity c(a)

    demands::Vector{Tuple{Int, Int}} # demand d = (s_d, t_d); d = JSON demand id + 1
    v::Vector{Float64}               # traffic volume v(d,0)
    maxSeg::Int

    dist::Matrix{Float64}            # d_ij(0)
    sigma::Matrix{Float64}           # σ_ij(0)
    R::Vector{Vector{Tuple{Int, Int, Float64}}}           # R[a] = [(i, j, r(i,j,a,0)) : r > 0]
    Rseg::Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}} # (i,j) => [(a, r(i,j,a,0)) : r > 0]
end

num_vertices(inst::InstanceT0) = nv(inst.network.graph)

"""
    load_instance_t0(dir, name)

Reads `dir/name-net.json`, `dir/name-tm.json` and `dir/name-scenario.json`
and computes every parameter of the period-0 model.
"""
function load_instance_t0(dir::AbstractString, name::AbstractString)
    network = json_to_graph(joinpath(dir, "$name-net.json"))
    tm = JSON3.read(read(joinpath(dir, "$name-tm.json"), String))
    scenario = JSON3.read(read(joinpath(dir, "$name-scenario.json"), String))

    # Arcs (no intervention at t = 0: G_0 = G)
    arcs = sort!(collect(keys(network.link_id)); by = e -> network.link_id[e])
    link_ids = [network.link_id[e] for e in arcs]
    omega = [network.metric[e] for e in arcs]
    c = [network.capacity[e] for e in arcs]

    # r(i,j,a,0)
    W = metric_matrix(network)
    dist, _ = all_pairs_shortest_data(network.graph, W)
    sigma = count_shortest_paths(network.graph, W, dist)
    R = split_coefficients(arcs, omega, dist, sigma)

    Rseg = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}()
    for (a, entries) in enumerate(R), (i, j, r) in entries
        push!(get!(Rseg, (i, j), Tuple{Int, Float64}[]), (a, r))
    end

    # Demands and v(d,0): time period 0 is the first entry of "v"
    demands = [(network.id_to_vertex[Int(dm["s"])], network.id_to_vertex[Int(dm["t"])])
               for dm in tm["demands"]]
    v = [Float64(dm["v"][1]) for dm in tm["demands"]]

    maxSeg = Int(scenario["max_segments"])

    return InstanceT0(name, network, arcs, link_ids, omega, c,
                      demands, v, maxSeg, dist, sigma, R, Rseg)
end

"""
    arc_loads(inst, paths)

Load of every arc when demand d follows the segment path `paths[d]`
(a node sequence s, w_1, …, t), using the coefficients r(i,j,a,0):
load(a) = Σ_d Σ_{(i,j) ∈ p_d} r(i,j,a,0) v(d,0). Independent of the MILP,
so it is used both for the no-waypoint baseline and to check solutions.
"""
function arc_loads(inst::InstanceT0, paths)
    load = zeros(length(inst.arcs))
    for (d, path) in enumerate(paths), k in 1:length(path)-1
        for (a, r) in get(inst.Rseg, (path[k], path[k+1]), Tuple{Int, Float64}[])
            load[a] += r * inst.v[d]
        end
    end
    return load
end

"Maximum link utilization max_a load(a) / c(a)."
mlu(inst::InstanceT0, paths) = maximum(arc_loads(inst, paths) ./ inst.c)

"Routing without waypoints: every demand follows the shortest paths s → t."
direct_paths(inst::InstanceT0) = [[s, t] for (s, t) in inst.demands]
