# Step 2 — JSON network file → Graphs.jl graph.
#
# Adapted from the course notebook `project/json_to_graph_pluto_v2.jl`.
# Only change: files flagged `"multigraph": true` are accepted. In SetA the
# instances 05, 06, 08, 10, 12, 13, 14, 16, 17 and 19 carry that flag but
# contain no parallel arcs, so a simple digraph represents them exactly.
# A real parallel arc is still rejected by `add_edge!` below.

struct NetworkData{G<:Graphs.AbstractGraph}
    graph::G

    # Node information
    id_to_vertex::Dict{Int, Int}
    vertex_to_id::Vector{Int}
    vertex_name::Vector{String}

    # Link information, indexed by Julia endpoints (u,v)
    link_id::Dict{Tuple{Int, Int}, Int}
    metric::Dict{Tuple{Int, Int}, Float64}
    capacity::Dict{Tuple{Int, Int}, Float64}
end

# Canonical key for an edge.
# For a directed graph:   (u,v)
# For an undirected graph: (min(u,v), max(u,v))
edge_key(directed::Bool, u::Int, v::Int) =
    directed ? (u, v) : (u <= v ? (u, v) : (v, u))

function json_to_graph(filename::AbstractString)
    data = JSON3.read(read(filename, String))

    # ---------- Check top-level JSON structure ----------
    required_top = ("directed", "multigraph", "nodes", "links")
    missing_top = [k for k in required_top if !haskey(data, k)]

    isempty(missing_top) ||
        error("Missing top-level JSON field(s): $(join(missing_top, ", "))")

    directed = Bool(data["directed"])

    nodes = data["nodes"]
    links = data["links"]

    # ---------- Nodes ----------
    isempty(nodes) && error("The JSON file contains no nodes.")

    for (i, node) in enumerate(nodes)
        haskey(node, "id") ||
            error("Node $i has no 'id' field.")
        haskey(node, "name") ||
            error("Node $i has no 'name' field.")
    end

    vertex_to_id = Int[Int(node["id"]) for node in nodes]

    length(unique(vertex_to_id)) == length(vertex_to_id) ||
        error("Node IDs must be unique.")

    id_to_vertex = Dict(
        json_id => julia_vertex
        for (julia_vertex, json_id) in enumerate(vertex_to_id)
    )

    vertex_name = String[String(node["name"]) for node in nodes]

    # ---------- Graph ----------
    n = length(nodes)

    g = if directed
        SimpleDiGraph(n)
    else
        SimpleGraph(n)
    end

    # ---------- Link attributes ----------
    link_id = Dict{Tuple{Int, Int}, Int}()
    metric = Dict{Tuple{Int, Int}, Float64}()
    capacity = Dict{Tuple{Int, Int}, Float64}()

    seen_link_ids = Set{Int}()

    for (i, link) in enumerate(links)
        for field in ("id", "from", "to", "metric", "capacity")
            haskey(link, field) ||
                error("Link $i has no '$field' field.")
        end

        e_id = Int(link["id"])
        from_id = Int(link["from"])
        to_id = Int(link["to"])

        e_id in seen_link_ids &&
            error("Duplicate link ID: $e_id")
        push!(seen_link_ids, e_id)

        haskey(id_to_vertex, from_id) ||
            error("Link $e_id refers to unknown node ID $from_id.")

        haskey(id_to_vertex, to_id) ||
            error("Link $e_id refers to unknown node ID $to_id.")

        u = id_to_vertex[from_id]
        v = id_to_vertex[to_id]

        # Graphs.jl SimpleGraph/SimpleDiGraph does not allow parallel edges.
        add_edge!(g, u, v) ||
            error("Duplicate/parallel link between node IDs $from_id and $to_id.")

        key = edge_key(directed, u, v)

        link_id[key] = e_id
        metric[key] = Float64(link["metric"])
        capacity[key] = Float64(link["capacity"])
    end

    return NetworkData(
        g,
        id_to_vertex,
        vertex_to_id,
        vertex_name,
        link_id,
        metric,
        capacity,
    )
end

function attribute_matrix(network::NetworkData, attribute::Symbol)
    values = if attribute === :metric
        network.metric
    elseif attribute === :capacity
        network.capacity
    else
        error("attribute must be :metric or :capacity")
    end

    directed = is_directed(network.graph)
    n = nv(network.graph)

    I = Int[]
    J = Int[]
    V = Float64[]

    for ((u, v), value) in values
        push!(I, u)
        push!(J, v)
        push!(V, value)

        if !directed && u != v
            push!(I, v)
            push!(J, u)
            push!(V, value)
        end
    end

    return sparse(I, J, V, n, n)
end

metric_matrix(network::NetworkData) =
    attribute_matrix(network, :metric)

capacity_matrix(network::NetworkData) =
    attribute_matrix(network, :capacity)
