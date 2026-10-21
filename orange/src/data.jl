# =============================================================================
#  data.jl — Data section of the Orange ROADEF 2026 project (T-ASR)
#
#  Reads the three JSON files of an instance and computes every model parameter:
#     c(a)          capacity of arc a                         (net file)
#     G(t)          graph at period t = G minus arcs in q(t)  (net + scenario)
#     r[i,j,a,t]    split coefficient of segment (i,j) on a   (computed, Dijkstra)
#     v(d,t)        traffic volume of demand d at period t    (tm file)
#     MaxSeg, κ(t)  max segments, budget                      (scenario file)
#
#  JSON node ids start at 0; Julia vertices are 1..n (mapping kept explicitly).
# =============================================================================

using JSON3, Graphs, SparseArrays

struct Arc
    id::Int          # JSON link id
    u::Int           # Julia tail vertex
    v::Int           # Julia head vertex
    metric::Float64  # ω(a)
    cap::Float64     # c(a)
end

struct Demand
    s::Int                    # Julia source vertex
    t::Int                    # Julia target vertex
    vol::Vector{Float64}      # vol[t+1] = v(d,t)
end

"""
Split coefficients for one period, stored sparsely:
  seg[(i,j)] = [(arc index k, r[i,j,a_k,t]), ...]   (only nonzero entries)
  reach[i,j] = true  if j is reachable from i in G(t) (segment (i,j) usable)
"""
struct PeriodData
    t::Int
    active::Vector{Bool}                                  # arc k available at t
    graph::SimpleDiGraph{Int}
    dist::Matrix{Float64}                                 # d_ij(t)
    sigma::Matrix{Float64}                                # σ_ij(t)
    reach::BitMatrix
    seg::Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}
end

struct Instance
    name::String
    n::Int
    names::Vector{String}
    id_to_vertex::Dict{Int,Int}
    vertex_to_id::Vector{Int}
    arcs::Vector{Arc}
    arc_index::Dict{Int,Int}          # JSON link id -> position k in `arcs`
    demands::Vector{Demand}
    T::Int                            # number of periods h
    maxseg::Int
    budget::Dict{Int,Int}             # κ(t), t ≥ 1
    down::Dict{Int,Vector{Int}}       # q(t): arc indices down at t ≥ 1
    periods::Vector{PeriodData}       # periods[t+1]
end

# -----------------------------------------------------------------------------
# Reading the three JSON files
# -----------------------------------------------------------------------------
function read_instance(dir::AbstractString, name::AbstractString;
                       periods = nothing, rmethod::Symbol = :pathcount)
    net = JSON3.read(read(joinpath(dir, "$name-net.json"), String))
    tm  = JSON3.read(read(joinpath(dir, "$name-tm.json"), String))
    sc  = JSON3.read(read(joinpath(dir, "$name-scenario.json"), String))

    Bool(net["directed"]) || error("network must be directed")

    # ---- nodes ----
    vertex_to_id = Int[Int(nd["id"]) for nd in net["nodes"]]
    names        = String[String(nd["name"]) for nd in net["nodes"]]
    id_to_vertex = Dict(id => k for (k, id) in enumerate(vertex_to_id))
    n = length(vertex_to_id)

    # ---- arcs: c(a), ω(a) ----
    # Some SetA files declare multigraph=true but contain no parallel arcs;
    # we check this explicitly, so one simple-digraph code path handles all.
    arcs = Arc[]
    seen = Set{Tuple{Int,Int}}()
    for l in net["links"]
        u = id_to_vertex[Int(l["from"])]; v = id_to_vertex[Int(l["to"])]
        (u, v) in seen && error("parallel arcs $u→$v not supported")
        push!(seen, (u, v))
        push!(arcs, Arc(Int(l["id"]), u, v, Float64(l["metric"]), Float64(l["capacity"])))
    end
    arc_index = Dict(a.id => k for (k, a) in enumerate(arcs))

    # ---- demands: v(d,t) ----
    h = Int(tm["num_time_slots"])
    demands = [Demand(id_to_vertex[Int(d["s"])], id_to_vertex[Int(d["t"])],
                      Float64.(collect(d["v"]))) for d in tm["demands"]]

    # ---- scenario: MaxSeg, κ(t), q(t) ----
    maxseg = Int(sc["max_segments"])
    budget = Dict(Int(b["t"]) => Int(b["value"]) for b in sc["budget"])
    down   = Dict{Int,Vector{Int}}()
    for iv in sc["interventions"]
        append!(get!(down, Int(iv["t"]), Int[]), [arc_index[Int(x)] for x in iv["links"]])
    end

    # ---- G(t) and r[i,j,a,t] for each requested period ----
    tlist = periods === nothing ? (0:h-1) : periods
    pdata = PeriodData[]
    for t in tlist
        active = trues(length(arcs))
        t > 0 && for k in get(down, t, Int[]); active[k] = false; end
        push!(pdata, period_data(n, arcs, Vector{Bool}(active), t; method = rmethod))
    end

    return Instance(name, n, names, id_to_vertex, vertex_to_id, arcs, arc_index,
                    demands, h, maxseg, budget, down, pdata)
end

# -----------------------------------------------------------------------------
# Graph G(t) and its arc-length (metric) matrix
# -----------------------------------------------------------------------------
function build_graph(n, arcs, active)
    g = SimpleDiGraph(n)
    I = Int[]; J = Int[]; W = Float64[]
    for (k, a) in enumerate(arcs)
        active[k] || continue
        add_edge!(g, a.u, a.v)
        push!(I, a.u); push!(J, a.v); push!(W, a.metric)
    end
    return g, sparse(I, J, W, n, n)
end

# all-pairs distances and shortest-path counts (calculating_r.pdf, §1)
function all_pairs_shortest_data(g, distmx)
    n = nv(g)
    d = fill(Inf, n, n)
    sigma = zeros(Float64, n, n)
    for i in vertices(g)
        st = dijkstra_shortest_paths(g, i, distmx; allpaths = true)
        d[i, :] .= st.dists
        sigma[i, :] .= st.pathcounts
    end
    return d, sigma
end

# tolerant equality: metrics may be non-integer floats
@inline sp_equal(x, y) = abs(x - y) <= 1e-9 * max(1.0, abs(y))

"""
    period_data(n, arcs, active, t; method)

method = :pathcount  -> r = σ_iu σ_vj / σ_ij           (calculating_r.pdf, course)
method = :ecmp       -> r = ECMP even split at nodes   (official challenge, Fig. 1)
"""
function period_data(n, arcs, active, t; method::Symbol = :pathcount)
    g, W = build_graph(n, arcs, active)
    d, sigma = all_pairs_shortest_data(g, W)
    reach = BitMatrix(isfinite.(d))
    seg = Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}()
    act = findall(active)

    if method === :pathcount
        # calculating_r.pdf §2–3: arc a=(u,v) is on a shortest i–j path iff
        # d_iu + c_a + d_vj = d_ij ; then r = σ_iu σ_vj / σ_ij
        for i in 1:n, j in 1:n
            (i == j || !reach[i, j]) && continue
            lst = Tuple{Int,Float64}[]
            dij = d[i, j]
            for k in act
                a = arcs[k]
                if sp_equal(d[i, a.u] + a.metric + d[a.v, j], dij)
                    push!(lst, (k, sigma[i, a.u] * sigma[a.v, j] / sigma[i, j]))
                end
            end
            seg[(i, j)] = lst
        end
    elseif method === :ecmp
        # Official ECMP rule: flow entering a node of FG(i,j) is split evenly
        # among its outgoing FG arcs. Process nodes by increasing d_i·.
        out = [Int[] for _ in 1:n]
        for k in act; push!(out[arcs[k].u], k); end
        flow = zeros(n)
        for i in 1:n
            order = sortperm(view(d, i, :))
            for j in 1:n
                (i == j || !reach[i, j]) && continue
                fill!(flow, 0.0); flow[i] = 1.0
                acc = Dict{Int,Float64}()
                for u in order
                    isfinite(d[i, u]) || break
                    flow[u] == 0.0 && continue
                    u == j && continue
                    nxt = [k for k in out[u] if
                           sp_equal(d[i, u] + arcs[k].metric + d[arcs[k].v, j], d[i, j]) &&
                           sp_equal(d[i, u] + arcs[k].metric, d[i, arcs[k].v])]
                    isempty(nxt) && continue
                    share = flow[u] / length(nxt)
                    for k in nxt
                        acc[k] = get(acc, k, 0.0) + share
                        flow[arcs[k].v] += share
                    end
                end
                seg[(i, j)] = [(k, r) for (k, r) in acc]
            end
        end
    else
        error("unknown method $method")
    end
    return PeriodData(t, active, g, d, sigma, reach, seg)
end

# convenience accessor r(i,j,a,t) by JSON link id (for checks / printing)
function r_value(inst::Instance, t::Int, i::Int, j::Int, link_id::Int)
    P = inst.periods[t+1]
    k = inst.arc_index[link_id]
    for (kk, r) in get(P.seg, (i, j), Tuple{Int,Float64}[])
        kk == k && return r
    end
    return 0.0
end