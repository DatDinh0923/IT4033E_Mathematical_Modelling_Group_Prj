# =============================================================================
#  orange.jl — Project Orange (ROADEF/EURO 2026, T-Adaptive Segment Routing)
#
#  Toàn bộ logic của project nằm ở file này:
#    Step 2  : json_to_graph / graph_at        (net.json  -> G(t))
#    Step 3  : all_pairs_shortest_data / compute_r   (r[i,j,a,t], calculating_r.pdf)
#    Step 4  : load_instance                    (Data section: r, c, v, maxSeg, budget)
#    Step 5,8: build_model                      (constraints (1)(2)(3) [+ (4)])
#    Step 6,9: solve!, decode_paths, verify_solution
#    Phụ     : baseline_paths (không waypoint), mcf_lower_bound (cận dưới hợp lệ)
# =============================================================================

using JSON3, Graphs, JuMP, HiGHS, Printf

const RTOL = 1e-9          # tolerance khi so sánh độ dài đường đi (metric không nguyên)

# ----------------------------------------------------------------------------
#  Step 2 — đọc net.json thành đồ thị (chuyển từ notebook json_to_graph_pluto_v2.jl)
# ----------------------------------------------------------------------------

struct Network
    n::Int
    vertex_name::Vector{String}
    vertex_to_id::Vector{Int}         # Julia vertex k  -> JSON node id
    id_to_vertex::Dict{Int,Int}       # JSON node id    -> Julia vertex 1:n
    arc_from::Vector{Int}             # arc k = (arc_from[k], arc_to[k]) theo Julia vertex
    arc_to::Vector{Int}
    metric::Vector{Float64}           # ω(a): trọng số dùng cho shortest path
    capacity::Vector{Float64}         # c(a)
    link_id::Vector{Int}              # arc k -> JSON link id
    id_to_arc::Dict{Int,Int}          # JSON link id -> arc k
end

function json_to_graph(filename::AbstractString)
    data = JSON3.read(read(filename, String))
    for k in ("directed", "multigraph", "nodes", "links")
        haskey(data, k) || error("$filename: thiếu trường '$k'")
    end
    Bool(data["directed"]) || error("$filename: đồ thị phải có hướng")
    # Notebook gốc báo lỗi khi multigraph=true. Trong setA cờ này bật ở vài instance
    # nhưng KHÔNG có cung song song nào -> ta bỏ lỗi đó; graph_at() vẫn kiểm tra
    # thật sự (add_edge! trả về false nếu cung trùng).

    nodes = data["nodes"]
    vertex_to_id = [Int(nd["id"]) for nd in nodes]
    vertex_name = [String(get(nd, "name", string(nd["id"]))) for nd in nodes]
    allunique(vertex_to_id) || error("$filename: node id trùng nhau")
    id_to_vertex = Dict(id => k for (k, id) in enumerate(vertex_to_id))

    links = data["links"]
    m = length(links)
    afrom = Vector{Int}(undef, m); ato = similar(afrom); lid = similar(afrom)
    metric = Vector{Float64}(undef, m); cap = similar(metric)
    for (k, l) in enumerate(links)
        afrom[k] = id_to_vertex[Int(l["from"])]
        ato[k] = id_to_vertex[Int(l["to"])]
        lid[k] = Int(l["id"])
        metric[k] = Float64(l["metric"])
        cap[k] = Float64(l["capacity"])
        afrom[k] == ato[k] && error("$filename: self-loop ở link $(lid[k])")
        metric[k] > 0 || error("$filename: metric phải > 0 (link $(lid[k]))")
        cap[k] > 0 || error("$filename: capacity phải > 0 (link $(lid[k]))")
    end
    allunique(lid) || error("$filename: link id trùng nhau")
    return Network(length(nodes), vertex_name, vertex_to_id, id_to_vertex,
                   afrom, ato, metric, cap, lid, Dict(id => k for (k, id) in enumerate(lid)))
end

"""
    graph_at(net, down_arcs) -> (g, W, alive)

G(t) = (V, A \\ q(t)): đồ thị Graphs.jl chỉ gồm các cung còn sống,
W = ma trận độ dài (metric), alive[k] = cung k còn hoạt động.
"""
function graph_at(net::Network, down_arcs::AbstractVector{Int} = Int[])
    alive = trues(length(net.arc_from))
    alive[down_arcs] .= false
    g = SimpleDiGraph(net.n)
    W = zeros(net.n, net.n)
    for k in eachindex(alive)
        alive[k] || continue
        u, v = net.arc_from[k], net.arc_to[k]
        add_edge!(g, u, v) || error("cung song song ($u,$v): cần multigraph")
        W[u, v] = net.metric[k]
    end
    return g, W, alive
end

# ----------------------------------------------------------------------------
#  Step 3 — r[i,j,a,t] = σ_iu σ_vj / σ_ij nếu d_iu + c_a + d_vj = d_ij  (calculating_r.pdf)
# ----------------------------------------------------------------------------

# Chép nguyên văn từ calculating_r.pdf
function all_pairs_shortest_data(g, distmx)
    n = nv(g)
    d = fill(Inf, n, n)          # d[i,j]     = độ dài đường ngắn nhất i -> j
    sigma = zeros(Float64, n, n) # sigma[i,j] = số đường ngắn nhất i -> j
    for i in vertices(g)
        state = dijkstra_shortest_paths(g, i, distmx; allpaths = true)
        d[i, :] .= state.dists
        sigma[i, :] .= state.pathcounts
    end
    # Graphs.jl dùng typemax(Float64)=Inf cho đỉnh không tới được; pathcounts = 0
    return d, sigma
end

"""
Cùng thuật toán ở mục 1 của calculating_r.pdf (Case 1/2/3), nhưng "ℓ = d(v)" được so
sánh với tolerance. Graphs.jl so sánh `alt == dists[v]` tuyệt đối, nên với metric không
nguyên (setA-08, setA-14) hai đường bằng nhau có thể lệch 1 ulp và bị đếm thiếu.
Dijkstra O(n²) (n ≤ 400 nên đủ nhanh).
"""
function all_pairs_shortest_data_tol(g, distmx; rtol = RTOL)
    n = nv(g)
    d = fill(Inf, n, n)
    sigma = zeros(Float64, n, n)
    visited = falses(n)
    for s in 1:n
        ds = view(d, s, :); σs = view(sigma, s, :)
        fill!(visited, false)
        ds[s] = 0.0; σs[s] = 1.0
        for _ in 1:n
            u = 0; best = Inf
            for w in 1:n
                (!visited[w] && ds[w] < best) && (best = ds[w]; u = w)
            end
            u == 0 && break
            visited[u] = true
            for v in outneighbors(g, u)
                ℓ = ds[u] + distmx[u, v]
                tol = rtol * max(1.0, ℓ)
                if ℓ < ds[v] - tol          # Case 1: đường tốt hơn hẳn
                    ds[v] = ℓ; σs[v] = σs[u]
                elseif abs(ℓ - ds[v]) <= tol # Case 2: thêm đường ngắn nhất
                    σs[v] += σs[u]
                end                          # Case 3: bỏ qua
            end
        end
    end
    return d, sigma
end

"""
Dữ liệu tỉ lệ chia luồng của MỘT chu kỳ t, lưu thưa kiểu CSR theo cặp (i,j):
với p = pid(n,i,j), các cung của forwarding graph FG(i,j) là
`arc[ptr[p]:ptr[p+1]-1]` với hệ số `r[ptr[p]:ptr[p+1]-1]`.
"""
struct SplitData
    n::Int
    dist::Matrix{Float64}
    sigma::Matrix{Float64}
    ptr::Vector{Int}
    arc::Vector{Int32}
    r::Vector{Float64}
end

pid(n, i, j) = (i - 1) * n + j
pair_of(n, p) = (div(p - 1, n) + 1, mod(p - 1, n) + 1)
fg_range(sd::SplitData, p) = sd.ptr[p]:(sd.ptr[p+1]-1)
reachable(sd::SplitData, i, j) = i != j && sd.sigma[i, j] > 0
reachable(sd::SplitData, p) = reachable(sd, pair_of(sd.n, p)...)

function compute_r(net::Network, g, W, alive; tolerant::Bool = false)
    n = net.n
    d, σ = tolerant ? all_pairs_shortest_data_tol(g, W) : all_pairs_shortest_data(g, W)
    arcs = findall(alive)
    ptr = Vector{Int}(undef, n * n + 1)
    arcv = Int32[]; rv = Float64[]
    ptr[1] = 1
    for i in 1:n, j in 1:n                      # thứ tự này đúng với pid = (i-1)n + j
        if i != j && σ[i, j] > 0
            dij = d[i, j]
            tol = RTOL * max(1.0, dij)
            for a in arcs
                u, v = net.arc_from[a], net.arc_to[a]
                dl = d[i, u] + net.metric[a] + d[v, j]   # Inf nếu không tới được
                if abs(dl - dij) <= tol                  # a nằm trên một đường ngắn nhất i->j
                    push!(arcv, a)
                    push!(rv, σ[i, u] * σ[v, j] / σ[i, j])
                end
            end
        end
        ptr[pid(n, i, j)+1] = length(arcv) + 1
    end
    return SplitData(n, d, σ, ptr, arcv, rv)
end

"""
Kiểm tra r: với mỗi cặp (i,j) tới được, r(i,j,·) phải là một luồng đơn vị i -> j
(ra - vào = 1 tại i, -1 tại j, 0 nơi khác). Trả về sai số lớn nhất.
Nếu Dijkstra đếm thiếu đường (do so sánh float), lỗi sẽ lộ ra ở đây.
"""
function check_r_conservation(net::Network, sd::SplitData)
    n = net.n
    bal = zeros(n)
    maxerr = 0.0
    for i in 1:n, j in 1:n
        reachable(sd, i, j) || continue
        fill!(bal, 0.0)
        for e in fg_range(sd, pid(n, i, j))
            a = sd.arc[e]
            bal[net.arc_from[a]] += sd.r[e]
            bal[net.arc_to[a]] -= sd.r[e]
        end
        bal[i] -= 1.0; bal[j] += 1.0
        maxerr = max(maxerr, maximum(abs, bal))
    end
    return maxerr
end

# ----------------------------------------------------------------------------
#  Step 4 — Data section: đọc 3 file JSON của một instance
# ----------------------------------------------------------------------------

struct Instance
    name::String
    net::Network
    src::Vector{Int}                 # s_d (Julia vertex)
    dst::Vector{Int}                 # t_d
    vol::Matrix{Float64}             # vol[d, t+1] = ν(d, t)
    maxseg::Int
    budget::Dict{Int,Int}            # t => κ(t)
    down::Dict{Int,Vector{Int}}      # t => các arc index bị intervention (q(t))
end

ndemands(inst::Instance) = length(inst.src)
narcs(inst::Instance) = length(inst.net.arc_from)

"prefix = đường dẫn không có hậu tố, vd. \".../setA/setA-01\""
function load_instance(prefix::AbstractString)
    net = json_to_graph(prefix * "-net.json")
    tm = JSON3.read(read(prefix * "-tm.json", String))
    sc = JSON3.read(read(prefix * "-scenario.json", String))
    T = Int(tm["num_time_slots"])
    dem = tm["demands"]
    src = [net.id_to_vertex[Int(x["s"])] for x in dem]
    dst = [net.id_to_vertex[Int(x["t"])] for x in dem]
    vol = zeros(length(dem), T)
    for (k, x) in enumerate(dem)
        length(x["v"]) == T || error("demand $k: số giá trị v khác num_time_slots")
        vol[k, :] .= Float64.(x["v"])
    end
    any(src .== dst) && error("demand có s = t")
    budget = Dict(Int(b["t"]) => Int(b["value"]) for b in get(sc, "budget", []))
    down = Dict{Int,Vector{Int}}()
    for it in get(sc, "interventions", [])
        append!(get!(down, Int(it["t"]), Int[]), [net.id_to_arc[Int(l)] for l in it["links"]])
    end
    return Instance(basename(prefix), net, src, dst, vol, Int(sc["max_segments"]), budget, down)
end

"""
Tính (alive, SplitData) cho chu kỳ t; t = 0 không có intervention.
Dùng hàm Graphs.jl của calculating_r.pdf; nếu r không bảo toàn luồng (đếm thiếu
đường do sai số float) thì tính lại bằng phiên bản có tolerance.
"""
function period_data(inst::Instance, t::Int)
    g, W, alive = graph_at(inst.net, get(inst.down, t, Int[]))
    sd = compute_r(inst.net, g, W, alive)
    if check_r_conservation(inst.net, sd) > 1e-6
        sd = compute_r(inst.net, g, W, alive; tolerant = true)
        err = check_r_conservation(inst.net, sd)
        err <= 1e-6 || error("$(inst.name) t=$t: r không bảo toàn luồng (err=$err)")
    end
    return alive, sd
end

# ----------------------------------------------------------------------------
#  Đánh giá một routing scheme (độc lập với mô hình MIP — dùng để kiểm tra)
# ----------------------------------------------------------------------------

"λ(a,t) cho routing `paths[d] = [s, w1, ..., t]`."
function arc_loads(inst::Instance, sd::SplitData, t::Int, paths)
    n = inst.net.n
    flow = zeros(narcs(inst))
    for d in 1:ndemands(inst)
        v = inst.vol[d, t+1]
        v == 0 && continue
        p = paths[d]
        for k in 1:length(p)-1
            for e in fg_range(sd, pid(n, p[k], p[k+1]))
                flow[sd.arc[e]] += sd.r[e] * v
            end
        end
    end
    return flow ./ inst.net.capacity
end

mlu(λ, alive) = maximum(λ[alive])

"Routing không waypoint (chỉ dùng shortest path s -> t): lời giải 'trivial' của đề."
baseline_paths(inst::Instance) = [[inst.src[d], inst.dst[d]] for d in 1:ndemands(inst)]

segments(path) = [(path[k], path[k+1]) for k in 1:length(path)-1]

"dist(P, P') của đề: số segment xuất hiện ở đúng một trong hai routing."
routing_distance(P0, P1) = sum(length(symdiff(Set(segments(a)), Set(segments(b)))) for (a, b) in zip(P0, P1))

# ----------------------------------------------------------------------------
#  Step 5 & 8 — Mô hình MIP
# ----------------------------------------------------------------------------

"""
Tập segment ứng viên cho mỗi demand.

* `:exact` — mọi cặp (i,j), i≠j, tới được ở ít nhất một chu kỳ: đúng mô hình của đề.
  Dùng chung MỘT vector cho mọi demand để tiết kiệm RAM.
* `:restricted` — (chỉ dùng cho instance quá lớn) giữ K waypoint ứng viên cho mỗi
  demand: các đỉnh w có độ vòng d(s,w)+d(w,t)-d(s,t) nhỏ nhất; segment là mọi cặp
  trong {s} ∪ W_d ∪ {t}, bỏ cung đi vào s và đi ra khỏi t (một đường đơn s->t không
  bao giờ dùng chúng). Đây là mô hình con: mọi nghiệm của nó khả thi cho bài gốc,
  nên objective là một cận trên (heuristic), KHÔNG chứng minh được tối ưu.
"""
function candidate_segments(inst::Instance, sds::Vector{SplitData}, mode::Symbol; K::Int = 10)
    n = inst.net.n
    D = ndemands(inst)
    if mode == :exact
        all_p = Int32[pid(n, i, j) for i in 1:n for j in 1:n if any(sd -> reachable(sd, i, j), sds)]
        return fill(all_p, D)          # cùng một object, không copy
    end
    mode == :restricted || error("mode không hợp lệ: $mode")
    d0 = sds[1].dist
    cand = Vector{Vector{Int32}}(undef, D)
    for d in 1:D
        s, t = inst.src[d], inst.dst[d]
        others = [w for w in 1:n if w != s && w != t]
        detour = [d0[s, w] + d0[w, t] - d0[s, t] for w in others]
        W = others[sortperm(detour)[1:min(K, length(others))]]
        N = vcat(s, W, t)
        cand[d] = Int32[pid(n, i, j) for i in N for j in N
                        if i != j && i != t && j != s && any(sd -> reachable(sd, i, j), sds)]
    end
    return cand
end

struct SRModel
    model::Model
    periods::Vector{Int}
    U::VariableRef
    dems::Vector{Int}                          # demand nằm trong mô hình (volume > 0)
    cand::Vector{Vector{Int32}}                # cand[k]: segment ứng viên của dems[k]
    x::Vector{Vector{Vector{VariableRef}}}     # x[ti][k][q] ↔ x^{d,t}_{ij}, (i,j)=cand[k][q]
    y::Vector{Vector{Pair{Int32,VariableRef}}} # y[ti]: (pid => y^t_{ij}) — chỉ có ở load_form=:aggregated
    z::Vector{Tuple{Int,Int,VariableRef}}      # budget: (k, q, z) với z ≥ x^1 - x^0 (chỉ khi có 2 chu kỳ)
    b::Vector{VariableRef}                     # budget: b[k] = số segment thay đổi của dems[k]
end

"""
    build_model(inst, sds, alives, periods, cand) -> SRModel

min U
 (1) Σ_j x^{dt}_{ij} - Σ_j x^{dt}_{ji} = 1 / -1 / 0          ∀d, ∀i, ∀t
 (2) Σ_{ij} x^{dt}_{ij} ≤ maxSeg                                ∀d, ∀t
 (3) Σ_d Σ_{ij} r(i,j,a,t) ν(d,t) x^{dt}_{ij} ≤ U c(a)          ∀t, ∀a ∈ A \\ q(t)
     — cài đặt dạng "tổng hợp theo segment":  y^t_{ij} = Σ_d ν(d,t) x^{dt}_{ij}
       và Σ_{ij} r(i,j,a,t) y^t_{ij} ≤ U c(a). Hai cách viết tương đương tuyệt đối
       (chỉ là đặt tên cho một biểu thức), nhưng số phần tử khác 0 giảm từ
       |D|·Σ|FG(i,j)| xuống |D|·|P| + Σ|FG(i,j)|.  (`load_form = :direct` = viết thẳng)
 (4) Σ_{d,ij} |x^{dt}_{ij} - x^{d,t-1}_{ij}| ≤ κ(t)                                ∀t ∈ T*
     — `budget_form = :two_rows` (cách tuyến tính hóa sách giáo khoa):
           z ≥ Δ, z ≥ -Δ, Σ z ≤ κ                                  (Δ = x^{dt} - x^{d,t-1})
     — `budget_form = :one_row` (mặc định): |Δ| = 2·Δ⁺ - Δ, nên
           z ≥ Δ, z ≥ 0, Σ (2z - Δ) ≤ κ
       Tương đương tuyệt đối (z nhỏ nhất = Δ⁺), nhưng chỉ một ràng buộc / cặp thay vì hai.
     — Tổng được tách theo demand: b_d = Σ_ij |Δ|, Σ_d b_d ≤ κ. Nếu viết một hàng Σ ≤ κ chứa
       MỌI biến z (hàng rất dày), presolve của HiGHS kẹt hàng chục phút và bỏ qua time limit
       (đo trên setA-02, setA-13); tách ra thì presolve còn ~5 s, setA-13 giải tối ưu trong 40 s.
"""
function build_model(inst::Instance, sds::Vector{SplitData}, alives::Vector{BitVector},
                     periods::Vector{Int}, cand::Vector{Vector{Int32}};
                     load_form::Symbol = :aggregated, budget_form::Symbol = :one_row)
    n = inst.net.n
    m = narcs(inst)
    model = direct_model(HiGHS.Optimizer())   # không giữ bản copy cache -> tiết kiệm RAM
    set_silent(model)

    # Demand có volume 0 ở mọi chu kỳ không ảnh hưởng tải: định tuyến thẳng s->t, bỏ khỏi mô hình.
    dems = [d for d in 1:ndemands(inst) if any(inst.vol[d, t+1] > 0 for t in periods)]
    candk = cand[dems]

    @variable(model, U >= 0)
    @objective(model, Min, 1.0 * U)     # AffExpr: direct_model HiGHS không nhận objective dạng biến đơn

    X = Vector{Vector{Vector{VariableRef}}}()
    Y = Vector{Vector{Pair{Int32,VariableRef}}}()
    Z = Tuple{Int,Int,VariableRef}[]
    B = VariableRef[]
    nodeexpr = [AffExpr(0.0) for _ in 1:n]
    for (ti, t) in enumerate(periods)
        sd = sds[ti]
        Xt = Vector{Vector{VariableRef}}(undef, length(dems))
        for (k, d) in enumerate(dems)
            c = candk[k]
            xs = @variable(model, [1:length(c)], Bin)
            for (q, p) in enumerate(c)
                reachable(sd, p) || fix(xs[q], 0.0; force = true)   # segment không tồn tại ở chu kỳ t
            end
            Xt[k] = xs
            # (1) bảo toàn luồng trên "đồ thị segment"
            foreach(e -> (e.constant = 0.0; empty!(e.terms)), nodeexpr)
            for (q, p) in enumerate(c)
                i, j = pair_of(n, p)
                add_to_expression!(nodeexpr[i], 1.0, xs[q])
                add_to_expression!(nodeexpr[j], -1.0, xs[q])
            end
            s, tt = inst.src[d], inst.dst[d]
            for w in 1:n
                rhs = w == s ? 1.0 : (w == tt ? -1.0 : 0.0)
                (isempty(nodeexpr[w].terms) && rhs == 0) && continue
                @constraint(model, nodeexpr[w] == rhs)
            end
            # (2) số segment
            @constraint(model, sum(xs) <= inst.maxseg)
        end
        push!(X, Xt)

        # (3) tải của cung
        load = [AffExpr(0.0) for _ in 1:m]
        Yt = Pair{Int32,VariableRef}[]
        if load_form == :direct
            for (k, d) in enumerate(dems)
                v = inst.vol[d, t+1]; v == 0 && continue
                for (q, p) in enumerate(candk[k]), e in fg_range(sd, p)
                    add_to_expression!(load[sd.arc[e]], sd.r[e] * v, Xt[k][q])
                end
            end
        else
            segexpr = Dict{Int32,AffExpr}()
            for (k, d) in enumerate(dems)
                v = inst.vol[d, t+1]; v == 0 && continue
                for (q, p) in enumerate(candk[k])
                    reachable(sd, p) || continue
                    add_to_expression!(get!(() -> AffExpr(0.0), segexpr, p), v, Xt[k][q])
                end
            end
            for (p, ex) in segexpr
                y = @variable(model, lower_bound = 0.0)          # y^t_{ij}: lưu lượng trên segment (i,j)
                push!(Yt, p => y)
                add_to_expression!(ex, -1.0, y)
                @constraint(model, ex == 0)
                for e in fg_range(sd, p)
                    add_to_expression!(load[sd.arc[e]], sd.r[e], y)
                end
            end
            segexpr = nothing
        end
        push!(Y, Yt)
        for a in 1:m
            alives[ti][a] || continue                           # a ∈ A \ q(t)
            isempty(load[a].terms) && continue                  # cung không bao giờ có luồng
            @constraint(model, (1 / inst.net.capacity[a]) * load[a] <= U)
        end
        load = nothing

        # (4) budget giữa t-1 và t
        if ti > 1
            κ = inst.budget[t]
            total = AffExpr(0.0)
            for k in eachindex(dems)
                bud = AffExpr(0.0)                              # số segment thay đổi của demand k
                for q in eachindex(candk[k])
                    x0, x1 = X[ti-1][k][q], Xt[k][q]
                    (is_fixed(x0) && is_fixed(x1)) && continue  # cả hai bằng 0
                    z = @variable(model, lower_bound = 0.0)
                    push!(Z, (k, q, z))
                    @constraint(model, z >= x1 - x0)
                    if budget_form == :two_rows
                        @constraint(model, z >= x0 - x1)
                        add_to_expression!(bud, 1.0, z)
                    else                                        # |Δ| = 2z - Δ
                        add_to_expression!(bud, 2.0, z)
                        add_to_expression!(bud, -1.0, x1)
                        add_to_expression!(bud, 1.0, x0)
                    end
                end
                # b_d = Σ_ij |Δ|: tách hàng budget dày (mọi biến) thành |D| hàng nhỏ + 1 hàng |D| phần tử
                b = @variable(model, lower_bound = 0.0)
                push!(B, b)
                add_to_expression!(bud, -1.0, b)
                @constraint(model, bud == 0)
                add_to_expression!(total, 1.0, b)
            end
            @constraint(model, total <= κ)
        end
        GC.gc()
    end
    return SRModel(model, periods, U, dems, candk, X, Y, Z, B)
end

"""
Warm start từ một routing `paths[ti][d] = [s, w1, ..., t]` có MLU = `U` (đã kiểm tra khả thi).
Mặc định = routing không waypoint (x^{dt}_{s_d t_d} = 1): luôn khả thi (budget dùng 0), nên
solver luôn có incumbent -> gap luôn xác định. Với mô hình exact, run_one.jl truyền vào
nghiệm của mô hình rút gọn (cũng khả thi cho mô hình đầy đủ).
"""
function set_warm_start!(srm::SRModel, inst::Instance, sds, U;
                         paths = [baseline_paths(inst) for _ in srm.periods])
    n = inst.net.n
    model = srm.model
    xval = [[zeros(length(srm.cand[k])) for k in eachindex(srm.dems)] for _ in srm.periods]
    for ti in eachindex(srm.periods)
        segvol = Dict{Int32,Float64}()          # y^t_{ij} = Σ_d ν(d,t) x^{dt}_{ij}
        for (k, d) in enumerate(srm.dems)
            used = Set(Int32(pid(n, i, j)) for (i, j) in segments(paths[ti][d]))
            for (q, p) in enumerate(srm.cand[k])
                p in used || continue
                xval[ti][k][q] = 1.0
                segvol[p] = get(segvol, p, 0.0) + inst.vol[d, srm.periods[ti]+1]
            end
            sum(xval[ti][k]) == length(used) || error("warm start: segment của demand $d không có trong mô hình")
            for q in eachindex(srm.cand[k])
                set_start_value(srm.x[ti][k][q], xval[ti][k][q])
            end
        end
        for (p, y) in srm.y[ti]
            set_start_value(y, get(segvol, p, 0.0))
        end
    end
    # budget: z = (x¹ - x⁰)⁺,  b_k = Σ (2z - x¹ + x⁰)
    bval = zeros(length(srm.b))
    for (k, q, z) in srm.z
        x0, x1 = xval[1][k][q], xval[2][k][q]
        set_start_value(z, max(x1 - x0, 0.0))
        bval[k] += abs(x1 - x0)
    end
    foreach((b, v) -> set_start_value(b, v), srm.b, bval)
    for v in all_variables(model)
        has_start_value(v) || set_start_value(v, 0.0)
    end
    set_start_value(srm.U, U)
    return
end

# ----------------------------------------------------------------------------
#  Step 6 & 9 — giải, đọc kết quả, decode và kiểm tra nghiệm
# ----------------------------------------------------------------------------

cpu_seconds() = ccall(:clock, Clong, ()) / 1e6   # CPU time của process (CLOCKS_PER_SEC = 1e6)

function solve!(srm::SRModel; time_limit::Real, seed::Int = 1, log_file::String = "",
                lp_solver::String = "")
    model = srm.model
    isempty(lp_solver) || set_attribute(model, "mip_lp_solver", lp_solver)
    set_attribute(model, "time_limit", Float64(time_limit))
    set_attribute(model, "random_seed", seed)
    set_attribute(model, "threads", 1)
    set_attribute(model, "mip_rel_gap", 1e-4)
    if !isempty(log_file)
        set_attribute(model, "log_file", log_file)
        set_attribute(model, "output_flag", true)
        set_attribute(model, "log_to_console", false)
    end
    cpu0 = cpu_seconds()
    wall = @elapsed optimize!(model)
    cpu = cpu_seconds() - cpu0
    ts = termination_status(model)
    has_sol = primal_status(model) == MOI.FEASIBLE_POINT
    obj = has_sol ? objective_value(model) : NaN
    bound = try objective_bound(model) catch; NaN end
    nodes = try node_count(model) catch; -1 end
    status = ts == MOI.OPTIMAL ? "Optimal" :
             ts == MOI.INFEASIBLE ? "Infeasible" :
             has_sol ? "Not Optimal" : "No solution"
    return (; status, termination = string(ts), has_sol, obj, bound, nodes,
              wall_time = wall, cpu_time = cpu)
end

"gap = 1 - LB/UB như readme; LB không hữu hạn được thay bằng 0 (U ≥ 0 luôn đúng)."
function rel_gap(lb, ub)
    lb = isfinite(lb) ? max(lb, 0.0) : 0.0
    ub > 0 || return 0.0
    return max(0.0, 1.0 - lb / ub)
end

"""
Đọc x từ nghiệm và dựng segment path cho từng demand: BFS trên các segment có x = 1
từ s tới t. BFS cho đường đơn, nên nếu nghiệm MIP có chu trình rời rạc (được phép
bởi (1) nhưng chỉ làm tăng tải) thì chu trình đó bị bỏ — tải không tăng.
"""
function decode_paths(srm::SRModel, inst::Instance)
    n = inst.net.n
    P = [baseline_paths(inst) for _ in srm.periods]      # demand volume 0: đi thẳng
    for ti in eachindex(srm.periods), (k, d) in enumerate(srm.dems)
        succ = Dict{Int,Vector{Int}}()
        for (q, p) in enumerate(srm.cand[k])
            value(srm.x[ti][k][q]) > 0.5 || continue
            i, j = pair_of(n, p)
            push!(get!(succ, i, Int[]), j)
        end
        s, t = inst.src[d], inst.dst[d]
        parent = Dict(s => 0); queue = [s]
        while !isempty(queue)
            u = popfirst!(queue)
            u == t && break
            for w in get(succ, u, Int[])
                haskey(parent, w) && continue
                parent[w] = u; push!(queue, w)
            end
        end
        haskey(parent, t) || error("decode: demand $d không có đường s->t ở chu kỳ $(srm.periods[ti])")
        path = [t]
        while path[1] != s
            pushfirst!(path, parent[path[1]])
        end
        P[ti][d] = path
    end
    return P
end

"""
Kiểm tra độc lập một routing: đường đơn s->t, số segment ≤ maxSeg, segment tồn tại ở
G(t), budget, và tính lại MLU bằng `arc_loads` (không dùng giá trị của MIP).
"""
function verify_solution(inst::Instance, sds, alives, periods, P)
    n = inst.net.n
    for (ti, t) in enumerate(periods), d in 1:ndemands(inst)
        p = P[ti][d]
        p[1] == inst.src[d] && p[end] == inst.dst[d] || error("demand $d: sai đầu/cuối")
        allunique(p) || error("demand $d: đường không đơn")
        length(p) - 1 <= inst.maxseg || error("demand $d: vượt maxSeg")
        all(reachable(sds[ti], p[k], p[k+1]) for k in 1:length(p)-1) ||
            error("demand $d: segment không tồn tại ở t=$t")
    end
    mlus = [mlu(arc_loads(inst, sds[ti], t, P[ti]), alives[ti]) for (ti, t) in enumerate(periods)]
    changes = Int[]
    for ti in 2:length(periods)
        dist = routing_distance(P[ti-1], P[ti])
        dist <= inst.budget[periods[ti]] || error("vượt budget ở t=$(periods[ti]): $dist")
        push!(changes, dist)
    end
    nwp = sum(length(p) - 2 for Pt in P for p in Pt)
    return (; mlu = maximum(mlus), mlu_per_period = mlus, changes, waypoints = nwp)
end

"Xuất nghiệm theo định dạng toy-srpaths.json của đề (node id gốc của JSON)."
function write_srpaths(file, inst::Instance, periods, P)
    ids = inst.net.vertex_to_id
    entries = [Dict("d" => d - 1, "t" => t, "w" => ids[P[ti][d][2:end-1]])
               for (ti, t) in enumerate(periods) for d in 1:ndemands(inst) if length(P[ti][d]) > 2]
    open(file, "w") do io
        JSON3.write(io, Dict("srpaths" => entries))
    end
end

# ----------------------------------------------------------------------------
#  Cận dưới hợp lệ: multi-commodity flow (luồng chia tùy ý, không ràng buộc ECMP/SR)
# ----------------------------------------------------------------------------

"""
Mọi segment routing đều sinh ra một luồng đa hàng hóa hợp lệ, nên min-MLU của bài
toán luồng đa hàng hóa (LP, gộp theo đỉnh nguồn) là một CẬN DƯỚI hợp lệ cho MLU của
T-ASR. Dùng để đánh giá chất lượng nghiệm cho các instance mà MIP không cho bound tốt.
"""
function mcf_lower_bound(inst::Instance, alives, periods; time_limit = 900.0, solver = "ipm")
    n = inst.net.n
    model = direct_model(HiGHS.Optimizer())
    set_silent(model)
    set_attribute(model, "time_limit", Float64(time_limit))
    set_attribute(model, "threads", 1)
    set_attribute(model, "solver", solver)     # LP lớn: interior point nhanh hơn simplex nhiều
    solver == "ipm" && set_attribute(model, "run_crossover", "off")   # crossover rất chậm ở LP lớn (setA-19/20)
    @variable(model, U >= 0)
    @objective(model, Min, 1.0 * U)     # AffExpr: direct_model HiGHS không nhận objective dạng biến đơn
    for (ti, t) in enumerate(periods)
        arcs = findall(alives[ti])
        srcs = unique(inst.src[d] for d in 1:ndemands(inst) if inst.vol[d, t+1] > 0)
        sidx = Dict(s => k for (k, s) in enumerate(srcs))
        supply = zeros(length(srcs), n)
        for d in 1:ndemands(inst)
            v = inst.vol[d, t+1]; v > 0 || continue
            k = sidx[inst.src[d]]
            supply[k, inst.src[d]] += v
            supply[k, inst.dst[d]] -= v
        end
        f = @variable(model, [1:length(srcs), 1:length(arcs)], lower_bound = 0.0)
        for k in eachindex(srcs)
            bal = [AffExpr(0.0) for _ in 1:n]
            for (q, a) in enumerate(arcs)
                add_to_expression!(bal[inst.net.arc_from[a]], 1.0, f[k, q])
                add_to_expression!(bal[inst.net.arc_to[a]], -1.0, f[k, q])
            end
            for w in 1:n
                @constraint(model, bal[w] == supply[k, w])
            end
        end
        for (q, a) in enumerate(arcs)
            @constraint(model, sum(f[k, q] for k in eachindex(srcs)) <= inst.net.capacity[a] * U)
        end
    end
    optimize!(model)
    termination_status(model) == MOI.OPTIMAL || return NaN
    return objective_value(model)
end
