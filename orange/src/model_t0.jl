# Step 5 — MILP for time period 0.
#
# The budget constraint (4) is not part of this model: it bounds the changes
# between the routings of periods t-1 and t, and only period 0 exists here.
#
#   min  U
#   s.t. Σ_{j≠i} x^d_ij − Σ_{j≠i} x^d_ji = 1 if i = s_d, −1 if i = t_d, 0 otherwise   ∀ d, i   (1)
#        Σ_{i≠j} x^d_ij ≤ maxSeg                                                    ∀ d      (2)
#        Σ_d Σ_{i≠j} r(i,j,a,0) v(d,0) x^d_ij ≤ U c(a)                               ∀ a      (3)
#        x^d_ij ∈ {0,1},  U ≥ 0
#
# U is the maximum link utilization max_a λ(a,0): minimizing it is the first,
# dominant level of the lexicographic objective (5).
#
# formulation = :aggregated (default) writes (3) through the total volume
# routed on each segment, f_ij = Σ_d v(d,0) x^d_ij:
#        Σ_{i≠j} r(i,j,a,0) f_ij ≤ U c(a)
# It is (3) after substitution, with far fewer non-zeros:
# |D|·Σ_a|{(i,j): r>0}| → Σ_a|{(i,j): r>0}| + |D|·n(n−1)  (28.5M → 7.9M for setA-07).
# formulation = :direct writes (3) literally.

"Segment set S = {(i,j) : i ≠ j}: segment p is `segs[p]`, `id[i,j] = p`."
struct Segments
    segs::Vector{Tuple{Int, Int}}
    id::Matrix{Int}
    out::Vector{Vector{Int}}         # segments leaving i
    in::Vector{Vector{Int}}          # segments entering j
end

function Segments(n::Int)
    segs = [(i, j) for i in 1:n for j in 1:n if i != j]
    id = zeros(Int, n, n)
    out = [Int[] for _ in 1:n]
    in = [Int[] for _ in 1:n]
    for (p, (i, j)) in enumerate(segs)
        id[i, j] = p
        push!(out[i], p)
        push!(in[j], p)
    end
    return Segments(segs, id, out, in)
end

"""
    add_segment_paths!(model, inst, D, S, sigma; base_name = "x")

Binary variables x[p,k] = x^{D[k]}_{segs[p]} of one time period with
constraints (1) and (2): every demand D[k] gets a segment path from s to t
with at most maxSeg segments. A segment (i,j) without any i→j path in the
period (σ_ij = 0) has r = 0 on every arc: its flow would vanish from the
network, so it is forbidden.
"""
function add_segment_paths!(model, inst, D, S::Segments, sigma; base_name = "x")
    n, P, K = num_vertices(inst), length(S.segs), length(D)
    x = @variable(model, [1:P, 1:K], Bin, base_name = base_name)

    # (1) flow conservation: one segment path from s_d to t_d
    rhs(i, k) = i == inst.demands[D[k]][1] ? 1 : i == inst.demands[D[k]][2] ? -1 : 0
    @constraint(model, [i in 1:n, k in 1:K],
        sum(x[p, k] for p in S.out[i]) - sum(x[p, k] for p in S.in[i]) == rhs(i, k))

    # (2) at most maxSeg segments per demand
    @constraint(model, [k in 1:K], sum(x[p, k] for p in 1:P) <= inst.maxSeg)

    for (p, (i, j)) in enumerate(S.segs)
        sigma[i, j] == 0 && fix.(x[p, :], 0; force = true)
    end
    return x
end

"""
    add_load_constraints!(model, x, D, v, R, arcs, c, S, U; base_name = "f")

Constraint (3) of one period for the arcs in `arcs`, written through the
volume on each segment f_p = Σ_k v(D[k]) x[p,k]:
Σ_{(i,j)} r(i,j,a) f_ij ≤ U c(a). Returns f.
"""
function add_load_constraints!(model, x, D, v, R, arcs, c, S::Segments, U; base_name = "f")
    P, K = size(x)
    f = @variable(model, [1:P], lower_bound = 0, base_name = base_name)
    @constraint(model, [p in 1:P], f[p] == sum(v[D[k]] * x[p, k] for k in 1:K if v[D[k]] > 0))
    @constraint(model, [a in arcs], sum(r * f[S.id[i, j]] for (i, j, r) in R[a]) <= c[a] * U)
    return f
end

"Starting point: every demand on its direct segment (s,t), i.e. no waypoint."
function direct_start(inst, D, S::Segments)
    X0 = zeros(length(S.segs), length(D))
    for (k, d) in enumerate(D)
        s, t = inst.demands[d]
        X0[S.id[s, t], k] = 1.0
    end
    return X0
end

struct ModelT0
    model::Model
    inst::InstanceT0
    D::Vector{Int}                   # demands in the model: v(d,0) > 0
    S::Segments
    x::Matrix{VariableRef}           # x[p,k] = x^{D[k]}_{segs[p]}
    U::VariableRef
end

"""
    build_model_t0(inst; formulation = :aggregated, warm_start = true)

Demands with v(d,0) = 0 are left out: they put no load on any arc, so any
segment path is optimal for them (they keep the shortest path s → t).
With `warm_start`, the routing without waypoints is given to the solver as
a feasible starting solution.
"""
function build_model_t0(inst::InstanceT0; formulation = :aggregated, warm_start = true)
    D = [d for d in eachindex(inst.demands) if inst.v[d] > 0]
    S = Segments(num_vertices(inst))
    A = 1:length(inst.arcs)

    model = Model(HiGHS.Optimizer)
    set_string_names_on_creation(model, length(S.segs) * length(D) <= 100_000)

    @variable(model, U >= 0)
    x = add_segment_paths!(model, inst, D, S, inst.sigma)              # (1), (2)

    if formulation == :aggregated                                      # (3)
        f = add_load_constraints!(model, x, D, inst.v, inst.R, A, inst.c, S, U)
    elseif formulation == :direct
        @constraint(model, [a in A],
            sum(r * inst.v[D[k]] * x[S.id[i, j], k] for (i, j, r) in inst.R[a], k in eachindex(D))
            <= inst.c[a] * U)
    else
        error("formulation must be :aggregated or :direct")
    end

    @objective(model, Min, U)

    if warm_start
        X0 = direct_start(inst, D, S)
        set_start_value.(x, X0)
        formulation == :aggregated && set_start_value.(f, X0 * inst.v[D])
        set_start_value(U, mlu(inst, direct_paths(inst)))
    end

    return ModelT0(model, inst, D, S, x, U)
end

"""
    solve_model!(model; time_limit, seed = nothing)

Runs the solver and returns what the result tables need. `cpu_time` is the
duration of the `optimize!` call only; `gap = 1 − (best lower bound)/(best
solution value)`.
"""
function solve_model!(model::Model; time_limit, seed = nothing)
    set_time_limit_sec(model, time_limit)
    seed === nothing || set_attribute(model, "random_seed", seed)

    cpu_time = @elapsed optimize!(model)

    term = termination_status(model)
    has_sol = primal_status(model) == FEASIBLE_POINT
    obj = has_sol ? objective_value(model) : NaN
    bound = try objective_bound(model) catch; NaN end
    gap = has_sol && obj > 0 ? 1.0 - bound / obj : NaN

    status = term == OPTIMAL ? "Optimal" :
             term == INFEASIBLE ? "Infeasible" :
             has_sol ? "Not Optimal" : "No solution ($term)"

    return (status = status, termination = string(term), objective = obj, bound = bound,
            gap = gap, cpu_time = cpu_time, solve_time = solve_time(model),
            nodes = node_count(model))
end

solve_t0!(m::ModelT0; time_limit = 1800.0, seed = nothing) =
    solve_model!(m.model; time_limit, seed)

"""
    walk_paths(X, inst, D, S)

Segment path (node sequence s, w_1, …, t) of every demand, read from the
values X[p,k] by walking the chosen segments from s to t. Flow conservation
also admits closed cycles of segments next to the s–t path; they only add
load and are dropped. Demands outside D keep the direct path. Returns the
paths and the number of dropped segments.
"""
function walk_paths(X, inst, D, S::Segments)
    paths = direct_paths(inst)
    dropped = 0

    for (k, d) in enumerate(D)
        s, t = inst.demands[d]
        used = Set{Int}()
        path = [s]
        while path[end] != t
            next = findfirst(p -> X[p, k] > 0.5 && !(p in used), S.out[path[end]])
            next === nothing && error("demand $d: no segment path from s to t in x")
            p = S.out[path[end]][next]
            push!(used, p)
            push!(path, S.segs[p][2])
        end
        paths[d] = path
        dropped += count(>(0.5), @view X[:, k]) - length(used)
    end

    return paths, dropped
end

extract_paths(m::ModelT0) = walk_paths(value.(m.x), m.inst, m.D, m.S)

"""
    write_srpaths(file, inst, paths; t = 0)
    write_srpaths(file, inst, [paths_0, paths_1, …])

Solution file in the challenge output format (subject, Section 5.2).
"""
function write_srpaths(file::AbstractString, inst::InstanceT0, paths::Vector{Vector{Int}}; t = 0)
    return write_srpaths(file, inst, [paths]; first_t = t)
end

function write_srpaths(file::AbstractString, inst::InstanceT0,
                       periods::Vector{Vector{Vector{Int}}}; first_t = 0)
    ids = inst.network.vertex_to_id
    entries = [(d = d - 1, t = first_t + k - 1, w = [ids[w] for w in path[2:end-1]])
               for (k, paths) in enumerate(periods) for (d, path) in enumerate(paths)
               if length(path) > 2]
    open(io -> JSON3.write(io, (srpaths = entries,)), file, "w")
    return file
end
