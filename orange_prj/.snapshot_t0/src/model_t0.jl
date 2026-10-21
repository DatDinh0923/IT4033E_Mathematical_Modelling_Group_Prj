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

struct ModelT0
    model::Model
    inst::InstanceT0
    D::Vector{Int}                   # demands in the model: v(d,0) > 0
    segs::Vector{Tuple{Int, Int}}    # segment p = (i,j), i ≠ j
    out_segs::Vector{Vector{Int}}    # segments leaving i
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
    n = num_vertices(inst)
    A = 1:length(inst.arcs)
    D = [d for d in eachindex(inst.demands) if inst.v[d] > 0]
    K = length(D)

    # Segments S = {(i,j) : i ≠ j}
    segs = [(i, j) for i in 1:n for j in 1:n if i != j]
    P = length(segs)
    seg_id = zeros(Int, n, n)
    out_segs = [Int[] for _ in 1:n]
    in_segs = [Int[] for _ in 1:n]
    for (p, (i, j)) in enumerate(segs)
        seg_id[i, j] = p
        push!(out_segs[i], p)
        push!(in_segs[j], p)
    end

    model = Model(HiGHS.Optimizer)
    set_string_names_on_creation(model, P * K <= 100_000)

    @variable(model, x[1:P, 1:K], Bin)
    @variable(model, U >= 0)

    # (1) flow conservation: one segment path from s_d to t_d
    rhs(i, k) = i == inst.demands[D[k]][1] ? 1 : i == inst.demands[D[k]][2] ? -1 : 0
    @constraint(model, flow[i in 1:n, k in 1:K],
        sum(x[p, k] for p in out_segs[i]) - sum(x[p, k] for p in in_segs[i]) == rhs(i, k))

    # (2) at most maxSeg segments per demand
    @constraint(model, nseg[k in 1:K], sum(x[p, k] for p in 1:P) <= inst.maxSeg)

    # (3) load of arc a at most U c(a)
    if formulation == :aggregated
        @variable(model, f[1:P] >= 0)
        @constraint(model, volume[p in 1:P], f[p] == sum(inst.v[D[k]] * x[p, k] for k in 1:K))
        @constraint(model, load[a in A],
            sum(r * f[seg_id[i, j]] for (i, j, r) in inst.R[a]) <= inst.c[a] * U)
    elseif formulation == :direct
        @constraint(model, load[a in A],
            sum(r * inst.v[D[k]] * x[seg_id[i, j], k] for (i, j, r) in inst.R[a], k in 1:K)
            <= inst.c[a] * U)
    else
        error("formulation must be :aggregated or :direct")
    end

    # A segment (i,j) without any i→j path has r = 0 on every arc: its flow
    # would vanish from the network, so it is forbidden.
    for (p, (i, j)) in enumerate(segs)
        inst.sigma[i, j] == 0 && fix.(x[p, :], 0; force = true)
    end

    @objective(model, Min, U)

    if warm_start
        X0 = zeros(P, K)
        for k in 1:K
            s, t = inst.demands[D[k]]
            X0[seg_id[s, t], k] = 1.0
        end
        set_start_value.(x, X0)
        formulation == :aggregated && set_start_value.(model[:f], X0 * inst.v[D])
        set_start_value(U, mlu(inst, direct_paths(inst)))
    end

    return ModelT0(model, inst, D, segs, out_segs, x, U)
end

"""
    solve_t0!(m; time_limit = 1800.0, seed = nothing)

Runs the solver and returns what the result tables need. `cpu_time` is the
duration of the `optimize!` call only; `gap = 1 − (best lower bound)/(best
solution value)`.
"""
function solve_t0!(m::ModelT0; time_limit = 1800.0, seed = nothing)
    model = m.model
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

"""
    extract_paths(m)

Segment path (node sequence s, w_1, …, t) of every demand, read from x by
walking the chosen segments from s to t. Flow conservation also admits
closed cycles of segments next to the s–t path; they only add load, are
never needed, and are dropped. Returns the paths and the number of dropped
segments.
"""
function extract_paths(m::ModelT0)
    inst = m.inst
    X = value.(m.x)
    paths = direct_paths(inst)
    dropped = 0

    for (k, d) in enumerate(m.D)
        s, t = inst.demands[d]
        used = Set{Int}()
        path = [s]
        while path[end] != t
            next = findfirst(p -> X[p, k] > 0.5 && !(p in used), m.out_segs[path[end]])
            next === nothing && error("demand $d: no segment path from s to t in x")
            p = m.out_segs[path[end]][next]
            push!(used, p)
            push!(path, m.segs[p][2])
        end
        paths[d] = path
        dropped += count(>(0.5), @view X[:, k]) - length(used)
    end

    return paths, dropped
end

"Solution file in the challenge output format (subject, Section 5.2)."
function write_srpaths(file::AbstractString, inst::InstanceT0, paths; t = 0)
    ids = inst.network.vertex_to_id
    entries = [(d = d - 1, t = t, w = [ids[w] for w in path[2:end-1]])
               for (d, path) in enumerate(paths) if length(path) > 2]
    open(io -> JSON3.write(io, (srpaths = entries,)), file, "w")
    return file
end
