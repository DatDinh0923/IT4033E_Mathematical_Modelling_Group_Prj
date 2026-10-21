# Step 8 — MILP for time periods 0 and 1, with the budget constraint (4).
#
#   min  U
#   s.t. (1), (2) for t = 0 and t = 1: a segment path per demand and period
#        Σ_{i≠j} r(i,j,a,t) f^t_ij ≤ U c(a),  f^t_ij = Σ_d v(d,t) x^dt_ij      ∀ t, a ∈ A \ q(t)   (3)
#        Σ_d Σ_{i≠j} |x^d1_ij − x^d0_ij| ≤ κ(1)                                                  (4)
#        x^dt_ij ∈ {0,1},  U ≥ 0
#
# One U for both periods: U = max_{a,t} λ(a,t), the first level of (5).
#
# (4) is linearized with |x1 − x0| = (x1 − x0) + 2 max(0, x0 − x1):
#        Σ_{d,i≠j} (x^d1_ij − x^d0_ij) + 2 Σ_{d,i≠j} w^d_ij ≤ κ(1),   w^d_ij ≥ x^d0_ij − x^d1_ij,  w ≥ 0
# which needs one row per (d,i,j) instead of the two rows of y ≥ ±(x1 − x0).
#
# Like the subject's formulation, (1) does not exclude closed cycles of
# segments beside the s–t path. They only add load, but with the budget they
# could make two routings look closer than they are, so every solution is
# checked afterwards on the extracted paths (routing_distance).

struct ModelT01
    model::Model
    inst::InstanceT0
    p1::PeriodT1
    D::Vector{Int}                   # demands with v(d,0) > 0 or v(d,1) > 0
    S::Segments
    x0::Matrix{VariableRef}          # x0[p,k] = x^{D[k],0}_{segs[p]}
    x1::Matrix{VariableRef}          # x1[p,k] = x^{D[k],1}_{segs[p]}
    w::Matrix{VariableRef}
    U::VariableRef
end

"""
    build_model_t01(inst, p1; warm_start = true)

With `warm_start`, the routing without waypoints in both periods (no change,
so within any budget) is given to the solver as a feasible starting solution.
"""
function build_model_t01(inst::InstanceT0, p1::PeriodT1; warm_start = true)
    D = [d for d in eachindex(inst.demands) if inst.v[d] > 0 || p1.v[d] > 0]
    S = Segments(num_vertices(inst))
    P, K = length(S.segs), length(D)

    model = Model(HiGHS.Optimizer)
    set_string_names_on_creation(model, P * K <= 50_000)

    @variable(model, U >= 0)

    # (1), (2) in each period
    x0 = add_segment_paths!(model, inst, D, S, inst.sigma; base_name = "x0")
    x1 = add_segment_paths!(model, inst, D, S, p1.sigma; base_name = "x1")

    # (3) in each period, on the arcs in service
    f0 = add_load_constraints!(model, x0, D, inst.v, inst.R, eachindex(inst.arcs), inst.c, S, U;
                               base_name = "f0")
    f1 = add_load_constraints!(model, x1, D, p1.v, p1.R, p1.up, inst.c, S, U; base_name = "f1")

    # (4) budget
    w = @variable(model, [1:P, 1:K], lower_bound = 0, base_name = "w")
    @constraint(model, [p in 1:P, k in 1:K], w[p, k] >= x0[p, k] - x1[p, k])
    @constraint(model, sum(x1[p, k] - x0[p, k] + 2 * w[p, k] for p in 1:P, k in 1:K) <= p1.budget)

    @objective(model, Min, U)

    if warm_start
        X0 = direct_start(inst, D, S)
        set_start_value.(x0, X0)
        set_start_value.(x1, X0)
        set_start_value.(w, 0.0)
        set_start_value.(f0, X0 * inst.v[D])
        set_start_value.(f1, X0 * p1.v[D])
        direct = direct_paths(inst)
        set_start_value(U, max(mlu(inst, direct), mlu(inst, p1, direct)))
    end

    return ModelT01(model, inst, p1, D, S, x0, x1, w, U)
end

solve_t01!(m::ModelT01; time_limit = 3600.0, seed = nothing) =
    solve_model!(m.model; time_limit, seed)

"Segment paths of every demand in periods 0 and 1."
function extract_paths(m::ModelT01)
    paths0, _ = walk_paths(value.(m.x0), m.inst, m.D, m.S)
    paths1, _ = walk_paths(value.(m.x1), m.inst, m.D, m.S)
    return paths0, paths1
end
