# =============================================================================
#  model.jl — MILP for T-Adaptive Segment Routing (T-ASR)
#
#  Sets      V (vertices), A (arcs), D (demands), Tm (periods modelled: [0] or [0,1])
#  Variables x[d,i,j,t] ∈ {0,1}  segment (i,j) used by demand d at period t
#            λ ≥ 0               maximum link utilisation over all (a,t)
#            y[i,j,t] ≥ 0        aggregated traffic on segment (i,j)  (optional)
#            z[d,i,j] ≥ 0        |x[d,i,j,1] - x[d,i,j,0]|            (budget)
#  Objective min λ   (first level of the lexicographic objective (5) = MLU)
#  Constraints
#   (1) flow conservation of every demand in the complete segment graph
#   (2) Σ_ij x[d,i,j,t] ≤ MaxSeg
#   (3) Σ_d Σ_ij r[i,j,a,t] v(d,t) x[d,i,j,t] ≤ λ c(a)     ∀ a ∈ A, a ∉ q(t)
#   (4) Σ_d Σ_ij |x[d,i,j,t] - x[d,i,j,t-1]| ≤ κ(t)          (only if t=1 modelled)
# =============================================================================

using JuMP, HiGHS

"""
    build_model(inst; periods=[0], budget=true, timelimit=1800.0)

`aggregate = true` writes (3) through y[i,j,t] = Σ_d v(d,t) x[d,i,j,t]:
an exact reformulation with fewer nonzeros, but its LP was much slower in
tests (setA-02 root LP > 120 s vs 2.9 s), so the default is `false`.
"""
function build_model(inst::Instance; periods = [0], aggregate::Bool = false,
                     budget::Bool = true, timelimit = 1800.0, threads = 0,
                     verbose = false, warmstart::Bool = true)
    m = Model(HiGHS.Optimizer)
    set_time_limit_sec(m, timelimit)
    (budget && length(periods) > 1) && set_optimizer_attribute(m, "presolve", "off")   # HiGHS presolve ignores the time limit on the budget model
    set_optimizer_attribute(m, "mip_rel_gap", 1e-6)
    threads > 0 && set_optimizer_attribute(m, "threads", threads)
    verbose || set_silent(m)

    n  = inst.n
    D  = 1:length(inst.demands)
    P  = Dict(t => inst.periods[t+1] for t in periods)

    # usable segments at period t: i ≠ j and j reachable from i in G(t)
    S = Dict(t => [(i, j) for i in 1:n for j in 1:n if i != j && P[t].reach[i, j]]
             for t in periods)

    # ---------------- variables ----------------
    x = Dict{Int,Any}()
    for t in periods
        x[t] = @variable(m, [d in D, s in S[t]], Bin, base_name = "x$t")
    end
    @variable(m, λ >= 0)

    # ---------------- (1) flow conservation ----------------
    for t in periods
        out = [Tuple{Int,Int}[] for _ in 1:n]; inn = [Tuple{Int,Int}[] for _ in 1:n]
        for (i, j) in S[t]; push!(out[i], (i, j)); push!(inn[j], (i, j)); end
        for d in D
            s, tg = inst.demands[d].s, inst.demands[d].t
            for i in 1:n
                rhs = i == s ? 1 : (i == tg ? -1 : 0)
                @constraint(m, sum(x[t][d, e] for e in out[i]; init = 0) -
                               sum(x[t][d, e] for e in inn[i]; init = 0) == rhs)
            end
        end
    end

    # ---------------- (2) maximum number of segments ----------------
    for t in periods, d in D
        @constraint(m, sum(x[t][d, e] for e in S[t]) <= inst.maxseg)
    end

    # ---------------- (3) load constraints ----------------
    for t in periods
        vol(d) = inst.demands[d].vol[t+1]
        # arc k -> list of (segment, r)
        onarc = [Tuple{Tuple{Int,Int},Float64}[] for _ in inst.arcs]
        for e in S[t], (k, r) in P[t].seg[e]
            r > 0 && push!(onarc[k], (e, r))
        end
        if aggregate
            y = @variable(m, [e in S[t]], lower_bound = 0, base_name = "y$t")
            for e in S[t]
                @constraint(m, y[e] == sum(vol(d) * x[t][d, e] for d in D if vol(d) > 0; init = 0))
            end
            for (k, a) in enumerate(inst.arcs)
                P[t].active[k] || continue
                isempty(onarc[k]) && continue
                @constraint(m, sum(r * y[e] for (e, r) in onarc[k]) <= a.cap * λ)
            end
        else
            for (k, a) in enumerate(inst.arcs)
                P[t].active[k] || continue
                isempty(onarc[k]) && continue
                @constraint(m, sum(r * vol(d) * x[t][d, e] for (e, r) in onarc[k], d in D) <= a.cap * λ)
            end
        end
    end

    # ---------------- (4) budget constraints (linearised |·|) ----------------
    if budget
        for t in periods
            t == 0 && continue
            (t - 1) in periods || continue
            prev = Set(S[t-1]); cur = Set(S[t])
            allseg = collect(union(prev, cur))
            z = @variable(m, [d in D, e in allseg], lower_bound = 0, base_name = "z$t")
            for d in D, e in allseg
                xc = e in cur  ? x[t][d, e]   : 0
                xp = e in prev ? x[t-1][d, e] : 0
                @constraint(m, z[d, e] >= xc - xp)
                @constraint(m, z[d, e] >= xp - xc)
            end
            @constraint(m, sum(z[d, e] for d in D, e in allseg) <= inst.budget[t])
            warmstart && set_start_value.(z, 0.0)
        end
    end

    @objective(m, Min, λ)

    # MIP start: the trivial routing (no waypoint, x[d,s,t] = 1) is always
    # feasible (budget: identical in both periods → 0 changes).
    if warmstart
        for t in periods, d in D
            s, tg = inst.demands[d].s, inst.demands[d].t
            for e in S[t]
                set_start_value(x[t][d, e], e == (s, tg) ? 1.0 : 0.0)
            end
        end
        W0 = Dict(d => Int[] for d in D)
        set_start_value(λ, maximum(mlu_from_waypoints(inst, W0, t) for t in periods))
    end
    return m, x, λ, S
end

# -----------------------------------------------------------------------------
# Solve + collect the statistics required by the readme
# -----------------------------------------------------------------------------
function solve!(m)
    cpu = @elapsed optimize!(m)          # time of the optimize! call only
    st  = termination_status(m)
    has = primal_status(m) == FEASIBLE_POINT
    obj = has ? objective_value(m) : NaN
    lb  = try objective_bound(m) catch; NaN end
    nodes = try MOI.get(m, MOI.NodeCount()) catch; -1 end
    gap = (has && obj > 0) ? max(0.0, 1.0 - lb / obj) : NaN
    status = st == MOI.OPTIMAL ? "Optimal" :
             st == MOI.INFEASIBLE ? "Infeasible" :
             has ? "Not Optimal" : "No solution"
    return (status = status, term = st, obj = obj, lb = lb, gap = gap,
            nodes = nodes, cpu = cpu)
end

# Waypoints of demand d at period t.
# The x = 1 segments form an s→t trail plus possibly cycles (constraint (1)
# does not forbid them). Hierholzer's algorithm from s gives an Euler trail
# s→t through every selected segment connected to s, so the waypoint list
# reproduces the model's loads exactly; cycles disconnected from s are dropped
# (they only add load, r ≥ 0).
function waypoints(inst, x, S, d, t)
    s = inst.demands[d].s
    succ = Dict{Int,Vector{Int}}()
    for (i, j) in S[t]
        value(x[t][d, (i, j)]) > 0.5 && push!(get!(succ, i, Int[]), j)
    end
    stack = [s]; trail = Int[]
    while !isempty(stack)
        u = stack[end]
        if !isempty(get(succ, u, Int[]))
            push!(stack, pop!(succ[u]))
        else
            push!(trail, pop!(stack))
        end
    end
    reverse!(trail)                         # s, w1, ..., wl, t
    trail[end] == inst.demands[d].t || error("extracted trail does not end at target")
    return trail[2:end-1]
end

# dist(P,P') of the subject: number of segments that differ
function scheme_distance(inst, W0::Dict, W1::Dict)
    dist = 0
    for (d, dem) in enumerate(inst.demands)
        seg(w) = (q = [dem.s; w; dem.t]; Set((q[k], q[k+1]) for k in 1:length(q)-1))
        dist += length(symdiff(seg(W0[d]), seg(W1[d])))
    end
    return dist
end

# independent check: recompute every load from the waypoint lists
function mlu_from_waypoints(inst, W::Dict, t)
    P = inst.periods[t+1]
    load = zeros(length(inst.arcs))
    for (d, dem) in enumerate(inst.demands)
        seq = [dem.s; W[d]; dem.t]
        for q in 1:length(seq)-1
            for (k, r) in P.seg[(seq[q], seq[q+1])]
                load[k] += r * dem.vol[t+1]
            end
        end
    end
    return maximum(load[k] / inst.arcs[k].cap for k in eachindex(inst.arcs) if P.active[k])
end

function write_srpaths(io, inst, Wt::Dict)   # Wt[(d,t)] = waypoints (Julia vertices)
    print(io, "{\n  \"srpaths\": [")
    first = true
    for ((d, t), w) in sort(collect(Wt))
        isempty(w) && continue
        print(io, first ? "\n" : ",\n"); first = false
        print(io, "    { \"d\": $(d-1), \"t\": $t, \"w\": [",
              join(inst.vertex_to_id[w], ", "), "] }")
    end
    print(io, "\n  ]\n}\n")
end