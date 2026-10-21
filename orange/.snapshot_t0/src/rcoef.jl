# Step 3 — split coefficients r(i,j,a,t), following `project/calculating_r.pdf`.
#
#   r(i,j,a,t) = σ_iu(t) σ_vj(t) / σ_ij(t)   if d_iu(t) + ω_a + d_vj(t) = d_ij(t),  a = (u,v)
#              = 0                            otherwise
#
# where d = shortest-path distance and σ = number of shortest paths, both
# w.r.t. the arc metric ω (the "c_a(t)" of calculating_r.pdf; not to be
# confused with the capacity c(a) of the subject).

"""
    all_pairs_shortest_data(g, distmx)

Dijkstra from every source (calculating_r.pdf, "Julia implementation"):
`d[i,j]` is the shortest i→j distance, `sigma[i,j]` the number of shortest
i→j paths.
"""
function all_pairs_shortest_data(g, distmx)
    n = nv(g)

    d = fill(Inf, n, n)
    sigma = zeros(Float64, n, n)

    for i in vertices(g)
        state = dijkstra_shortest_paths(g, i, distmx; allpaths = true)
        d[i, :] .= state.dists
        sigma[i, :] .= state.pathcounts
    end

    return d, sigma
end

# Equality of path lengths. Some SetA instances (05, 08, 12, 14) have
# non-integer metrics, where d_iu + ω_a + d_vj and d_ij can differ by
# rounding even when the arc lies on a shortest path.
same_length(x, y; rtol = 1e-9) = abs(x - y) <= rtol * max(1.0, abs(y))

"""
    count_shortest_paths(g, distmx, d)

σ_ij recomputed from the distances `d` with the tie test `same_length`
(Case 2 of calculating_r.pdf, §1, with a tolerance): σ_ii = 1 and
σ_ij = Σ_{(u,j) ∈ A : d_iu + ω_uj = d_ij} σ_iu, taking j by increasing d_ij.

Needed because Dijkstra in Graphs.jl compares lengths exactly: with
non-integer metrics it misses real ties, e.g. 0.0188 + 0.9585 ≠ 0.9773 in
floating point (setA-08, setA-14), and r would then sum to 2 out of i.
On integer metrics the result equals `state.pathcounts`.
"""
function count_shortest_paths(g, distmx, d)
    n = nv(g)
    sigma = zeros(Float64, n, n)

    for i in vertices(g)
        sigma[i, i] = 1.0
        for j in sortperm(d[i, :])
            (j == i || !isfinite(d[i, j])) && continue
            for u in inneighbors(g, j)
                if same_length(d[i, u] + distmx[u, j], d[i, j])
                    sigma[i, j] += sigma[i, u]
                end
            end
        end
    end

    return sigma
end

"""
    split_coefficients(arcs, omega, d, sigma)

Sparse storage of r: `R[a]` lists every `(i, j, r_ija)` with `r_ija > 0`.
A dense n×n×|A| array would need 2.5 GB for setA-20, while the non-zeros
number only a few per pair (i,j).
"""
function split_coefficients(arcs, omega, d, sigma)
    n = size(d, 1)
    R = [Tuple{Int, Int, Float64}[] for _ in arcs]

    for (a, (u, v)) in enumerate(arcs)
        for j in 1:n, i in 1:n
            i == j && continue
            sigma[i, j] > 0 || continue          # no i→j path: r = 0
            if same_length(d[i, u] + omega[a] + d[v, j], d[i, j])
                push!(R[a], (i, j, sigma[i, u] * sigma[v, j] / sigma[i, j]))
            end
        end
    end

    return R
end

"""
    check_split_coefficients(n, arcs, R, sigma)

For every connected pair (i,j), one unit of flow leaves i and enters j:
Σ_{a ∈ δ⁺(i)} r_ija = Σ_{a ∈ δ⁻(j)} r_ija = 1. Returns the largest deviation.
"""
function check_split_coefficients(n, arcs, R, sigma)
    out_sum = zeros(n, n)
    in_sum = zeros(n, n)

    for (a, (u, v)) in enumerate(arcs), (i, j, r) in R[a]
        u == i && (out_sum[i, j] += r)
        v == j && (in_sum[i, j] += r)
    end

    worst = 0.0
    for j in 1:n, i in 1:n
        (i == j || sigma[i, j] == 0) && continue
        worst = max(worst, abs(out_sum[i, j] - 1), abs(in_sum[i, j] - 1))
    end

    return worst
end
