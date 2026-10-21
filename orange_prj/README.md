# Orange "Keep the flow!" — segment routing (readme Steps 2–9)

| Step | File | Content |
|---|---|---|
| 2 | `src/graph.jl` | JSON net file → `Graphs.jl` digraph (course Pluto script, multigraph flag accepted) |
| 3 | `src/rcoef.jl` | `r(i,j,a,t) = σ_iu σ_vj / σ_ij` on shortest paths (`calculating_r.pdf`), stored sparsely |
| 4 | `src/data.jl` | Data section: `c(a)`, `r` from net file, `v(d,0)` from tm file, `maxSeg` from scenario file |
| 5 | `src/model_t0.jl` | MILP for t = 0: constraints (1), (2), (3), objective min MLU; no budget constraint (4) |
| 6 | `run.jl --model t0`, `run_batch.sh` | Experiments, period 0, 1800 s |
| 7 | `src/data_t1.jl` | `G_1 = G − q(1)`, `r(i,j,a,1)`, `v(d,1)`, budget `κ(1)` |
| 8 | `src/model_t01.jl` | MILP for t = 0 and 1 with the budget constraint (4) |
| 9 | `run.jl --model t01`, `run_batch.sh` | Experiments, periods 0 and 1, 3600 s |

## Setup

```bash
curl -fsSL https://install.julialang.org | sh      # Julia (juliaup)
cd orange_prj
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Run

```bash
julia --project=. test/runtests.jl                     # toy instance, r checks on all SetA, models
MODEL=t0  LP=choose ./run_batch.sh 1                   # Step 6, default LP solver
MODEL=t0  LP=ipm    ./run_batch.sh 2 3 4               # Step 6, interior point
MODEL=t01 LP=ipm PRESOLVE=off MEM_MAX=26G ./run_batch.sh 4 5   # Step 9, memory cap
```

`run_batch.sh` runs each instance once (seed 0); when that run ends Optimal
within the time limit it is repeated with seeds 1–9, so the tables can give
the average of 10 runs. Runs that reach the time limit are not repeated.
Each run is a separate Julia process; with `MEM_MAX`, a run that exceeds
the cap is killed and listed in `results/failures.csv`.

Results: `results/t0_results.csv`, `results/t01_results.csv` (one row per
run); solver logs: `results/logs/`; routings in the challenge output format
(Section 5.2): `results/solutions/`. `results/first_simplex_runs/` holds the
first period-0 runs, all with the default LP solver (dual simplex).

`cpu_time_s` is the duration of the `optimize!` call only (it includes
JuMP passing the model to HiGHS); `highs_time_s` is HiGHS's own run time.
`gap = 1 − best_bound / objective`. `mlu_check` (`mlu0_check`, `mlu1_check`)
recomputes the MLU from the extracted waypoints, independently of the MILP;
`budget_used` is dist(P⁰, P¹) of the extracted routings.

## LP solver inside branch-and-bound

`--mip-lp-solver ipm` (`LP=ipm`) makes HiGHS use interior point for its LPs.
With the default (dual simplex) the root LP of setA-03…09 is not solved
within 1800 s; with `ipm` the LP relaxation of setA-03 solves in 3 s and of
setA-04 in 78 s. On setA-01 the default is faster (6 s vs 38 s).

## Presolve on the periods 0–1 model

The budget row (4) has one term per (d,i,j) in both periods, and HiGHS's
presolve time grows roughly quadratically with it: 11 s on setA-01, 83 s on
setA-02, 97 s on setA-03, and on setA-04 (1.47M columns) it was still in
presolve after 67 min, ignoring the time limit. Presolve removes less than
1 % of the rows and columns here, so setA-04 onwards run with
`--presolve off` (`PRESOLVE=off`); it then reaches branch-and-bound in 7 s.

## Modelling notes

- **Objective.** `min U` with `load(a,t) ≤ U c(a)`: U is the maximum link
  utilization over all arcs and periods, the first level of the
  lexicographic objective (5).
- **Constraint (3)** is written through `f_ij = Σ_d v(d,t) x^dt_ij`, the total
  volume on segment (i,j). Same constraint after substitution, far fewer
  non-zeros. `--formulation direct` writes (3) literally (period 0); both
  give the same optimum (test on setA-01).
- **Budget (4)** uses `|x1 − x0| = (x1 − x0) + 2 max(0, x0 − x1)`: one
  variable `w ≥ x0 − x1, w ≥ 0` and one row per (d,i,j), plus
  `Σ (x1 − x0) + 2 Σ w ≤ κ(1)`.
- **Cycles.** Flow conservation (1) admits closed cycles of segments beside
  the s–t path. They only add load and are dropped when the paths are read;
  as they could make two routings look closer under (4), the budget is
  checked again on the extracted paths (`budget_used`).
- **Zero-volume demands** put no load on any arc: period 0 leaves them out;
  periods 0–1 keep a demand when it has volume in either period.
- **Warm start.** The routing without waypoints (in both periods: no change,
  so within any budget) is always feasible and is given to HiGHS as a MIP
  start (`--no-warm-start` disables it).
- **r on the toy instance** differs from Table 2 of the subject when an
  intermediate node splits: path counting gives r(v0,v5,a13) = 1/3, ECMP
  gives 1/2. In SetA fewer than 1.1 % of node pairs have several shortest paths.
- **Ties with non-integer metrics** (setA-05, 08, 12, 14): Graphs.jl's
  Dijkstra counts paths with exact float equality and misses ties such as
  0.0188 + 0.9585 vs 0.9773, so σ is recounted with a tolerance
  (`count_shortest_paths`). On integer metrics it equals Dijkstra's counts.
