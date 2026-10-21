# Orange "Keep the flow!": project overview, work done, results

This document is for the group. It explains the problem, the model, what was
implemented and run, and how to read every number in the result tables.
`README.md` is the code manual (how to install and run).

Contents:
1. [The problem](#1-the-problem)
2. [The mathematical model](#2-the-mathematical-model)
3. [What we did](#3-what-we-did)
4. [How to read the results: every metric explained](#4-how-to-read-the-results-every-metric-explained)
5. [Results](#5-results)
6. [What the results mean](#6-what-the-results-mean)
7. [Caveats to state in the report](#7-caveats-to-state-in-the-report)
8. [Possible next steps](#8-possible-next-steps)
9. [Files](#9-files)

---

## 1. The problem

### 1.1 Context

The challenge EURO/ROADEF 2026 is proposed by Orange (a French internet
provider). An internet backbone is a set of routers (nodes) connected by
links (arcs). Each link `a` has

- a **capacity** `c(a)`: how much traffic it can carry;
- a **metric** (weight) `ω(a)`: used by the routing protocol (IGP, e.g. OSPF).

By default, traffic from `s` to `t` follows the **shortest paths** for the
metric `ω`. When several shortest paths tie, the traffic is split between
them (ECMP, equal-cost multi-path). The weakness: shortest paths ignore
capacity, so some links get overloaded while others stay almost empty.

**Segment routing** fixes this without changing any weight. The operator can
tell a demand to travel through intermediate routers, called **waypoints**:
`s → w1 → w2 → … → t`. Each piece `(i, j)` is a **segment**, and inside a
segment the traffic still follows the shortest paths from `i` to `j`. Choosing
waypoints steers traffic away from busy links. A path may contain at most
`maxSeg` segments (6 in all SetA instances, i.e. at most 5 waypoints).

**Time periods.** Period `t = 0` is the normal network. At `t = 1` a
maintenance intervention takes some links down (one link in each SetA
instance), and the traffic volumes change. Traffic that used the missing link
moves elsewhere and may overload other links, so the routing should be
adapted. But every change is a manual reconfiguration of routers, with a cost
and a risk, so the number of changes from `t = 0` to `t = 1` is limited by a
**budget** `κ(1)`.

### 1.2 Data of an instance (3 JSON files in `../setA/`)

| File | Content |
|---|---|
| `setA-XX-net.json` | nodes; links with capacity `c(a)` and metric `ω(a)` |
| `setA-XX-tm.json` | traffic matrix: demands `d = (s, t)` with a volume `v(d, t)` for `t = 0` and `t = 1` |
| `setA-XX-scenario.json` | `maxSeg`, interventions (links down at `t = 1`), budget `κ(1)` |

The 9 instances we ran:

| Instance | Nodes | Links (directed) | Demands | with volume > 0 at t=0 / t=1 | Link down at t=1 | Budget κ(1) |
|---|---|---|---|---|---|---|
| setA-01 | 20 | 80 | 40 | 39 / 40 | 22 | 51 |
| setA-02 | 30 | 150 | 45 | 45 / 45 | 59 | 63 |
| setA-03 | 50 | 250 | 20 | 20 / 20 | 52 | 53 |
| setA-04 | 50 | 250 | 200 | 193 / 199 | 106 | 44 |
| setA-05 | 100 | 396 | 100 | 100 / 100 | 337 | 1 |
| setA-06 | 100 | 500 | 500 | 500 / 500 | 53 | 13 |
| setA-07 | 100 | 500 | 800 | 793 / 791 | 220 | 90 |
| setA-08 | 150 | 654 | 200 | 200 / 200 | 584 | 13 |
| setA-09 | 150 | 750 | 200 | 198 / 199 | 152 | 18 |

Instances 10–20 are larger still: their period-0 models alone would have 15 million
to 958 million binary variables and cannot even be built in 30 GB of memory.

### 1.3 Key quantities

- **Split coefficient** `r(i, j, a, t)`: the fraction of the traffic sent from
  `i` to `j` (on shortest paths, at period `t`) that crosses arc `a`. Following
  `calculating_r.pdf`: for an arc `a = (u, v)` on a shortest `i–j` path,
  `r = σ_iu · σ_vj / σ_ij`, where `σ_xy` is the number of shortest paths from
  `x` to `y`. Otherwise `r = 0`.
  *Note:* this "path counting" rule splits by number of paths. The subject's
  Table 2 uses per-router ECMP (each router splits equally among its next
  hops). They differ only when an intermediate router splits again: on the
  toy example `r(v0, v5, a13)` is 1/3 with path counting and 1/2 with ECMP. We
  follow the teacher's PDF. In SetA fewer than 1.1 % of node pairs have
  several shortest paths, so the difference is small.
- **Load** (link utilization) `λ(a, t) = (traffic on a at t) / c(a)`, with the
  traffic on `a` = Σ over demands and over the segments they use of
  `r(i, j, a, t) · v(d, t)`. `λ = 0.8` means the link is 80 % full; `λ > 1`
  means it is **overloaded** (congestion, packet loss).
- **MLU** (maximum link utilization): the largest load over all links. This is
  the quantity we minimize.
- **dist(P, P′)**: the number of changes between two routings, i.e. the number
  of segments a demand uses in one routing but not in the other, summed over
  demands. Examples from the subject (Table 1): removing one waypoint
  `⟨s,w,t⟩ → ⟨s,t⟩` costs 3 (2 segments removed, 1 added); swapping two
  waypoints costs 6; replacing a waypoint by another costs 4. The smallest
  real change therefore costs 3.

### 1.4 The objective

The subject's objective (5) is **lexicographic**. Over all links and both
periods, first minimize the largest load, then the second largest, and so on.
We implement its **first level**: minimize `U`, the largest load over all
links (and over both periods when there are two). This is the readme's model,
and the first level is the dominant criterion. Doing the next levels needs a
sequence of further solves; see Section 8.

---

## 2. The mathematical model

A **MILP** (mixed-integer linear program), as written in the subject §3.

**Variables**
- `x^dt_ij ∈ {0, 1}`: 1 if demand `d` uses segment `(i, j)` at period `t`,
  for every ordered pair of nodes `i ≠ j`;
- `U ≥ 0`: the maximum link utilization.

**Constraints**
- **(1) Flow conservation.** For each demand and period, the chosen segments
  form a path from `s` to `t`: one segment leaves `s`, one enters `t`, and every
  other node has as many chosen segments in as out.
- **(2) maxSeg.** At most `maxSeg` segments per demand and period.
- **(3) Capacity.** For every link `a` in service at period `t`:
  `Σ_d Σ_ij r(i,j,a,t) · v(d,t) · x^dt_ij ≤ U · c(a)`.
- **(4) Budget** (periods 0+1 only): `Σ_d Σ_ij |x^d1_ij − x^d0_ij| ≤ κ(1)`.

**Objective:** `min U`.

**Size.** There is one binary per demand, ordered pair of nodes and period:
`|D| · n(n−1)` per period. This grows very fast: 14,820 binaries for setA-01
but 7.85 million for setA-07 (period 0 only).

**Implementation choices** (all exact: they do not change the optimum)

| Choice | Why |
|---|---|
| (3) written with `f_ij = Σ_d v(d,t) x^dt_ij` (total volume on segment `(i,j)`) | Same constraint after substitution, far fewer nonzeros, so a smaller model. A test checks that it gives the same optimum as the literal (3). |
| `|x1 − x0|` in (4) via `w ≥ x0 − x1, w ≥ 0` and `Σ(x1 − x0) + 2Σw ≤ κ` | Linear; needs one extra row per term instead of the usual two. |
| Zero-volume demands removed (period 0) | They put no load anywhere. |
| `x = 0` fixed when `i` cannot reach `j` | Never usable. |
| Warm start: the routing without waypoints | Always feasible (in both periods it means "no change", so within any budget), so HiGHS starts with a solution. |
| Shortest-path counting with a tolerance | With non-integer metrics (setA-05, 08, 12, 14), Graphs.jl's Dijkstra misses ties such as 0.0188 + 0.9585 vs 0.9773, which gave wrong `r`. |

**Why this method and this solver.** A MILP solved by branch-and-bound
(lectures 3–5) is an *exact* method. It returns a solution and also a
**lower bound**, i.e. a proof of how far that solution can be from the best
possible. It is the model the readme asks for. We use **HiGHS**, the best
open-source MILP solver, because Gurobi and CPLEX are not available to us.

---

## 3. What we did

### 3.1 Implementation (readme steps)

| Step | What | File |
|---|---|---|
| 2 | JSON network → graph (teacher's Pluto script; its "multigraph" error removed: 10 instances carry the flag but have no parallel links) | `src/graph.jl` |
| 3 | Split coefficients `r(i,j,a,0)` | `src/rcoef.jl` |
| 4 | Data section: `c(a)`, `r`, `v(d,0)`, `maxSeg` | `src/data.jl` |
| 5 | Model for period 0: (1), (2), (3), min U, no budget | `src/model_t0.jl` |
| 6 | Experiments, period 0, 1800 s | `run.jl --model t0`, `run_batch.sh` |
| 7 | Period 1: graph without the downed link, `r(i,j,a,1)`, `v(d,1)`, `κ(1)` | `src/data_t1.jl` |
| 8 | Model for periods 0 and 1, with the budget (4) | `src/model_t01.jl` |
| 9 | Experiments, periods 0 and 1, 3600 s | `run.jl --model t01`, `run_batch.sh` |

**Step 9's wording.** The readme says "for only time period 0" in Step 9.
That is a copy-paste error: Steps 7–8 build period 1 and the budget, so
Step 9 runs periods 0 **and** 1.

### 3.2 Checks (`test/runtests.jl`)

- **Toy example of the subject (§4–5):** `r` values, a hand-computed load, the
  model solves to optimality.
- **Period 1 on the toy:** links down, new `r` (matches Figure 3), the
  `dist` examples of Table 1, and the budget respected.
- **All 20 SetA instances:** graphs strongly connected; for every pair
  `(i, j)` the `r` values on the arcs leaving `i` (and on those entering `j`)
  sum to 1; and path counts equal Graphs.jl's on integer metrics.
- **Formulations:** the aggregated and literal constraint (3) give the same
  optimum on setA-01.
- **Every run:** the MLU is recomputed from the extracted paths, independently
  of the solver (`mlu_check`), and the budget is recomputed on the extracted
  routings (`budget_used`).

### 3.3 Experimental protocol

- **Time limits:** 1800 s for period 0 (Step 6), 3600 s for periods 0+1 (Step 9).
- **CPU time** = duration of the `optimize!` call only, as the readme asks.
  Julia's compilation is excluded: each run first solves the toy instance.
- **10-run rule (agreed in the group):** each instance is run once (seed 0).
  If that run ends Optimal within the limit, it is repeated with seeds 1–9 and
  the tables give the average of 10 runs. A run that reaches the time limit is
  not repeated: 10 runs of 30–60 min each would cost a day per instance.
- **One Julia process per run**, each with a memory cap. A run that exceeds the
  cap is killed and listed in `results/failures.csv`, instead of freezing the
  machine.
- **Machine:** AMD Ryzen 7 9800X3D (8 cores), 30 GB RAM, Ubuntu; Julia 1.13.2,
  JuMP 1.32.1, HiGHS 1.15.1. HiGHS's branch-and-bound is essentially serial
  here ("parallel search off").

### 3.4 Timeline of the experiments

1. **First period-0 runs, default LP solver (dual simplex)**, 01–09 (kept in
   `results/first_simplex_runs/`). Only 01 and 02 were solved. On 03–09 HiGHS
   never finished the first LP relaxation in 30 min, so it had no lower bound
   (gap 100 %).
2. **Switched to the interior point method** (`--mip-lp-solver ipm`) for
   02–09. It solves those LP relaxations (setA-03 in 3 s, setA-04 in 78 s in
   a first test). The final period-0 table uses it, except for 01 where
   the default is faster (6 s vs 38 s).
3. **Periods 0+1 runs** (01–03), in parallel with step 2.
4. **Presolve problem.** On setA-04 (periods 0+1), HiGHS's presolve hung for
   67 min and ignored the time limit. Presolve time grows roughly with the
   square of the budget row (one row with a term for every demand and segment
   in both periods): 11 s on 01, 83 s on 02, 97 s on 03, hours on 04. We added
   `--presolve off` for 04 onwards. Presolve removed less than 1 % of the model
   anyway.
5. **Reboot.** The machine had to restart. Runs cannot be paused (HiGHS's time
   limit is wall-clock), so they were stopped cleanly and only the missing ones
   were re-run.
6. **Overnight queue (21 Oct, 01:29–05:32), one run at a time** so that the
   big runs get the whole memory: periods 0+1 for 04, 05; period 0 for 06, 08,
   09, 07; periods 0+1 for 06, 08, 09, 07.

---

## 4. How to read the results: every metric explained

### 4.1 First, how branch-and-bound works

HiGHS keeps two numbers during the search:

```
          0        lower bound                     best solution (incumbent)
          |------------[=========== gap ===========]------------------> U
                 "no routing can do         "this routing exists and
                  better than this"          achieves this MLU"
```

- The **best solution** (incumbent) is a real routing it has found. It can
  only go down over time.
- The **lower bound** comes from the **LP relaxation**: the same model with
  `0 ≤ x ≤ 1` instead of `x ∈ {0,1}`, plus cuts. It is a proof that no routing
  can do better. It can only go up over time.
- To close the gap, HiGHS **branches**: it fixes some `x` to 0 or 1, creating
  sub-problems (nodes of the search tree).
- When the two numbers meet (relative gap ≤ 0.01 %), the solution is proven
  optimal.

### 4.2 Metric by metric

| Metric | What it is | Why we report it / how to read it |
|---|---|---|
| **Vertices, links, demands** | Size of the instance (Section 1.2). | Asked by the readme; explains the model size. |
| **Binary vars** (`n_binary_vars`) | Number of 0/1 variables `x` in the model: demands × n(n−1) × periods. | The main measure of difficulty: the solver's work grows quickly with it. It explains why 01–02 are solved and 06–09 are not, and why memory runs out. |
| **Constraints** (`n_constraints`) | Number of rows of the model. | Secondary size measure. |
| **Runs** | How many times the instance was solved (different seeds). "10 (8 Optimal)" = 10 runs, 8 of them proven optimal. | The solver's running time varies with its random seed. Averaging 10 runs gives a fair time; the spread shows how reliable it is. Only instances solved within the limit get 10 runs (Section 3.3). |
| **Seed** | The random seed given to HiGHS. | Changes its tie-breaking and heuristics, not the problem. Same seed means the same run: our setA-01 re-run on an idle machine took the same time per seed (3.05 s, 16.2 s, …). |
| **Status** | `Optimal`: proven best (gap ≤ 0.01 %). `Not Optimal`: time limit reached; the result is the best found, not proven. `Infeasible`: no routing exists (never happens here). **Killed**: the run exceeded the memory cap, no result. | Says whether the MLU in the table is *the* optimum or only an upper bound on it. |
| **No waypoints** (`baseline_mlu`) | The MLU when every demand goes straight from `s` to `t` on shortest paths, with no waypoint. That is the network's default behaviour. For periods 0+1: the larger of the two periods. | The **reference point**. "Best MLU" below it shows what segment routing gains. It is also HiGHS's starting solution, so "Best = No waypoints" means the solver found nothing better. |
| **Best MLU** / **Best U** (`objective`) | The MLU of the best routing found. Period 0: max load over all links. Periods 0+1: `U` = max load over all links **and both periods**. | The answer. Lower is better; above 1 means a link is overloaded. |
| **Lower bound** (`best_bound`) | Proven lower bound on the optimal MLU (Section 4.1). "none" = the LP relaxation was not solved in time. | Tells how good the answer could still become: the true optimum lies between the lower bound and the best MLU. |
| **Gap** | `1 − lower bound / best MLU` (the readme's definition). 0 = proven optimal; 33 % = the best MLU may be up to 33 % above the optimum; 100 % = no bound at all. | The quality guarantee of a non-optimal result. Large gaps can mean the solution is poor or only that the bound is weak; we cannot tell which. |
| **CPU time** (`cpu_time_s`) | Duration of the `optimize!` call, in seconds. Includes JuMP handing the model to HiGHS (90–170 s for the largest models), so time-limited runs show slightly more than 1800/3600 s. `highs_time_s` = HiGHS's own time only. | Asked by the readme. For Optimal runs it is the time to prove optimality. For time-limited runs it just equals the limit. |
| **B&B nodes** (`bb_nodes`) | Number of nodes of the branch-and-bound tree explored. 1 = solved at the root (the first LP plus cuts and heuristics were enough). 0 = the search never got past the root. | Asked by the readme. Shows how much search was done. Here, 0 nodes means the solver spent the whole time on the first LP and its cuts. |
| **LP solver** (`mip_lp_solver`) | Algorithm for the LP relaxations: `default` (HiGHS's choice = dual simplex) or `ipm` (interior point). | Dual simplex stalls on our large, degenerate LPs; interior point solves them (Section 5.1). |
| **Presolve** | HiGHS's model simplification before solving: on (default) or off. | Off for periods 0+1 from setA-04 on, because it hangs on the budget row (Section 3.4). |
| **Peak RAM** (`peak_rss_gb`) | Maximum memory used by the run. | Shows why the larger two-period models could not run: they need more than the 26 GB cap. |
| **MLU period 0 / period 1** (`mlu0_check`, `mlu1_check`) | For periods 0+1: the MLU of each period separately, recomputed from the paths. | `Best U` is the max of the two, and the two show which period is the bottleneck. |
| **Budget κ(1)** | Maximum number of routing changes allowed between periods (Section 1.3, `dist`). | A small budget forces period 1 to reuse period 0's routing. setA-05 has κ = 1 < 3, so no change at all is possible. |
| **Budget used** (`budget_used`) | `dist(P⁰, P¹)` of the routings found, recomputed from the paths. | Must be ≤ κ(1) (it always is). Shows whether the budget is tight. 0 = the same routing in both periods. |
| **mlu_check** | MLU recomputed from the extracted waypoints with our own code. | Independent check that the solver's answer is a real routing with that MLU. |
| **dropped_segments** | Segments in closed cycles that the model allows beside the s–t path; they are dropped when reading the paths. | Cycles only add load, so dropping them never makes the MLU worse (hence `mlu_check ≤ objective`). |

---

## 5. Results

### 5.1 Step 6: period 0 only, time limit 1800 s

Instances with 10 runs show averages over the 10 runs; the others show one run.

| Instance | Vertices | Links | Demands | Binary vars | LP | Runs (Optimal) | Status | No waypoints | Best MLU | Lower bound | Gap | CPU time (s) | B&B nodes | Peak RAM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 | 20 | 80 | 40 | 14,820 | default | 10 (10) | Optimal | 1.0000 | **0.9294** | 0.9294 | 0 % | 7.5 (3.1–16.2) | 1 | 0.7 GB |
| 02 | 30 | 150 | 45 | 39,150 | ipm | 10 (8) | Optimal ×8 | 1.0000 | **0.5486** (mean 0.5550) | 0.5486 | mean 1.1 % | 933.8 (60.5–1800.7) | 360 (1–1031) | 1.1 GB |
| 03 | 50 | 250 | 20 | 49,000 | ipm | 1 | Not Optimal | 1.0000 | **0.9435** | 0.6337 | 32.8 % | 1800.9 | 690 | 1.1 GB |
| 04 | 50 | 250 | 200 | 472,850 | ipm | 1 | Not Optimal | 0.6918 | 0.6918 | 0.4591 | 33.6 % | 1808.5 | 0 | 3.0 GB |
| 05 | 100 | 396 | 100 | 990,000 | ipm | 1 | Not Optimal | 0.3333 | 0.3333 | 0.1155 | 65.3 % | 1818.5 | 0 | 6.3 GB |
| 06 | 100 | 500 | 500 | 4,950,000 | ipm | 1 | Not Optimal | 1.0000 | 1.0000 | 0.0986 | 90.1 % | 1897.3 | 0 | 12.9 GB |
| 07 | 100 | 500 | 800 | 7,850,700 | ipm | 1 | Not Optimal | 0.9080 | 0.9080 | none | 100 % | 1970.4 | 0 | 19.2 GB |
| 08 | 150 | 654 | 200 | 4,470,000 | ipm | 1 | Not Optimal | 1.0000 | 1.0000 | 0.3189 | 68.1 % | 1887.9 | 0 | 14.1 GB |
| 09 | 150 | 750 | 200 | 4,425,300 | ipm | 1 | Not Optimal | 0.6111 | 0.6111 | 0.5054 | 17.3 % | 1887.2 | 0 | 11.6 GB |

setA-02 in detail: the 8 Optimal runs took 60.5–1357.8 s (mean 717.0 s). The
2 others (seeds 8, 9) hit 1800 s with MLU 0.5873 and 0.5738 (gaps 6.6 %,
4.4 %); their lower bound had already reached 0.5486.

**Dual simplex vs interior point** (period 0, first run of each):

| Instance | Dual simplex (first runs) | Interior point (final runs) |
|---|---|---|
| 01 | Optimal, 5.7 s | Optimal, 38 s (single test) → we keep simplex |
| 02 | Optimal, 1352.7 s, 103 nodes | 8/10 Optimal, mean 717 s over those 8 |
| 03 | 1.0000, no bound (gap 100 %) | **0.9435**, gap 32.8 % |
| 04 | 0.6918, no bound | 0.6918, gap 33.6 % |
| 05 | 0.3333, no bound | 0.3333, gap 65.3 % |
| 06 | 1.0000, no bound | 1.0000, gap 90.1 % |
| 07 | 0.9080, no bound | 0.9080, no bound |
| 08 | 1.0000, no bound | 1.0000, gap 68.1 % |
| 09 | 0.6111, no bound | 0.6111, gap 17.3 % |

### 5.2 Step 9: periods 0 and 1 with the budget, time limit 3600 s

| Instance | Binary vars | κ(1) | LP / presolve | Runs (Optimal) | Status | No waypoints | Best U | Lower bound | Gap | CPU time (s) | B&B nodes | MLU t=0 | MLU t=1 | Budget used | Peak RAM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 | 30,400 | 51 | default / on | 10 (10) | Optimal | 1.0000 | **0.9294** | 0.9294 | 0 % | 279.2 (19.2–947.1) | 11.5 (1–106) | 0.9294 | 0.831 (0.758–0.918) | 31.2 (18–47) | 1.9 GB |
| 02 | 78,300 | 63 | ipm / on | 1 | Not Optimal | 1.0981 | 1.0981 | 0.6243 | 43.1 % | 3601.5 | 0 | 1.0000 | 1.0981 | 0 | 2.9 GB |
| 03 | 98,000 | 53 | ipm / on | 1 | Not Optimal | 1.0000 | 1.0000 | 0.6337 | 36.6 % | 3601.9 | 0 | 1.0000 | 0.6019 | 0 | 3.8 GB |
| 04 | 980,000 | 44 | ipm / off | 1 | Not Optimal | 1.0000 | 1.0000 | none | 100 % | 3620.1 | 0 | 0.6918 | 1.0000 | 0 | 5.7 GB |
| 05 | 1,980,000 | 1 | ipm / off | 1 | **Killed** at 1903 s (> 26 GB) | 1.0000 | (1.0000) | (0.2050, from the log) | – | – | – | – | – | – | > 26 GB |
| 06 | 9,900,000 | 13 | ipm / off | 1 | **Killed** at 344 s (> 26 GB) | – | – | – | – | – | – | – | – | – | > 26 GB |
| 07 | 15,840,000 | 90 | ipm / off | 1 | **Killed** at 401 s (> 26 GB) | – | – | – | – | – | – | – | – | – | > 26 GB |
| 08 | 8,940,000 | 13 | ipm / off | 1 | **Killed** at 311 s (> 26 GB) | – | – | – | – | – | – | – | – | – | > 26 GB |
| 09 | 8,940,000 | 18 | ipm / off | 1 | **Killed** at 303 s (> 26 GB) | – | – | – | – | – | – | – | – | – | > 26 GB |

For 01 the values in brackets are min–max over the 10 runs; the others are
means. setA-05 had solved its first LP relaxation (bound 0.2050 at 1829 s)
and was killed while HiGHS added cuts, so it has no CSV row. 06, 08 and 09
were killed while solving the first LP relaxation, and 07 before HiGHS
started solving.

---

## 6. What the results mean

**1. The exact MILP works only on small instances.** HiGHS proves optimality
on 01 and 02 for period 0 (15k and 39k binaries) and only on 01 for both
periods (30k). From 04 on (and already on 02–03 with both periods), the search
never gets past the root node: the best
solution stays the starting routing without waypoints, and the solver only
contributes a lower bound. The cause is size: the model has one binary per
demand and *pair of nodes*, i.e. millions of variables, and even the first LP
relaxation takes most of the time limit.

**2. Segment routing helps a lot when it is found.**
- 02: from 1.0000 (a link full) to 0.5486, the busiest link drops from 100 % to
  55 %.
- 01: from 1.0000 to 0.9294.
- 03: from 1.0000 to 0.9435.

**3. Interior point is the right LP solver for these models.** Dual simplex,
HiGHS's default, stalls on these large degenerate LPs and gave no bound at all
on 03–09. Interior point gives bounds on all of them except 07, and finds the
improvement on 03. Simplex remains faster on the tiny 01.

**4. The lower bounds still say something useful.**
- setA-09: no routing can beat 0.5054, and we have 0.6111, so at most 17 % can
  still be gained.
- setA-06: the bound (0.0986) is very weak, so we know almost nothing there.
- A large gap does not prove the solution is bad. It only shows the relaxation
  is too weak to tell.

**5. Adding period 1 and the budget makes the problem much harder.**
- 01 takes 279 s on average instead of 7.5 s. The two periods are coupled by
  the budget, and the model is twice as big plus one extra variable and row
  per term.
- 01's optimum is the same as for period 0 alone (0.9294): the bottleneck is
  in period 0, and period 1 has room to spare.
- Within the optimal solutions, period 1's MLU varies from 0.758 to 0.918
  between seeds: any period-1 routing below 0.9294 counts as optimal for
  `min U`. The second level of the lexicographic objective (5) is what would
  choose among them.

**6. Maintenance without segment routing overloads the network.**
- 02: after its link goes down at t = 1, the no-waypoint routing loads a link
  at 109.8 %, so that link is congested.
- 04 and 05: period 1 reaches 100 % (from 69 % and 33 % in period 0).
- This is exactly the situation the challenge is about. Our solver could not
  fix it within an hour on 02–04. On 02 the optimum is somewhere between
  0.6243 and 1.0981.

**7. Practical limits of the exact approach.**
- Presolve hangs on the dense budget row (solved by turning it off).
- The two-period models of 05–09 need more than 26 GB of memory.
- These limits are a result in themselves: they are why the challenge
  instances call for decomposition (column generation) or heuristics.

---

## 7. Caveats to state in the report

- Time-limited runs were run once (10-run rule), so their numbers are not
  averages.
- Some runs shared the machine with one other job (two at a time on 8 cores):
  01–05 period 0, 01–03 periods 0+1, and the first simplex runs. The overnight
  runs (06–09 period 0, 04–09 periods 0+1) ran alone. Re-running setA-01 on an idle machine gave the same time for every
  seed, so the effect looks small.
- CPU times include the transfer of the model from JuMP to HiGHS.
- Presolve was on for 01–03 (periods 0+1) and off for 04–09. On 02–03 presolve
  cost about 90 s of the hour.
- The LP solver is the default (dual simplex) for 01 and interior point for
  02–09.
- setA-05 (periods 0+1): the bound 0.2050 comes from the solver log, not the CSV.
- `r` is computed by path counting (`calculating_r.pdf`), not per-router ECMP
  (subject Table 2). See Section 1.3.

---

## 8. Possible next steps

- **Better starting point for periods 0+1:** start from the period-0 optimum,
  used in both periods. It costs 0 budget, and on 02/03 it would immediately
  give a better `U` than the no-waypoint routing the solver is stuck on.
- **Column generation** (lecture 6). Use one variable per *segment path*
  instead of per segment, and generate only useful paths. The pricing problem
  is a shortest path with at most `maxSeg` segments (Bellman–Ford on segment
  costs `Σ_a π_a r(i,j,a)`). This avoids building millions of variables, and is
  the standard approach in the literature (e.g. Jadin et al., INFOCOM 2019).
- **Heuristics:** allow at most one waypoint per demand (2-segment routing,
  much smaller), or large neighbourhood search (re-optimize a few demands at a
  time).
- **Second level of the lexicographic objective:** fix `U` at its optimum and
  minimize the next largest load. This applies to 01, where period 1 has slack.

---

## 9. Files

| Path | Content |
|---|---|
| `results/t0_results.csv` | Step 6, one row per run (final runs) |
| `results/t01_results.csv` | Step 9, one row per run |
| `results/first_simplex_runs/` | First period-0 runs with dual simplex (for the comparison) |
| `results/t0_setA-01_quiet.csv` | setA-01 re-run on an idle machine (check of timing variation) |
| `results/failures.csv` | Runs that did not finish: presolve hang, memory cap |
| `results/logs/` | Full HiGHS log of every run: `setA-XX-<model>-<lp>[-presolve-off]-seed<k>.log` |
| `results/solutions/` | Best routing of each instance, in the challenge output format (§5.2) |
| `.snapshot_t0/` | Frozen copy of the code that produced the period-0 runs (same model as `src/`) |
| `README.md` | How to install and run |
