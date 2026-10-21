# IT4033E Mathematical Modelling: Orange EURO/ROADEF Challenge 2026

Group project for the Mathematical Modelling course (HUST, School of Information and Communications Technology), on the EURO/ROADEF Challenge 2026 subject proposed by Orange, *"Keep The Flow!"*.

**Team:** Le Tuan Huy (20242461M) · Dinh Quoc Dat (20252745M) · Nguyen Manh Hung (20252276M)

## The problem

For every traffic demand in an IP network, choose a short list of **waypoints** (Segment Routing) so that the **maximum link utilisation (MLU)** is as low as possible. Traffic between consecutive waypoints follows all shortest paths (ECMP). There are two periods: period 0 is the nominal network, and in period 1 one link is down for maintenance. The routing may change between the two periods, but the number of changed segments must stay within a budget κ(1).

We model it as a mixed-integer linear program (MILP) in Julia/JuMP and solve it with HiGHS:

| Readme step | What | Where |
|---|---|---|
| 2–4, 7 | Data section: graph G(t), split coefficients r(i,j,a,t), c(a), v(d,t), MaxSeg, κ(1) | `src/data.jl` |
| 5, 6 | Model for period 0 (no budget), 1800 s limit | `src/model.jl`, `run_experiments.jl p0` |
| 8, 9 | Joint model for periods 0 and 1 with the budget constraint, 3600 s limit | `src/model.jl`, `run_experiments.jl p01` |
| (extra) | Two-stage method: fix the period-0 routing, solve period 1 with a single budget row | `two_stage.jl` |

The split coefficients follow the course note `calculating_r.pdf`: r = σ_iu·σ_vj / σ_ij, from Dijkstra with path counting. The challenge's own ECMP rule (an even split at each node) is also implemented and reproduces Table 2 of the subject exactly. Select it with `RMETHOD=ecmp`.

## Repository layout

```
.
├── Project.toml, Manifest.toml   Julia environment (exact package versions)
├── src/
│   ├── data.jl                   reads the 3 JSON files, builds G(t), computes r
│   └── model.jl                  MILP (period 0, or periods 0–1 with budget), solve statistics, checks
├── run_experiments.jl            Steps 6 and 9: one CSV line per run + solution file
├── two_stage.jl                  two-stage method for Step 9 (needs the Step 6 results)
├── summarize.jl                  averages runs per instance → Markdown and LaTeX table rows
├── test_toy.jl                   checks on the toy instance of the subject
├── toy/                          toy instance (Section 4 of the subject)
├── results/                      result CSVs and solutions (see below)
├── logs/                         one log per run + kill logs of the guards
├── run_all.sh, run_p01_rest.sh,  batch scripts used for the experiments
│   run_p01_all.sh, run_ts.sh
├── watchdog.sh, timeguard.sh,    guards: kill a run above 26 GB of memory or past its time budget
│   guard_ts.sh
├── test.sh                       how the first watchdog was started
└── info.txt                      notes on the experiment runs (stopped runs, reruns)
```

## Setup

1. Install Julia (we used **1.13.2**).
2. Get the SetA instances from the course repository and place it **next to** this repository:
   ```
   work/
   ├── SOICT_mathematical_modelling/      git clone https://github.com/vhnguyen97/SOICT_mathematical_modelling.git
   └── IT4033E_Mathematical_Modelling_Group_Prj/   (this repo)
   ```
   Or point to the data folder anywhere with `export DATADIR=/path/to/setA`.
3. Install the packages (from the repository root):
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
4. Check everything on the toy instance:
   ```bash
   julia --project=. test_toy.jl
   ```
   Expected output with ECMP: `r(v0,v5,a35,0) = 0.75`, MLU without waypoints `0.4375`, optimum `0.375`. With path counts: `0.6667`, `0.4167`, `0.375`.

## Running the experiments

```bash
# Step 6: period 0, 1800 s limit, 10 runs
julia --project=. run_experiments.jl p0 1800 10 setA-01

# Step 9: joint two-period model with budget, 3600 s limit
julia --project=. run_experiments.jl p01 3600 1 setA-01

# Step 9: two-stage method (run Step 6 for the instance first)
julia --project=. two_stage.jl 3600 1 setA-01

# Table rows (averages over the runs of each instance)
julia --project=. summarize.jl results/p0-pathcount.csv
```

`run_experiments.jl <mode> <timelimit_s> <runs> [instances...]` runs all 12 instances setA-01 to 12 if no instance is given. Lines are appended to the CSV, so reruns add runs. The batch scripts `run_all.sh` and `run_ts.sh` run each instance in its own process with its own log. Start them with `nohup bash run_all.sh > /dev/null 2>&1 &`, together with a guard script in the same way.

## Results files

| File | Content |
|---|---|
| `results/p0-pathcount.csv` | Step 6, period 0. The first line is a 120 s test run (see `logs/2.log`); 10 runs for setA-01 to 03, 1 run for 04 to 09 |
| `results/p01-pathcount-presolve-on.csv` | Step 9, joint model with default HiGHS settings (setA-01 to 04) |
| `results/p01-pathcount.csv` | Step 9, joint model with HiGHS presolve off (setA-01 to 05) |
| `results/p01ts-pathcount.csv` | Step 9, two-stage method (setA-01 to 06 and 08) |
| `results/<inst>-<mode>-pathcount-srpaths.json` | Solution of run 1 in the challenge output format (waypoints per demand and period) |
| `logs/watchdog.log`, `logs/guard_ts.log` | Runs killed for memory (26 GB) or time (100 min) |

Every CSV line includes the readme's columns (`V`, `A`, `D`, `obj`, `lb`, `gap` = 1 − LB/UB, `nodes`, `cpu` = time of `optimize!` only, `status`) and three independent checks:

- `mlu_check`: the MLU recomputed from the extracted waypoints.
- `dist_check`: the number of changes between the periods, which must be ≤ budget.
- `maxrss_MB`: the peak memory of the run.

## Main results

Environment: Lenovo ThinkBook 14 G6 IMH, Intel Core Ultra 7 155H, 32 GB RAM (30.4 GB usable), Ubuntu 26.04.1 LTS, Julia 1.13.2, JuMP 1.32.1, HiGHS 1.15.1.

**Step 6, period 0 (1800 s):**

| Instance | MLU | Gap | CPU (s) | Status |
|---|---|---|---|---|
| setA-01 | 0.9294 | 0 | 3.8 | Optimal (10 runs) |
| setA-02 | 0.5486 | 0 | 33.9 | Optimal (10 runs) |
| setA-03 | 0.9435 | 0 | 14.9 | Optimal (10 runs) |
| setA-04 | 0.6918 | 33.6% | 1812 | Not optimal |
| setA-05 | 0.1656 | 30.3% | 1827 | Not optimal |
| setA-06 | 1.0000 | 100% | 1955 | Not optimal |
| setA-07 | 0.9080 | 100% | 2075 | Not optimal |
| setA-08 | 1.0000 | 68.1% | 1937 | Not optimal |
| setA-09 | 0.6111 | 17.3% | 3494 | Not optimal |
| setA-10 to 12 | – | – | – | Out of memory (killed at 26–27 GB while building) |
| setA-13 to 20 | – | – | – | Not run (37 to 958 million binary variables) |

**Step 9, two periods (3600 s):**

| Instance | Joint model | Two-stage method |
|---|---|---|
| setA-01 | 0.9294, optimal (20 s) | 0.9294, proven optimal |
| setA-02 | 0.9031, optimal (1736 s) | 0.9031, same optimum |
| setA-03 | 0.9435, optimal (203 s) | 0.9435, proven optimal |
| setA-04 | 1.0000 (no improvement) | **0.6918** (gap 16.0%) |
| setA-05 | 1.0000 (no improvement) | **0.2050** (gap 43.6%) |
| setA-06 | out of memory | 1.0000 (limited by the period-0 solution) |
| setA-07 | out of memory | out of memory |
| setA-08 | out of memory | 1.0000 (limited by the period-0 solution) |
| setA-09 | out of memory | stopped after 100 min |

The joint-model column uses the presolve-on run for setA-01 to 03 and the presolve-off run for 04 and 05.

## Known issues and lessons

- **Model size.** The model has |D|·n(n−1) binary variables per period. About 9 million fit in 30 GB, so setA-10 and above cannot be built.
- **HiGHS presolve ignores the time limit** on the joint two-period model from setA-04 upward: it ran for hours. `src/model.jl` therefore switches presolve off for the joint model. That keeps the time limit, but the LP relaxation becomes much harder. The two-stage method avoids both problems.
- **Budget constraint construction.** The budget sum must be written as a generator (`sum(z[d, e] for d in D, e in allseg)`), not `sum(z)`. Summing the JuMP array directly is quadratic and took 2570 s to build on setA-04, against 6 s with the generator.
- The `p01` srpaths files for setA-01 to 05 come from the presolve-off rerun.

## Report

The full report, with the model, explanations, all tables and the discussion, is the group's LaTeX report *"Mathematical Modelling – EURO/ROADEF Challenge 2026: Keep The Flow!"*.

## References

- Challenge subject: `Challenge_Orange_ROADEF_2026_Subject.pdf` and `calculating_r.pdf` in the course repository.
- C. Filsfils et al., *Segment Routing Architecture*, RFC 8402, 2018.
- M. Lubin et al., *JuMP 1.0*, Mathematical Programming Computation, 2023.
- Q. Huangfu, J. A. J. Hall, *Parallelizing the dual revised simplex method* (HiGHS), Mathematical Programming Computation, 2018.
