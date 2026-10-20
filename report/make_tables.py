#!/usr/bin/env python3
"""Generate the LaTeX tables and the bar-chart data of the report from ../results.

    python3 make_tables.py        (run from orange/report)

Every number in tables/*.tex comes from results/instances/*.json and results/raw/*/*.json,
so the report can be regenerated after a new experimental run.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
OUT = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)
NAMES = [f"setA-{i:02d}" for i in range(1, 21)]


def load(path):
    with open(path) as f:
        return json.load(f)


inst = {nm: load(os.path.join(RES, "instances", nm + ".json")) for nm in NAMES}
raw = {(nm, s): load(os.path.join(RES, "raw", f"step{s}", nm + "_s1.json")) for nm in NAMES for s in (6, 9)}


def f4(x):
    return f"{x:.4f}"


def model_label(r):
    return "exact" if r["mode"] == "exact" else f"$K{{=}}{r['K']}$"


def certified(nm, step):
    """Optimal for the original problem: exact+Optimal, objective equal to the MCF lower bound,
    or (step 9) objective equal to the proven step-6 optimum (step 6 is a relaxation of step 9)."""
    r = raw[(nm, step)]
    obj = r["best_obj"]
    if r["mode"] == "exact" and r["status"] == "Optimal":
        return True
    lb = inst[nm]["mcf_lb_t0" if step == 6 else "mcf_lb_t01"]
    if lb is not None and abs(obj - lb) <= 1e-6 * max(1.0, obj):
        return True
    if step == 9:
        r6 = raw[(nm, 6)]
        if r6["mode"] == "exact" and r6["status"] == "Optimal" and abs(obj - r6["best_obj"]) <= 1e-6:
            return True
    return False


# ---------------------------------------------------------------- instances
rows = []
for nm in NAMES:
    I = inst[nm]
    rows.append(" & ".join([
        nm[-2:], str(I["n"]), str(I["arcs"]), str(I["demands"]), str(I["budget_t1"]),
        str(I["down_links_t1"][0]),
        f4(I["baseline_mlu_t0"]), f4(I["baseline_mlu_t1"]),
        f4(I["mcf_lb_t0"]), f4(I["mcf_lb_t01"]),
        f"{I['x_vars_step9'] / 1e6:.2f}",
    ]) + r" \\")
with open(os.path.join(OUT, "instances.tex"), "w") as f:
    f.write("\n".join(rows) + "\n")


# ---------------------------------------------------------------- results
def result_rows(step):
    out = []
    for nm in NAMES:
        r = raw[(nm, step)]
        I = inst[nm]
        status = "Opt." if r["status"] == "Optimal" else "Not opt."
        gap = "--" if r["status"] == "Optimal" else f"{100 * r['gap']:.1f}"
        cpu = f"{r['cpu_time']:.1f}"
        if r["cpu_time"] > r["time_limit"] + 5:
            cpu += r"$^\ddagger$"
        obj = f4(r["best_obj"]) + (r"$^\star$" if certified(nm, step) else "")
        lb = I["mcf_lb_t0" if step == 6 else "mcf_lb_t01"]
        cells = [nm[-2:], str(I["n"]), str(I["arcs"]), str(I["demands"]), model_label(r),
                 status, obj, gap, cpu, str(r["nodes"]), f4(r["baseline_mlu"]), f4(lb)]
        if step == 9:
            cells.append(f"{r['changes'][0]}/{r['budget']}")
        out.append(" & ".join(cells) + r" \\")
    return out


for step in (6, 9):
    with open(os.path.join(OUT, f"step{step}.tex"), "w") as f:
        f.write("\n".join(result_rows(step)) + "\n")

# ---------------------------------------------------------------- IPM supplement
rows = []
for nm in ("setA-06", "setA-08"):
    for step in (6, 9):
        rs = raw[(nm, step)]
        ri = load(os.path.join(RES, "raw", f"step{step}_ipm", nm + "_s1.json"))

        def cell(r):
            b = r["best_bound"]
            return f4(r["best_obj"]), ("--" if b is None else f4(max(b, 0.0))), f"{r['cpu_time']:.0f}"

        rows.append(" & ".join([nm[-2:], str(step), f"$K{{=}}{rs['K']}$", *cell(rs), *cell(ri)]) + r" \\")
with open(os.path.join(OUT, "ipm.tex"), "w") as f:
    f.write("\n".join(rows) + "\n")

# ---------------------------------------------------------------- bar chart (step 9)
lines = []
for k, nm in enumerate(NAMES):
    r = raw[(nm, 9)]
    lines.append(f"{k}/{nm[-2:]}/{r['baseline_mlu']:.4f}/{r['best_obj']:.4f}")
with open(os.path.join(OUT, "bars9.tex"), "w") as f:
    f.write(r"\def\barsdata{" + ",".join(lines) + "}\n")

# ---------------------------------------------------------------- summary numbers used in the text
summary = {}
for step in (6, 9):
    rs = [raw[(nm, step)] for nm in NAMES]
    summary[step] = dict(
        optimal=sum(r["status"] == "Optimal" for r in rs),
        certified=[nm[-2:] for nm in NAMES if certified(nm, step)],
        improved=sum(r["improvement"] > 1e-6 for r in rs),
        mean_impr=100 * sum(r["improvement"] for r in rs) / len(rs),
        cpu_h=sum(r["cpu_time"] for r in rs) / 3600,
    )
print(json.dumps(summary, indent=1))
