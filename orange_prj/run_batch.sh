#!/usr/bin/env bash
# Experiments of Steps 6 and 9, one Julia process per run.
#
#   MODEL=t0  LP=ipm ./run_batch.sh 2 3 4      # Step 6: period 0, 1800 s
#   MODEL=t01 LP=ipm ./run_batch.sh 1 2        # Step 9: periods 0 and 1, 3600 s
#
# Each instance is run once (seed 0). If that run ends Optimal within the time
# limit, it is repeated with seeds 1..RUNS-1, so the tables can report the
# average of RUNS runs. A run that reaches the time limit is not repeated.
#
# Environment: MODEL (t0), LP (choose), TIME_LIMIT (1800 for t0, 3600 for t01),
# RUNS (10), OUT (results/<MODEL>_results.csv), MEM_MAX (e.g. 20G: memory cap
# per run; a run killed by it is listed in results/failures.csv), PRESOLVE
# (on|off, HiGHS presolve; default: HiGHS's choice), CODE_DIR (directory
# holding run.jl, default: this directory).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CODE_DIR="${CODE_DIR:-$HERE}"
MODEL="${MODEL:-t0}"
LP="${LP:-choose}"
if [ "$MODEL" = t0 ]; then TIME_LIMIT="${TIME_LIMIT:-1800}"; else TIME_LIMIT="${TIME_LIMIT:-3600}"; fi
RUNS="${RUNS:-10}"
OUT="${OUT:-$HERE/results/${MODEL}_results.csv}"
mkdir -p "$HERE/results/logs"

run_once() {  # instance seed; succeeds when the run ends Optimal
    local name log start code wall
    name=$(printf "setA-%02d" "$1")
    log="$HERE/results/logs/$name-$MODEL-$LP${PRESOLVE:+-presolve-$PRESOLVE}-seed$2.log"
    local cmd=(julia --project="$CODE_DIR" "$CODE_DIR/run.jl" --model "$MODEL"
               --time-limit "$TIME_LIMIT" --mip-lp-solver "$LP" ${PRESOLVE:+--presolve "$PRESOLVE"} --seed "$2"
               --data-dir "$HERE/../setA" --results-dir "$HERE/results" --out "$OUT" "$1")
    if [ -n "${MEM_MAX:-}" ]; then
        cmd=(systemd-run --user --scope --quiet -p MemoryMax="$MEM_MAX" -p MemorySwapMax=0 "${cmd[@]}")
    fi

    start=$(date +%s)
    "${cmd[@]}" > "$log" 2>&1
    code=$?
    wall=$(( $(date +%s) - start ))
    echo "$name $MODEL lp=$LP${PRESOLVE:+ presolve=$PRESOLVE} seed=$2: exit $code, wall $wall s"

    if [ $code -ne 0 ]; then
        [ -f "$HERE/results/failures.csv" ] ||
            echo "instance,model,mip_lp_solver,seed,exit_code,wall_s,mem_max,date" > "$HERE/results/failures.csv"
        echo "$name,$MODEL,$LP${PRESOLVE:+ presolve=$PRESOLVE},$2,$code,$wall,${MEM_MAX:-none},$(date -Is)" >> "$HERE/results/failures.csv"
        return 1
    fi
    grep "^\[$name\] " "$log" | tail -1
    grep -q "^\[$name\] Optimal" "$log"
}

batch_start=$(date +%s)
for k in "$@"; do
    if run_once "$k" 0; then
        for seed in $(seq 1 $((RUNS - 1))); do run_once "$k" "$seed"; done
    fi
done
echo "batch MODEL=$MODEL LP=$LP${PRESOLVE:+ PRESOLVE=$PRESOLVE} $*: total wall $(( $(date +%s) - batch_start )) s"
