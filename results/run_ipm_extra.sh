#!/usr/bin/env bash
# Thử nghiệm bổ sung: setA-06 và setA-08 (mô hình rút gọn, cùng K) giải LP gốc bằng IPM thay vì simplex.
cd "$(dirname "$0")/.."
run() { ORANGE_LP=ipm ORANGE_TAG=_ipm systemd-run --user --scope --collect --slice=orange.slice -p OOMPolicy=continue \
        --unit="orange-ipm-$1-$2" julia --project=. run_one.jl "$1" "$2" 1 restricted "$3" "$4"; }
run setA-06 6 1800 23 & run setA-08 6 1800 30 & wait
run setA-06 9 3600 23 & run setA-08 9 3600 30 & wait
