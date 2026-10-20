#!/usr/bin/env bash
# Chạy toàn bộ project trong một cgroup giới hạn cứng 15 GB RAM (không swap).
# OOMPolicy=continue: nếu kernel OOM-kill một job thì chỉ job đó chết (bộ điều phối chạy lại nó),
# không dừng cả scope như mặc định của systemd.
#   ./run.sh                    # mặc định: 1 seed / instance
#   ORANGE_NSEEDS=10 ./run.sh   # như readme (10 lần chạy)
cd "$(dirname "$0")"
exec systemd-run --user --scope --collect --slice=orange.slice -p MemoryMax=15000000000 -p MemorySwapMax=0 -p OOMPolicy=continue --unit="orange-$(date +%s)" \
     julia --project=. run_all.jl
