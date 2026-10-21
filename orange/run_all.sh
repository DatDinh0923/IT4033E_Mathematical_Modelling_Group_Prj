# readme Step 6: period 0, 1800 s limit
for i in 01 02 03; do
  julia --project=. run_experiments.jl p0 1800 10 setA-$i > logs/p0-setA-$i.log 2>&1
done
for i in 04 05 06 07 08 09 10 11 12; do
  julia --project=. run_experiments.jl p0 1800 1 setA-$i > logs/p0-setA-$i.log 2>&1
done

# readme Step 9: periods 0+1 with budget, 3600 s limit
for i in 01 02 03 04 05 06 07 08 09 10 11 12; do
  julia --project=. run_experiments.jl p01 3600 1 setA-$i > logs/p01-setA-$i.log 2>&1
done