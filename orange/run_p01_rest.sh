for i in 05 06 07 08 09 10 11 12; do
     julia --project=. run_experiments.jl p01 3600 1 setA-$i > logs/p01-setA-$i.log 2>&1
   done
