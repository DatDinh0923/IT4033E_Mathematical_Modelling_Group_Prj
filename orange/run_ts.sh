for i in 01 02 03 04 05 06 07 08 09; do
     julia --project=. two_stage.jl 3600 1 setA-$i > logs/p01ts-setA-$i.log 2>&1
   done
