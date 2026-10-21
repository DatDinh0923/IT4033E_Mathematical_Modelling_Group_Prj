while true; do
  ps -eo pid,etimes,args | awk '/julia/ && /run_experiments/ && !/awk/ && $2 > 4800 {print $1, $2}' |
  while read pid secs; do
    echo "$(date '+%F %T') killed pid $pid after $((secs/60)) min: $(ps -o args= -p $pid)" >> logs/timeguard.log
    kill -9 $pid
  done
  sleep 30
done
