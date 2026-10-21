LIMIT_KB=$((26*1024*1024)); LIMIT_S=6000
while true; do
  ps -eo pid,etimes,rss,args | awk -v L=$LIMIT_KB -v S=$LIMIT_S '/julia/ && /two_stage/ && !/awk/ && ($3 > L || $2 > S) {print $1, $2, $3}' |
  while read pid secs rss; do
    echo "$(date '+%F %T') killed pid $pid after $((secs/60)) min at $((rss/1024/1024)) GB: $(ps -o args= -p $pid)" >> logs/guard_ts.log
    kill -9 $pid
  done
  sleep 10
done
