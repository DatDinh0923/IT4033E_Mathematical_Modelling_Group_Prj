LIMIT_KB=$((26*1024*1024))   # 26 GB
while true; do
  ps -eo pid,rss,args | awk -v L=$LIMIT_KB '/julia/ && /run_experiments/ && $2 > L {print $1, $2, $0}' |
  while read pid rss rest; do
    echo "$(date '+%F %T') killed pid $pid at $((rss/1024/1024)) GB: $rest" >> logs/watchdog.log
    kill -9 $pid
  done
  sleep 5
done
