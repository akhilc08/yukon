#!/bin/bash
# driver.sh: self-driving queue for y56/tors on cores 28-31 (<=4 threads total, one job at a time).
# Runs each line of jobs.txt (once; done keys in jobs.done), then: hits.py -> mkcand.sh for every new hit ->
# ec_verify on the candidate -> sync results (no sieve.txt) and cands to s3://ecrank-y56-115135126148/results/tors/.
# Append lines to jobs.txt at any time; the driver polls every 60 s when idle.
cd /work/y56/tors
FAMD_E1="Elkies disc -163 K3 Z/4 rank-4 family E1 (fam/z4_elkies163_E1.gp)"
post() {
  python3 hits.py > hits.new 2>>driver.log
  while IFS=$'\t' read key cell h tg job p q iso name; do
    [ -z "$key" ] && continue
    echo "$(date -u +%FT%TZ) HIT $cell H=$h target=$tg $job t=$p/$q iso=$iso" >> driver.log
    fd=$(grep -h "^#FAMD" $job/famd 2>/dev/null | head -1 | cut -c7-); [ -z "$fd" ] && fd="$job"
    ./mkcand.sh $job $p $q $name "$fd" $iso >> driver.log 2>&1 && \
      (cd /work/ec-rank && timeout 600 python3 -m ec_verify /work/y56/tors/cands/$name.ec > /work/y56/tors/cands/$name.verify.json 2>&1)
    echo "$key" >> hits.seen
    aws s3 sync /work/y56/tors/cands s3://ecrank-y56-115135126148/results/tors/cands --quiet || true
    aws s3 cp hits.seen s3://ecrank-y56-115135126148/results/tors/hits.seen --quiet || true
  done < hits.new
  aws s3 sync /work/y56/tors s3://ecrank-y56-115135126148/results/tors --exclude '*sieve.txt' --exclude '*.tgz' --exclude 'rin_0?' --quiet || true
}
while true; do
  ran=0
  while IFS= read -r line; do
    [ -z "$line" ] && continue; case "$line" in \#*) continue;; esac
    key=$(echo "$line" | md5sum | cut -c1-12)
    grep -q "^$key" jobs.done 2>/dev/null && continue
    echo "$(date -u +%FT%TZ) START $key $line" >> driver.log
    bash -c "$line" >> driver.log 2>&1
    echo "$key $(date -u +%FT%TZ) $line" >> jobs.done
    aws s3 cp jobs.done s3://ecrank-y56-115135126148/results/tors/jobs.done --quiet || true
    echo "$(date -u +%FT%TZ) END $key" >> driver.log
    post; ran=1; break
  done < jobs.txt
  [ $ran = 0 ] && { post; sleep 60; }
done
