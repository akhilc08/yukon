#!/bin/bash
# driver.sh [ULIST]: resumable queue over u slices (cores 0-27). Per slice: full Hlow<HB sieve, rank top MAXR by S(30000).
# Slices with run/TAG/done (restored from S3) are skipped. Hits -> cand/*.ec + .verify -> S3.
W=/work/y56/z2ek; cd $W
UL=${1:-$W/ulist.txt}; HB=${HB:-17.6}; MAXR=${MAXR:-3000}; SB=${SB:-88}
S3=s3://ecrank-y56-115135126148/results/z2ek
mkdir -p $W/cand $W/logs
certify() {
  for h in $W/run/*/hits/*.txt; do
    [ -f "$h" ] || continue
    b=$(basename $h .txt); [ -f $W/cand/$b.ec ] && continue
    python3 $W/mkec.py $h $W/cand/$b.ec > /dev/null 2>&1
    $W/verify.sh $W/cand/$b.ec > $W/cand/$b.verify 2>&1
  done
  aws s3 sync $W/cand $S3/cand --quiet
}
while read TAG U M; do
  [ -z "$TAG" ] && continue
  [ -f $W/run/$TAG/done ] && continue
  echo "$(date -u +%FT%TZ) start $TAG u=$U m=$M" >> logs/driver.log
  ./pipe.sh $TAG $U $M 1000000 100000000 -1e9 -1e9 -1e9 $HB $SB $MAXR >> logs/pipe.log 2>&1
  python3 stats.py run/$TAG/rank run/$TAG/torank > run/$TAG/stats.txt 2>/dev/null
  echo "$(date -u +%FT%TZ) done $TAG $(cat run/$TAG/nsurv 2>/dev/null) surv; $(grep -E '^rank 1[3-9]' run/$TAG/stats.txt | tr '\n' ';')" >> logs/driver.log
  certify
  aws s3 sync $W/logs $S3/logs --quiet
done < $UL
echo "$(date -u +%FT%TZ) QUEUE DONE $UL" >> logs/driver.log; aws s3 sync $W/logs $S3/logs --quiet
