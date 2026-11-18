#!/bin/bash
# driver.sh: self-driving queue on cores 12-21 (<=10 threads at any time).
# 1) low-height pass (cells 2/16, 2/17); 2) interleaved X808/X552 chunks (sieve -> rank); 3) certify hits; sync to S3.
W=/work/y56/z2r20; cd $W
S3=s3://ecrank-y56-115135126148/results/z2r20
mkdir -p $W/cand $W/logs
certify() {
  for h in $W/run_*/hits/*.txt; do
    [ -f "$h" ] || continue
    b=$(basename $h .txt); [ -f $W/cand/$b.ec ] && continue
    python3 $W/mkec.py $h $W/cand/$b.ec > /dev/null 2>&1
    $W/verify.sh $W/cand/$b.ec > $W/cand/$b.verify 2>&1
  done
  aws s3 sync $W/cand $S3/cand --quiet
}
runqueue() {  # run any queued jobs (each must keep to cores 12-21, <=10 threads), then mark done
  mkdir -p $W/queue/done
  for j in $W/queue/*.sh; do [ -f "$j" ] || continue; bash $j > $W/logs/$(basename $j).log 2>&1; mv $j $W/queue/done/; certify; done
}
lowh() {  # FAM QMAX
  F=$1; QM=$2; R=$W/run_lowh; mkdir -p $R
  [ -f $R/X$F.done ] && return
  TABLES=$W/tab_X${F}_30000.bin OMP_NUM_THREADS=10 taskset -c 12-21 ./ksieve X$F.fam 18.92 1 $QM 3000000 -1e9 -1e9 95 200 1000 30000 > $R/X$F.surv 2> $R/X$F.err
  ./rankpar.sh $F $R/X$F.surv 18.92 -1e9 lowh$F
  cat rk_lowh$F/shard_*.out > $R/X$F.rank; mkdir -p $R/hits; cp rk_lowh$F/hits/* $R/hits/ 2>/dev/null; rm -rf rk_lowh$F
  touch $R/X$F.done; aws s3 sync $R $S3/run_lowh --quiet
}
runqueue
lowh 552 20000; certify
lowh 808 2000000; certify
# interleaved main queue: X808 q-chunks of 20000 up to 2e6, X552 q-chunks of 5000 up to 1e5
i=0
while :; do
  runqueue
  q8=$((1 + i*20000)); q5=$((1 + i*5000)); did=0
  if [ $q8 -le 2000000 ]; then ./pipe.sh 808 $q8 $((q8+19999)) 20000 400000 17 35 81 27.55 105 x808a >> logs/pipe808.log 2>&1; did=1; certify; fi
  if [ $q5 -le 100000 ]; then ./pipe.sh 552 $q5 $((q5+4999)) 5000 400000 17 35 81 27.55 105 x552a >> logs/pipe552.log 2>&1; did=1; certify; fi
  [ $did -eq 0 ] && break
  i=$((i+1))
done
echo ALLDONE > $W/logs/DRIVER_DONE; aws s3 sync $W/logs $S3/logs --quiet
