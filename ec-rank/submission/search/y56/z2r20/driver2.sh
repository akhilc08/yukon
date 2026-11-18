#!/bin/bash
# driver2.sh: resumable queue on cores 12-21 (<=10 threads). Chunks with a .done marker (restored from S3) are skipped.
# X552 q 100001..130000 (step 5000), then X808 q 1260001..2000000 (step 20000). ellrank cutoff S(1e5) >= 100.
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
for ((q=100001; q<=130000; q+=5000)); do
  ./pipe.sh 552 $q $((q+4999)) 5000 400000 17 35 81 27.55 100 x552a >> logs/pipe552.log 2>&1; certify
done
for ((q=1260001; q<=2000000; q+=20000)); do
  ./pipe.sh 808 $q $((q+19999)) 20000 400000 17 35 81 27.55 100 x808a >> logs/pipe808.log 2>&1; certify
done
echo ALLDONE > $W/logs/DRIVER2_DONE; aws s3 sync $W/logs $S3/logs --quiet
