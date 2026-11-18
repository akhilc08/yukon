#!/bin/bash
# pipe.sh FAMK QMIN QMAX QSTEP PMAX T1 T2 T3 HB SB TAG
# chunked: C sieve (10 threads on cores 12-21) -> gp rank stage (10 procs) -> sync to S3
FAMK=$1; QMIN=$2; QMAX=$3; QSTEP=$4; PMAX=$5; T1=$6; T2=$7; T3=$8; HB=$9; SB=${10}; TAG=${11}
W=/work/y56/z2r20; R=$W/run_$TAG; mkdir -p $R
cd $W
q=$QMIN
while [ $q -le $QMAX ]; do
  q2=$((q+QSTEP-1)); [ $q2 -gt $QMAX ] && q2=$QMAX
  C=$R/c_${q}_${q2}
  if [ ! -f $C.done ]; then
    TABLES=$W/tab_X${FAMK}_30000.bin OMP_NUM_THREADS=10 taskset -c 12-21 ./ksieve X$FAMK.fam 27.55 $q $q2 $PMAX $T1 $T2 $T3 200 1000 30000 > $C.surv 2> $C.err
    ./rankpar.sh $FAMK $C.surv $HB $SB ${TAG}_${q}
    cat rk_${TAG}_${q}/shard_*.out > $C.rank 2>/dev/null
    mkdir -p $R/hits; cp rk_${TAG}_${q}/hits/* $R/hits/ 2>/dev/null
    rm -rf rk_${TAG}_${q}
    touch $C.done
    aws s3 sync $R s3://ecrank-y56-115135126148/results/z2r20/run_$TAG --quiet
  fi
  q=$((q2+1))
done
echo finished > $R/FINISHED
aws s3 sync $R s3://ecrank-y56-115135126148/results/z2r20/run_$TAG --quiet
