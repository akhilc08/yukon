#!/bin/bash
# pipe.sh TAG U M QMAX PMAX T1 T2 T3 HB SB [MAXRANK]: build fam, C sieve (cores 0-27), sort by S, rank top MAXRANK in gp (28 procs), sync.
TAG=$1; U=$2; M=$3; QMAX=$4; PMAX=$5; T1=$6; T2=$7; T3=$8; HB=$9; SB=${10}; MAXR=${11:-100000}
W=/work/y56/z2ek; R=$W/run/$TAG; mkdir -p $R; cd $W
S3=s3://ecrank-y56-115135126148/results/z2ek
[ -f $R/done ] && exit 0
[ -s fam/$TAG.fam ] || { mkdir -p fam; U=$U M=$M FN=fam/$TAG.fam /work/ec-rank/tools/pari/bin/gp -q -s 500000000 mkfam.gp > $R/mkfam.log 2>&1; }
grep -q '^SP' fam/$TAG.fam 2>/dev/null || { echo "bad fam $TAG" > $R/done; exit 0; }
if [ ! -f $R/sieve.done ]; then
  mkdir -p tab
  SH2=${SH2:-1.5} TABLES=$W/tab/$TAG.bin OMP_NUM_THREADS=28 taskset -c 0-27 ./ksieve fam/$TAG.fam $HB 1 $QMAX $PMAX $T1 $T2 $T3 200 1000 30000 > $R/surv 2> $R/sieve.err
  sort -k6,6gr $R/surv > $R/surv.sorted; head -n $MAXR $R/surv.sorted > $R/torank; touch $R/sieve.done
fi
[ -f $R/surv ] && wc -l < $R/surv > $R/nsurv; [ -f $R/surv.sorted ] && head -20000 $R/surv.sorted | gzip > $R/surv.top.gz; rm -f $R/surv $R/surv.sorted
./rankpar.sh $U $M $R/torank $HB $SB $TAG 28 0
cat rk/$TAG/shard_*.out > $R/rank 2>/dev/null
mkdir -p $R/hits; cp rk/$TAG/hits/* $R/hits/ 2>/dev/null
rm -rf rk/$TAG; rm -f tab/$TAG.bin
touch $R/done
aws s3 sync $R $S3/run/$TAG --quiet
