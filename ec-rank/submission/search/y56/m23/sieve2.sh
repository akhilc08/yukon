#!/bin/bash
# m23 sieve loop (self-driving): 2D-sweep ten-root sieve (m23f MODE=3), rotating root spread M and pairing,
# fresh seeds every round; output out/s3_M*_P*_seed.txt, synced to S3 after each chunk. Cores 0-10.
W=/work/y56/m23; cd $W; OD=${OUTD:-out}; LG=${LOGF:-log.txt}; mkdir -p $OD
export OMP_NUM_THREADS=${NT:-11}
round=${ROUND0:-10}
while [ ! -f $W/STOP ]; do
  for cfg in ${MLIST:-100:3 120:3 140:3 170:3 120:0 200:3 140:4 240:3}; do
    M=${cfg%%:*}; NP=${cfg##*:}; seed=$((round*1000+M+NP))
    NB=$(( ${CHUNK:-40000000000} * 2 / (M*M) ))
    echo "$(date -u +%FT%T) start MODE=3 M=$M NPAIR=$NP seed=$seed NB=$NB" >> $LG
    MODE=3 NPAIR=$NP G=8 PB=47 GATE=10.5 C97=16 C199=27 C1000=56 T2=${T2:-82} PMAX=3000 M=$M NB=$NB SEED=$seed HMAX=17.88 \
      taskset -c ${CORES:-0-10} ./m23f > $OD/s3_M${M}_P${NP}_${seed}.txt 2>> $LG
    aws s3 sync --quiet $W/$OD s3://ecrank-y56-115135126148/results/m23/$OD/
    aws s3 cp --quiet $LG s3://ecrank-y56-115135126148/results/m23/$LG
    [ -f $W/STOP ] && break
  done
  round=$((round+1))
done
