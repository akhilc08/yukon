#!/bin/bash
# twelve-root family sieve: f list from fsol (R), then m12f over v in [V0,V], u <= UMAX.  cores via CORES.
W=/work/y56/m23/r12; cd $W
gcc -O3 -mcpu=neoverse-v1 -fopenmp -Wno-empty-body -o m12f m12f.c -lm
gcc -O3 -o fsol2 fsol2.c
[ -f f${R:-100}.txt ] || ./fsol2 ${R:-100} > f${R:-100}.txt
wc -l f${R:-100}.txt >> log12.txt
for v in $(seq ${V0:-1} ${V:-6}); do
  echo "$(date -u +%FT%T) v=$v start" >> log12.txt
  OMP_NUM_THREADS=${NT:-5} V0=$v V=$v UMAX=${UMAX:-4000} PB=47 GATE=10.5 C97=16 C199=27 C1000=54 T2=${T2:-78} PMAX=3000 \
    taskset -c ${CORES:-6-10} ./m12f < f${R:-100}.txt > s12_v$v.txt 2>> log12.txt
  aws s3 sync --quiet $W s3://ecrank-y56-115135126148/results/m23/r12/ --exclude "m12f" --exclude "fsol"
  [ -f /work/y56/m23/STOP ] && break
done
