#!/bin/bash
# runA.sh: sieve one family with psieve, then ellrank survivors sorted by S, one gp per core (28-31).
# usage: runA.sh TAG FAM QMAX QFULL HMAX T1 T2 T3 NRANK [W]
set -e
TAG=$1; FAM=$2; QM=$3; QF=$4; HM=$5; T1=$6; T2=$7; T3=$8; NR=$9; W=${10:-60}
GP="/work/ec-rank/tools/pari/bin/gp -q -s 1G -D parisizemax=4G"
cd /work/y56/tors; D=$PWD/$TAG; mkdir -p $D
OMP_NUM_THREADS=4 taskset -c 28-31 ./psieve $FAM $QM $QF $HM $T1 $T2 $T3 1 $W > $D/sieve.txt 2> $D/sieve.log
sort -k5 -gr $D/sieve.txt | head -n $NR | awk -v F=$FAM '{print "["F","$1","$2"]"}' > $D/rin.gp
split -n r/4 -d $D/rin.gp $D/rin_
i=28
for c in $D/rin_0?; do (echo "INF=\"$c\"; OUTF=\"$c.out\"; read(\"rankfam.gp\");" | taskset -c $i $GP > $c.log 2>&1) & i=$((i+1)); done
wait
echo DONE >> $D/status
