#!/bin/bash
# runR.sh: re-rank an existing sieve output. usage: runR.sh SRCDIR TAG "GPSETUP" HAPPROXMAX NR HCUT
set -e
SRC=$1; TAG=$2; SETUP=$3; HA=$4; NR=$5; HCUT=$6
GP="/work/ec-rank/tools/pari/bin/gp -q -s 1G -D parisizemax=4G"
cd /work/y56/tors; D=$PWD/$TAG; mkdir -p $D
awk -v h=$HA '$3<h' $SRC/sieve.txt | sort -k5 -gr | head -n $NR | awk '{print "["$1","$2"]"}' > $D/rin.gp
split -n r/4 -d $D/rin.gp $D/rin_
i=28
for c in $D/rin_0?; do (echo "$SETUP INF=\"$c\"; OUTF=\"$c.out\"; HCUT=$HCUT; MINR=9; read(\"rankg.gp\");" | taskset -c $i $GP > $c.log 2>&1) & i=$((i+1)); done
wait
echo DONE >> $D/status
