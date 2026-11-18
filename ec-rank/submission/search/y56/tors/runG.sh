#!/bin/bash
# runG.sh: gsieve a spec family, then exact-H filter + seeded ellrank of the top NR survivors, one gp per core (28-31).
# usage: runG.sh TAG SPEC "GPSETUP" QM QF HAPPROX T1 T2 T3 NR HCUT [PMUL P0 W SYM]
#   GPSETUP: gp code defining FF (t->ainvs) and SX (t->section x list), e.g. 'read("../families/x.gp"); FF=fam; SX=secx;'
set -e
TAG=$1; SPEC=$2; SETUP=$3; QM=$4; QF=$5; HA=$6; T1=$7; T2=$8; T3=$9; NR=${10}; HCUT=${11}; PMUL=${12:-20}; P0=${13:-500}; W=${14:-60}; SYM=${15:-0}
GP="/work/ec-rank/tools/pari/bin/gp -q -s 1G -D parisizemax=4G"
cd /work/y56/tors; D=$PWD/$TAG; mkdir -p $D
# resumable: if the ranked list rin.gp already exists (restored from s3), skip the sieve; rankg.gp resumes from .prog
if [ ! -s $D/rin.gp ]; then
  OMP_NUM_THREADS=4 taskset -c 28-31 ./gsieve $SPEC $QM $QF $HA $T1 $T2 $T3 $PMUL $P0 $W $SYM > $D/sieve.txt 2> $D/sieve.log
  sort -k5 -gr $D/sieve.txt | head -n $NR | awk '{print "["$1","$2"]"}' > $D/rin.gp.tmp && mv $D/rin.gp.tmp $D/rin.gp
  aws s3 sync $D s3://ecrank-y56-115135126148/results/tors/$(basename $D) --exclude 'sieve.txt' --exclude 'rin_0?' --quiet || true
fi
split -n r/4 -d $D/rin.gp $D/rin_
i=28
for c in $D/rin_0?; do (echo "$SETUP INF=\"$c\"; OUTF=\"$c.out\"; HCUT=$HCUT; MINR=9; read(\"rankg.gp\");" | taskset -c $i $GP > $c.log 2>&1) & i=$((i+1)); done
wait
echo DONE >> $D/status
