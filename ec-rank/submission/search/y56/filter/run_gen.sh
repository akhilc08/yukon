#!/bin/bash
# Dataset generation on the box: bash run_gen.sh <outdir> <cores>   (cores e.g. 22-27)
# Each job = (family, band, seed) chunk; 6 in parallel pinned to the given cores.
OUT=${1:-/work/y56/filter/data}; CORES=${2:-22-27}; NCH=${NCH:-6}; NPER=${NPER:-10000}; NMES=${NMES:-2500}
GP=/work/ec-rank/tools/pari/bin/gp; DIR=$(cd $(dirname $0); pwd); mkdir -p $OUT
jobs=()
for s in $(seq 1 $NCH); do
  for f in mestre; do for b in 1 2; do jobs+=("$f $b $s $NMES"); done; done
  for f in Z6 Z2Z2 rand Z5 Z2Z4 Z7 Z8 Z2Z6; do for b in 1 2; do jobs+=("$f $b $s $NPER"); done; done
done
printf '%s\n' "${jobs[@]}" | taskset -c $CORES xargs -P 6 -L 1 bash -c '
f=$0; b=$1; s=$2; n=$3; o='$OUT'/$f.$b.$s.txt; [ -f $o.done ] && exit 0
echo "GENDIR=\"'$DIR'\";FAM=\"$f\";BAND=$b;N=$n;SEED=$((1000*s+b));OUT=\"$o\";read(\"'$DIR'/gen_feats.gp\");go();" | '$GP' -q -s 1G > $o.log 2>&1 && touch $o.done'
