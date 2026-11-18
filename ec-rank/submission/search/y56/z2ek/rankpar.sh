#!/bin/bash
# rankpar.sh U M SURV HB SB TAG NCORE FIRSTCORE : shard survivors (already sorted by S desc) over cores, one gp per core
U=$1; M=$2; SURV=$3; HB=$4; SB=$5; TAG=$6; NC=${7:-28}; C0=${8:-0}
W=/work/y56/z2ek; D=$W/rk/$TAG; mkdir -p $D/hits
split -n r/$NC -d -a 3 $SURV $D/shard_
i=$C0
for f in $D/shard_???; do
  U=$U M=$M IN=$f OUT=$f.out HB=$HB SB=$SB HITS=$D/hits TAG=$TAG nohup taskset -c $i /work/ec-rank/tools/pari/bin/gp -q -s 2000000000 $W/rankek.gp < /dev/null > $f.log 2>&1 &
  i=$((i+1))
done
wait
