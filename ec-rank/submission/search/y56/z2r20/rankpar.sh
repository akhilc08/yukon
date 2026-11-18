#!/bin/bash
# rankpar.sh FAMK SURVIVORS HB SB TAG : shard survivors over cores 12-21, one gp per core
FAMK=$1; SURV=$2; HB=$3; SB=$4; TAG=$5
D=/work/y56/z2r20/rk_$TAG; mkdir -p $D/hits
split -n r/10 -d $SURV $D/shard_
i=12
for f in $D/shard_*; do
  FAMK=$FAMK IN=$f OUT=$f.out HB=$HB SB=$SB HITS=$D/hits nohup taskset -c $i /work/ec-rank/tools/pari/bin/gp -q -s 2000000000 /work/y56/z2r20/rankrun.gp < /dev/null > $f.log 2>&1 &
  i=$((i+1))
done
wait
