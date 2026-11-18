#!/bin/bash
# S(X) up to 1e5 for the first NMAX curves of chunks 1-2 of selected families, then the per-family refine eval.
# bash run_s5.sh   (cores 22-27; outputs /work/y56/filter/data_s5, /work/y56/filter/runs/refine)
DIR=$(cd $(dirname $0); pwd); GP=/work/ec-rank/tools/pari/bin/gp; DATA=/work/y56/filter/data; S5=/work/y56/filter/data_s5
NMAX=${NMAX:-7500}; FAMS=${FAMS:-"Z6 Z2Z2 rand Z2Z4"}; mkdir -p $S5
for f in $FAMS; do for b in 1 2; do for s in 1 2; do echo "$f.$b.$s"; done; done; done | taskset -c 22-27 xargs -P 6 -I{} bash -c '
[ -f '$S5'/{}.s5.done ] && exit 0
echo "IN=\"'$DATA'/{}.txt\";OUT=\"'$S5'/{}.s5\";NMAX='$NMAX';read(\"'$DIR'/add_s.gp\");go();" | '$GP' -q -s 256M && touch '$S5'/{}.s5.done'
mkdir -p /work/y56/filter/runs/refine/snap
for f in $FAMS; do for b in 1 2; do for s in 1 2; do ln -sf $DATA/$f.$b.$s.txt /work/y56/filter/runs/refine/snap/; done; done; done
cd $DIR; OMP_NUM_THREADS=6 taskset -c 22-27 /work/mlenv/bin/python train_refine.py /work/y56/filter/runs/refine/snap $S5 /work/y56/filter/runs/refine > /work/y56/filter/runs/refine/report.txt 2>&1
aws s3 sync /work/y56/filter/runs/refine s3://ecrank-y56-115135126148/results/filter/runs/refine --exclude 'snap/*' --quiet
aws s3 sync $S5 s3://ecrank-y56-115135126148/results/filter/data_s5 --quiet
echo DONE
