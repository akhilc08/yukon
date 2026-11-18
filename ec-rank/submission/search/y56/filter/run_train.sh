#!/bin/bash
# Train + evaluate + export + bench on the box:  bash run_train.sh <tag> [K ...]
# Uses only completed chunks (*.done) from /work/y56/filter/data; outputs /work/y56/filter/runs/<tag>.
set -e
TAG=$1; shift; KS=${@:-25 50 100 168}
DIR=$(cd $(dirname $0); pwd); PY=/work/mlenv/bin/python; DATA=/work/y56/filter/data; OUT=/work/y56/filter/runs/$TAG
mkdir -p $OUT/snap; cd $DATA
for d in *.done; do ln -sf $DATA/${d%.done} $OUT/snap/; done
cd $DIR
export OMP_NUM_THREADS=6
taskset -c 22-27 $PY train.py $OUT/snap $OUT $KS > $OUT/train.log 2>&1
$PY report.py $OUT > $OUT/report.txt
$PY export.py $OUT/final_models.pkl $OUT/model $KS
$PY known_check.py /work/y56/filter/known_ap.txt $OUT/model/models.json > $OUT/known_check.txt 2>&1 || true
$PY -c "
import sys; sys.path.insert(0,'$DIR'); import data as D, numpy as np
ds=D.load('$OUT/snap'); i=np.random.default_rng(0).permutation(len(ds.r))[:50000]
np.savetxt('$OUT/ap_sample.txt', np.where(ds.bad,0,ds.ap)[i], fmt='%d')"
for K in $KS; do gcc -O2 -mcpu=native -I$OUT/model -DBK=$K bench.c -o $OUT/bench$K -lm && taskset -c 22 $OUT/bench$K < $OUT/ap_sample.txt; done > $OUT/bench.txt
rm -f $OUT/ap_sample.txt
aws s3 sync $OUT s3://ecrank-y56-115135126148/results/filter/runs/$TAG --exclude 'snap/*' --exclude 'bench*[0-9]' --quiet
echo DONE
