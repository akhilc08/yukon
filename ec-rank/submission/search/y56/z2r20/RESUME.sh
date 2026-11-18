#!/bin/bash
# RESUME.sh: one-command restart of the z2r20 search on a fresh box (run as root or ubuntu).
# Restores code + completed-chunk markers from S3, recompiles for this CPU, relaunches driver2.sh in tmux.
set -u
W=/work/y56/z2r20; S3=s3://ecrank-y56-115135126148
mkdir -p $W && cd $W
aws s3 cp $S3/code/z2r20/z2r20_all.tgz . --quiet && tar xzf z2r20_all.tgz 2>/dev/null
# completed chunks, ranks, hits, candidates (results are small); tables are rebuilt automatically
for d in run_x552a run_x808a run_lowh run_lowh2 run_bc552 cand; do aws s3 sync $S3/results/z2r20/$d $W/$d --quiet; done
gcc -O3 -mcpu=native -fopenmp -o ksieve ksieve.c -lm
gcc -O3 -mcpu=native -o nagao nagao.c -lm
chmod +x *.sh
rm -rf rk_* ; rm -f run_*/FINISHED
# drop partial (not .done) chunk outputs so they are redone cleanly
for f in run_x*/c_*.surv; do [ -f "${f%.surv}.done" ] || rm -f "$f" "${f%.surv}.err" "${f%.surv}.rank"; done
chown -R ubuntu:ubuntu /work/y56
if ! sudo -u ubuntu tmux has-session -t z2r20drv 2>/dev/null; then
  sudo -u ubuntu tmux new-session -d -s z2r20drv "cd $W && ./driver2.sh > logs_driver2.txt 2>&1"
fi
echo "z2r20 driver running: $(ls $W/run_x552a/*.done 2>/dev/null | wc -l) X552 + $(ls $W/run_x808a/*.done 2>/dev/null | wc -l) X808 chunks already done"
