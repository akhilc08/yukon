#!/bin/bash
# RESUME.sh: restart the z2ek search on a fresh/rebooted box (run as root). Restores code + finished slices from S3,
# recompiles, relaunches driver.sh in tmux session z2ek-drv (cores 0-27).
set -u
W=/work/y56/z2ek; S3=s3://ecrank-y56-115135126148
mkdir -p $W && cd $W
aws s3 cp $S3/code/z2ek/z2ek.tgz . --quiet && tar xzf z2ek.tgz
for d in run cand logs; do aws s3 sync $S3/results/z2ek/$d $W/$d --quiet; done
gcc -O3 -mcpu=native -fopenmp -o ksieve ksieve.c -lm 2>/dev/null
chmod +x *.sh; rm -rf rk tab
# drop partial slices (no done marker) so they are redone cleanly
for d in run/*/; do [ -f "$d/done" ] || rm -rf "$d"; done
chown -R ubuntu:ubuntu $W
if ! sudo -u ubuntu tmux has-session -t z2ek-drv 2>/dev/null; then
  sudo -u ubuntu tmux new-session -d -s z2ek-drv "cd $W && ./driver.sh $W/ulist_run.txt > $W/logs/driver.out 2>&1"
fi
echo "z2ek driver running; $(ls -d run/*/done 2>/dev/null | wc -l) slices already done"
