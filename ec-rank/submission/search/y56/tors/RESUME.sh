#!/bin/bash
# RESUME.sh: one-command restart of the y56/tors pipeline on a fresh/rebooted box (run as root or ubuntu).
# 1. code from the repo copy, 2. results/state (jobs.done, hits.seen, rin.gp, .prog, .out) from s3,
# 3. compile for this CPU, 4. start syncloop.sh and driver.sh on cores 28-31 (<=4 threads) unless already running.
set -e
W=/work/y56/tors; S3=s3://ecrank-y56-115135126148
mkdir -p $W && cd $W
SRC=/work/ec-rank/submission/search/y56/tors
# state first (results), then code on top (repo copy, then the newest code tarball), so old code in results never wins
aws s3 sync $S3/results/tors $W --quiet --exclude '*sieve.txt' || true
[ -d $SRC ] && cp -r $SRC/. $W/ 2>/dev/null || true
aws s3 cp $S3/code/tors/latest.tgz /tmp/tors_latest.tgz --quiet 2>/dev/null && tar xzf /tmp/tors_latest.tgz -C $W 2>/dev/null || true
[ -d fam ] || cp -r /work/ec-rank/submission/search/y56/families fam
gcc -O3 -mcpu=native -fopenmp -o gsieve gsieve.c -lm
gcc -O3 -mcpu=native -fopenmp -o psieve psieve.c -lm
chmod +x *.sh *.py; chown -R ubuntu:ubuntu /work/y56 2>/dev/null || true
U="sudo -u ubuntu"; [ "$(id -un)" = ubuntu ] && U=""
pgrep -f "^/bin/bash ./syncloop.sh" >/dev/null || $U setsid nohup ./syncloop.sh > syncloop.log 2>&1 < /dev/null &
if pgrep -f "^/bin/bash ./driver.sh" | while read p; do [ "$(readlink /proc/$p/cwd)" = $W ] && echo y; done | grep -q y; then echo "driver already running"; else
  $U setsid nohup ./driver.sh > driver.nohup 2>&1 < /dev/null & echo "driver started"; fi
