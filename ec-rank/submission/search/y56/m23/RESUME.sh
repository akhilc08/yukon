#!/bin/bash
# One-command resume of the y56/m23 pipeline from S3 state (run as root, e.g. via SSM):
#   aws s3 cp s3://ecrank-y56-115135126148/results/m23/code/RESUME.sh /tmp/ && bash /tmp/RESUME.sh
# Restores code + survivors + ranked/certified lists, rebuilds done lists, continues seed rounds after the
# last logged ones, rebuilds m23f for this CPU, keeps the absolute DEADLINE from S3, and starts (if not
# already running, as tmux sessions m23-*) sieve.sh (cores 0-10), sieve2.sh (cores 22-27) and rankloop.sh (+certify.sh, core 11).
set -u
S3=s3://ecrank-y56-115135126148/results/m23
W=/work/y56/m23; REPO=/work/ec-rank/submission/search/y56/m23
mkdir -p $W/out $W/out2 $W/cands $W/code; cd $W
# --- code: S3 code/ wins over the repo copy
aws s3 sync --quiet $S3/code/ $W/code/
for f in m23f.c lib.gp rank.gp cert.gp certify.sh rankloop.sh sieve.sh sieve2.sh; do
  if [ -f $W/code/$f ]; then cp $W/code/$f $W/$f; else cp $REPO/$f $W/$f; fi; done
mkdir -p $REPO; cp $W/lib.gp $REPO/lib.gp     # rank.gp/certify read lib.gp from the repo path
chmod +x *.sh
# --- state
aws s3 sync --quiet $S3/out/ $W/out/; aws s3 sync --quiet $S3/out2/ $W/out2/; aws s3 sync --quiet $S3/cands/ $W/cands/
for f in ranked.txt ranked.txt.pts hits.txt log.txt log2.txt; do aws s3 cp --quiet $S3/$f $W/$f 2>/dev/null || touch $W/$f; done
{ echo dummy; cat ranked.txt; } > done.txt          # rankloop keys on the first 10 fields (= roots)
: > certdone.txt
while IFS= read -r l; do k=$(echo "${l%%|*}" | md5sum | cut -c1-8); grep -q " $k " hits.txt && echo "$l" >> certdone.txt; done < ranked.txt.pts
R1=$(grep -o 'seed=[0-9]*' log.txt | cut -d= -f2 | awk '$1<500000{r=int($1/1000); if(r>m)m=r} END{print (m?m+1:10)}')
R2=$(grep -o 'seed=[0-9]*' log2.txt | cut -d= -f2 | awk '{r=int($1/1000); if(r>m)m=r} END{print (m?m+1:500)}')
# --- build for this CPU
gcc -O3 -mcpu=native -fopenmp -Wno-empty-body -o m23f m23f.c -lm || exit 1
chown -R ubuntu:ubuntu /work/y56/m23 $REPO
# --- deadline (absolute, shared through S3)
if aws s3 cp --quiet $S3/code/DEADLINE $W/DEADLINE 2>/dev/null; then :; else echo $(( $(date +%s) + ${HOURS:-37}*3600 )) > $W/DEADLINE; aws s3 cp --quiet $W/DEADLINE $S3/code/DEADLINE; fi
DL=$(cat $W/DEADLINE); rm -f $W/STOP; [ $(date +%s) -ge $DL ] && { echo "deadline passed"; touch $W/STOP; exit 0; }
# launch detached in tmux sessions owned by ubuntu (they survive the SSM command that started them)
tl(){ sudo -u ubuntu tmux has-session -t "$1" 2>/dev/null && { echo "$1 already running"; return; }; sudo -u ubuntu tmux new-session -d -s "$1" "cd $W; $2"; }
tl m23-watch "while [ \$(date +%s) -lt $DL ]; do sleep 300; done; touch $W/STOP"
tl m23-sieve "ROUND0=$R1 NT=11 CORES=0-10 CHUNK=40000000000 T2=82 ./sieve.sh > sieve.out 2>&1"
tl m23-sieve2 "OUTD=out2 LOGF=log2.txt ROUND0=$R2 NT=6 CORES=22-27 CHUNK=20000000000 T2=82 MLIST='110:3 130:3 150:3 160:4 110:0 180:3 130:4 210:3' ./sieve2.sh > sieve2.out 2>&1"
tl m23-rank "RT=82 S2MIN=118 BATCH=30 ./rankloop.sh > rank.out 2>&1"
echo "resumed: R1=$R1 R2=$R2 deadline=$(date -u -d @$DL +%FT%TZ) done=$(wc -l < done.txt) certdone=$(wc -l < certdone.txt)"
