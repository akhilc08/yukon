#!/bin/bash
# m23 certify: for each new line of ranked.txt.pts ("roots|ainvs|points"), if the point count r and the
# height H beat the trivial-torsion baseline of cell 1/r (or r>=23 and H<17.88), run cert.gp (bigger point
# search + ellrank effort 2 + saturation), write cands/<cell>-<id>.ec and verify it with the pinned verifier
# (novelty ON). Results go to hits.txt and S3. Runs on core 11 between rank batches.
W=${W:-/work/y56/m23}; cd $W; mkdir -p cands; touch certdone.txt hits.txt
declare -A BASE=([16]=9.608854 [17]=11.073288 [18]=11.964586 [19]=13.902455 [20]=14.26261 [21]=14.799306 [22]=16.917582 [23]=17.880415 [24]=19.980041 [25]=21.769212)
[ -f ranked.txt.pts ] || exit 0
while IFS= read -r line; do
  grep -qxF "$line" certdone.txt && continue; echo "$line" >> certdone.txt
  roots=${line%%|*}; key=$(echo "$roots" | md5sum | cut -c1-8)
  # H and lb for this root set from ranked.txt
  rl=$(grep -F "$roots H=" ranked.txt | head -1); H=$(echo "$rl" | sed 's/.*H=\([0-9.]*\).*/\1/'); lb=$(echo "$rl" | sed 's/.*lb=\([0-9-]*\).*/\1/'); ub=$(echo "$rl" | sed 's/.*ub=\([0-9-]*\).*/\1/')
  go=0; b=${BASE[$lb]:-0}; [ "$lb" -ge 23 ] && b=17.880415
  awk -v h=$H -v b=$b 'BEGIN{exit !(h<b)}' && go=1
  [ "$ub" -ge 23 ] && awk -v h=$H 'BEGIN{exit !(h<17.880415)}' && go=1
  echo "$(date -u +%FT%T) $key lb=$lb ub=$ub H=$H go=$go" >> hits.txt
  [ -n "$FORCE" ] && go=1
  [ $go = 1 ] || continue
  # cheap novelty/duplicate gate first: minimal model only, no points, against a temp copy of the board
  ai=$(echo "ROOT=\"/work/ec-rank\"; read(\"/work/ec-rank/submission/search/y56/m23/lib.gp\"); a=eval(Str(\"[\",strjoin(strsplit(\"$roots\",\" \"),\",\"),\"]\")); if(#a==12||#a==10, E=mesE(a); print(E.a1,\" \",E.a2,\" \",E.a3,\" \",E.a4,\" \",E.a6))" | taskset -c ${RCORES:-11} /work/ec-rank/tools/pari/bin/gp -q -s 256M 2>/dev/null | tail -1)
  tb=$(mktemp -d); cp /work/ec-rank/submission/board/*.ec $tb/; echo "E $ai" > $tb/zz-m23-$key.ec
  dv=$(cd /work/ec-rank && taskset -c ${RCORES:-11} python3 -m ec_verify $tb 2>&1 | head -c 300 | tr '\n' ' '); rm -rf $tb
  case "$dv" in REJECTED*) echo "$(date -u +%FT%T) $key DROP $dv" >> hits.txt; aws s3 cp --quiet hits.txt s3://ecrank-y56-115135126148/results/m23/hits.txt; continue;; esac
  fn=$W/cands/cand-$key.ec
  echo "L=\"$line\"; pp=strsplit(L,\"|\"); a=eval(Str(\"[\",strjoin(strsplit(pp[1],\" \"),\",\"),\"]\")); EXTRA=eval(pp[3]); ROOT=\"/work/ec-rank\"; B=${CB:-30000}; TARGET=$lb; FN=\"$fn\"; HDR=Str(\"# Mestre ten-root construction (own search, y56/m23 m23f 2D sieve): p(x)=prod(x-a_i)=q^2-r, E=Jac(y^2=r(x)), minimal model\\n# roots a=\",a,\"\\n# points: hyperellratpoints on the quartic + ellrank + ellsaturation; search log s3://ecrank-y56-115135126148/results/m23/\\n# PARI/GP \",strjoin(apply(z->Str(z),version()[1..3]),\".\"),\", $(date -u +%F)\"); read(\"$W/cert.gp\")" | \
    taskset -c ${RCORES:-11} /work/ec-rank/tools/pari/bin/gp -q -D parisizemax=8G >> cert.log 2>&1
  if [ -f $fn ]; then
    v=$(cd /work/ec-rank && taskset -c ${RCORES:-11} python3 -m ec_verify $fn 2>&1 | python3 -c 'import json,sys; t=sys.stdin.read();
try:
  d=json.loads(t); print(d.get("cell"), d.get("naive_height"), d.get("contribution"), d.get("certified_rank"))
except Exception: print("ERR", t[:300].replace("\n"," "))')
    echo "$(date -u +%FT%T) $key VERIFY $v" >> hits.txt
    cell=$(echo $v | awk '{print $1}' | tr '/' '-'); [ -n "$cell" ] && cp $fn $W/cands/$cell-$key.ec
    contrib=$(echo $v | awk '{print $3}')
    if awk -v c="$contrib" 'BEGIN{exit !(c+0>0)}'; then   # below baseline: board duplicate check on a temp copy
      tb=$(mktemp -d); cp /work/ec-rank/submission/board/*.ec $tb/; cp $fn $tb/zz-m23-$key.ec
      bv=$(cd /work/ec-rank && taskset -c ${RCORES:-11} python3 -m ec_verify $tb 2>&1 | python3 -c 'import sys,json; t=sys.stdin.read().strip()
try: print("OK board_total=%s" % json.loads(t)["metrics"]["board_total"])
except Exception: print(t[:300].replace("\n"," "))')
      echo "$(date -u +%FT%T) $key BOARDCHECK $bv" >> hits.txt; rm -rf $tb
    fi
  fi
  aws s3 sync --quiet $W/cands s3://ecrank-y56-115135126148/results/m23/cands/
  aws s3 cp --quiet hits.txt s3://ecrank-y56-115135126148/results/m23/hits.txt
done < ranked.txt.pts
