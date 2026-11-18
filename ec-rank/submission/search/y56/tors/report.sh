#!/bin/bash
# report.sh: one-line-per-job summary of /work/y56/tors job dirs: sieve stats, ranked count, rank histogram, best H per rank>=9.
cd /work/y56/tors
for d in */; do d=${d%/}; [ -f $d/sieve.log ] || continue
  s=$(tail -1 $d/sieve.log | sed 's/scanned/sc/' | cut -c1-110)
  n=$(cat $d/rin_0?.out 2>/dev/null | wc -l); t=$(cat $d/rin.gp 2>/dev/null | wc -l)
  # rank column: psieve/rankfam outputs col5, rankg col4
  if head -c 3 $d/rin.gp 2>/dev/null | grep -q '^\[[0-9],' && [ "$(head -1 $d/rin.gp | tr -cd ',' | wc -c)" = 2 ]; then rc=5; hc=4; else rc=4; hc=3; fi
  h=$(cat $d/rin_0?.out 2>/dev/null | awk -v r=$rc -v h=$hc '$r>=9{k=$r; if(!(k in b)||$h<b[k]) b[k]=$h} END{for(k in b) printf "r%s:%.3f ",k,b[k]}')
  echo "$d | $s | ranked $n/$t | $h"
done
