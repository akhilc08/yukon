#!/bin/bash
# m23 rank worker: periodically rank the best unranked sieve survivors on core 11.
W=/work/y56/m23; cd $W; touch done.txt ranked.txt; [ -s done.txt ] || echo "dummy" > done.txt
while [ ! -f $W/STOP ]; do
  cat out/s*.txt out2/s*.txt r12/s12_*.txt 2>/dev/null | awk -v t=${RT:-82} '(NF==14 && $14>=t){print $14" "$0} (NF==11 && $11>=t){print $11" "$0}' | sort -gr | cut -d" " -f2- | awk '!s[$0]++' > all.txt
  awk 'NR==FNR{d[$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10]=1; next} !(($1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10) in d)' done.txt all.txt | head -${BATCH:-20} > in.txt
  if [ -s in.txt ]; then
    cat in.txt >> done.txt
    echo "ROOT=\"/work/ec-rank\";IN=\"$W/in.txt\";OUT=\"$W/ranked.txt\";B=${B:-3000};EFF=1;S2MIN=${S2MIN:-118};read(\"$W/rank.gp\")" | \
      taskset -c ${RCORES:-11} /work/ec-rank/tools/pari/bin/gp -q -D parisizemax=8G 2>> rank.err
    aws s3 cp --quiet ranked.txt s3://ecrank-y56-115135126148/results/m23/ranked.txt
    [ -f ranked.txt.pts ] && aws s3 cp --quiet ranked.txt.pts s3://ecrank-y56-115135126148/results/m23/ranked.txt.pts
    ./certify.sh
  else sleep 120; fi
done
