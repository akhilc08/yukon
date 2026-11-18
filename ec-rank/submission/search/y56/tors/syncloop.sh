#!/bin/bash
# syncloop.sh: push /work/y56/tors results to s3 every 10 minutes (spot boxes vanish). Light: one aws process.
cd /work/y56/tors
while true; do
  aws s3 sync /work/y56/tors s3://ecrank-y56-115135126148/results/tors --exclude '*sieve.txt' --exclude '*.tgz' --exclude '*/rin_0?' --exclude '_s3restore/*' --quiet || true
  sleep 600
done
