#!/bin/bash
# nagrun.sh N : Nagao twist scan for both families on core 21 (nice), results to nag_*.txt and S3
cd /work/y56/z2r20
N=$1
python3 gend.py $N 1  | nice -n 10 taskset -c 21 ./nagao X552.fam 1000 0.45 > nag_X552_N$N.txt
python3 gend.py $N 40 | nice -n 10 taskset -c 21 ./nagao X808.fam 1000 0.45 > nag_X808_N$N.txt
aws s3 cp nag_X552_N$N.txt s3://ecrank-y56-115135126148/results/z2r20/ --quiet
aws s3 cp nag_X808_N$N.txt s3://ecrank-y56-115135126148/results/z2r20/ --quiet
