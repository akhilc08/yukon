# second low-height pass with looser S(30000) threshold (cells 2/16, 2/17)
W=/work/y56/z2r20; cd $W; R=$W/run_lowh2; mkdir -p $R
for F in 552 808; do
  QM=20000; [ $F = 808 ] && QM=2000000
  TABLES=$W/tab_X${F}_30000.bin OMP_NUM_THREADS=10 taskset -c 12-21 ./ksieve X$F.fam 18.92 1 $QM 3000000 -1e9 -1e9 78 200 1000 30000 > $R/X$F.surv 2> $R/X$F.err
  ./rankpar.sh $F $R/X$F.surv 18.92 -1e9 lowh2$F
  cat rk_lowh2$F/shard_*.out > $R/X$F.rank; mkdir -p $R/hits; cp rk_lowh2$F/hits/* $R/hits/ 2>/dev/null; rm -rf rk_lowh2$F
done
aws s3 sync $R s3://ecrank-y56-115135126148/results/z2r20/run_lowh2 --quiet
