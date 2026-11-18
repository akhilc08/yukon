# stats.py RANKFILE [SIEVEFILE] : per-rank counts, best H, mean S5 (and sieve S3 if given); ranked fibres only
import sys, collections
S3 = {}
if len(sys.argv) > 2:
    for l in open(sys.argv[2]):
        w = l.split()
        if len(w) >= 6: S3[(w[0], w[1])] = float(w[5])
rows = []
for l in open(sys.argv[1]):
    w = l.split()
    if 'R' not in w: continue
    i = w.index('R'); rows.append((int(w[i+1]), int(w[i+2]), float(w[3]), float(w[5]), S3.get((w[0], w[1]), float('nan')), w[0] + '/' + w[1]))
by = collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
for k in sorted(by):
    L = by[k]; b = min(L, key=lambda x: x[2])
    print(f"rank {k}: n={len(L)} bestH={b[2]:.3f} (t={b[5]}) meanS5={sum(x[3] for x in L)/len(L):.1f} minS5={min(x[3] for x in L):.1f} minS3={min(x[4] for x in L):.1f} inexact={sum(1 for x in L if x[0]!=x[1])}")
