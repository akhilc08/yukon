"""Aggregate eval_rows.json into precision@k tables.   python report.py OUT_DIR [K]

For each experiment (A = held-out band, B = held-out family), K, target offset dt (rank >= median+dt) and k-fraction,
sums hits over pools and prints hits per ranker plus the ratio to S_K (S over the same K primes = same a_p cost).
"""
from __future__ import annotations

import collections
import json
import os
import sys

NAMES = ["S_K", "S1e3", "S1e4", "lin", "quad", "gbm", "mlp", "meta", "meta_rn", "lin_f", "quad_f", "gbm_f"]


def main(out_dir, Kf=None):
    rows = json.load(open(os.path.join(out_dir, "eval_rows.json")))
    agg = collections.OrderedDict()
    for r in rows:
        exp, K, fam, split = r["pool"].split("|")
        if Kf and K not in (f"K={Kf}", Kf):
            continue
        dt = r["target"] - r["generic"]
        kf = r.get("kf", round(r["k"] / r["N"], 4))
        key = (exp, K, dt, kf)
        a = agg.setdefault(key, collections.Counter())
        a["pools"] += 1
        a["k"] += r["k"]
        a["base_hits"] += r["base"] * r["k"]
        for nm in NAMES:
            if nm in r:
                a[nm] += r[nm]
    print(f"{'exp':3} {'K':6} {'dt':>2} {'k/N':>6} {'pools':>5} {'k':>6} {'rand':>7} " + " ".join(f"{n:>6}" for n in NAMES))
    for (exp, K, dt, kf), a in agg.items():
        print(f"{exp:3} {K:6} {dt:>2} {kf:>6} {a['pools']:>5} {a['k']:>6} {a['base_hits']:>7.1f} "
              + " ".join(f"{a[n]:>6}" for n in NAMES))


def per_pool(out_dir, K, dt, kf):
    rows = json.load(open(os.path.join(out_dir, "eval_rows.json")))
    for r in rows:
        exp, KK, fam, split = r["pool"].split("|")
        if KK == f"K={K}" and r["target"] - r["generic"] == dt and abs(r["k"] / r["N"] - kf) < 1e-4:
            print(f"{exp} {fam:7} {split:18} N={r['N']:6} tgt={r['target']:2} base={r['base']:.4f} k={r['k']:5} "
                  + " ".join(f"{n}={r.get(n, '-')}" for n in NAMES))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
