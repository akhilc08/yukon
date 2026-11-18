"""Sanity check on the frozen known-curve ledger (selection-biased: these are record/table curves).
python known_check.py known_ap.txt models.json
known_ap.txt lines: torsion_key rank_claimed certified_rank naive_height ap(168, comma; 0 at bad primes)
Reports Spearman(score, rank) per torsion group with >= 40 curves, for S_K and each model at K.
"""
from __future__ import annotations

import collections
import sys

import numpy as np

from scorer import Scorer


def spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b + 1e-9 * np.arange(len(b))))
    return float(np.corrcoef(ra, rb)[0, 1])


def main(fn, mj):
    rows = [l.split() for l in open(fn)]
    tk = np.array([r[0] for r in rows])
    rank = np.array([int(r[2]) if int(r[2]) >= 0 else int(r[1]) for r in rows])
    ap = np.array([[int(x) for x in r[4].split(",")] for r in rows])
    sc = Scorer(mj)
    names = sc.models()
    Ks = sorted({sc.m[n]["K"] for n in names})
    out = collections.OrderedDict()
    for t in sorted(set(tk), key=lambda s: -int((tk == s).sum())):
        m = tk == t
        if m.sum() < 40:
            continue
        rec = {"n": int(m.sum()), "ranks": f"{rank[m].min()}-{rank[m].max()}"}
        for K in Ks:
            rec[f"S{K}"] = spearman(sc.S(ap[m], K), rank[m])
        for n in names:
            rec[n] = spearman(sc.score_ap(ap[m], n), rank[m])
        out[t] = rec
    keys = list(next(iter(out.values())).keys())
    print("tors " + " ".join(f"{k:>8}" for k in keys))
    for t, rec in out.items():
        print(f"{t:4} " + " ".join(f"{v:>8.3f}" if isinstance(v, float) else f"{v:>8}" for v in rec.values()))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
