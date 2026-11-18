"""Per-family 'refine S' models vs S(1e5) at matched a_p cost.   python train_refine.py DATA_DIR S5_DIR OUT_DIR

Features (all from a_p, p <= 1e5, plus the discriminant -- the same work S(1e5) already does):
  S(X) for X in 100..1e5 (10 cutoffs), quadratic fit of S against log X (level, slope, curvature),
  raw a_p for p <= 31, #bad primes <= 1000, log|disc|.     Variant "+rn" also gets the root number (needs factoring).
Splits per family: held-out band (train b1 -> test b2, b2 -> b1) and a random 50/50 split of the pooled bands.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

import data as D
from train import eval_pool, target

XS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]


def load_s5(d):
    m = {}
    for fn in glob.glob(os.path.join(d, "*.s5")):
        for line in open(fn):
            t = line.split()
            if len(t) == 16:
                m[(t[0],) + tuple(t[1:6])] = [float(x) for x in t[6:]]
    return m


def feats(ds, S5, rn=False):
    lx = np.log(np.array(XS, dtype=float))
    V = np.vstack([np.ones_like(lx), lx - lx.mean(), (lx - lx.mean()) ** 2]).T
    fit = np.linalg.lstsq(V, S5.T, rcond=None)[0].T          # (n,3): level, slope, curvature
    x = [S5, fit, ds.ap[:, :11].astype(float), ds.nbad[:, None], ds.logD[:, None]]
    if rn:
        x.append(ds.rootno[:, None])
    return np.hstack(x)


def gbm():
    import lightgbm as lgb
    return lgb.LGBMRegressor(n_estimators=300, num_leaves=15, learning_rate=0.04, min_child_samples=40,
                             subsample=0.8, subsample_freq=1, colsample_bytree=0.8, verbose=-1, n_jobs=6)


def main(data_dir, s5_dir, out_dir):
    from sklearn.linear_model import Ridge
    os.makedirs(out_dir, exist_ok=True)
    s5 = load_s5(s5_dir)
    ds = D.load(data_dir)
    keys = [(f,) + tuple(str(v) for v in a) for f, a in zip(ds.fam, ds.ainv)]
    have = np.array([k in s5 for k in keys]) & (ds.r >= 0)
    ds = ds.take(have)
    S5 = np.array([s5[(f,) + tuple(str(v) for v in a)] for f, a in zip(ds.fam, ds.ainv)])
    print("rows with S(1e5):", len(ds), {f: int((ds.fam == f).sum()) for f in np.unique(ds.fam)}, flush=True)
    rows = []
    rng = np.random.default_rng(1)
    for f in np.unique(ds.fam):
        mf = ds.fam == f
        splits = [(mf & (ds.band == 1), mf & (ds.band == 2), "b1->b2"), (mf & (ds.band == 2), mf & (ds.band == 1), "b2->b1")]
        half = rng.random(len(ds)) < 0.5
        splits += [(mf & half, mf & ~half, "rand_a"), (mf & ~half, mf & half, "rand_b")]
        for trm, tem, nm in splits:
            tr, te = ds.take(trm), ds.take(tem)
            if len(tr) < 2000 or len(te) < 2000:
                continue
            Str, Ste = S5[trm], S5[tem]
            y = target(tr)
            sc = {"S1e3": Ste[:, 3], "S1e4": Ste[:, 6], "S1e5": Ste[:, 9]}
            sc["ridge"] = Ridge(alpha=1.0).fit(feats(tr, Str), y).predict(feats(te, Ste))
            sc["gbm"] = gbm().fit(feats(tr, Str), y).predict(feats(te, Ste))
            sc["gbm_rn"] = gbm().fit(feats(tr, Str, True), y).predict(feats(te, Ste, True))
            # S(1e5) + learned correction only from small-prime a_p / disc data (keeps S as the backbone)
            res = y - np.polyval(np.polyfit(Str[:, 9], y, 1), Str[:, 9])
            corr = gbm().fit(feats(tr, Str)[:, 13:], res)
            sc["S+corr"] = np.polyval(np.polyfit(Str[:, 9], y, 1), Ste[:, 9]) + corr.predict(feats(te, Ste)[:, 13:])
            for b in (1, 2):
                mb = te.band == b
                if mb.sum() >= 1000:
                    eval_pool(te.take(mb), {k: v[mb] for k, v in sc.items()}, rows, f"R|{f}|{nm}|b{b}")
        print("done", f, flush=True)
    json.dump(rows, open(os.path.join(out_dir, "refine_rows.json"), "w"))
    names = ["S1e3", "S1e4", "S1e5", "ridge", "gbm", "gbm_rn", "S+corr"]
    print(f"{'fam':6} {'split':7} {'dt':>2} {'k':>5} {'rand':>6} " + " ".join(f"{n:>7}" for n in names))
    tot = {}
    for r in rows:
        _, f, nm, b = r["pool"].split("|")
        dt = r["target"] - r["generic"]
        if r["kf"] in (0.005, 0.01) and dt in (1, 2, 3):
            print(f"{f:6} {nm:7} {dt:>2} {r['k']:>5} {r['base']*r['k']:>6.1f} " + " ".join(f"{r[n]:>7}" for n in names))
            a = tot.setdefault((dt, r["kf"], nm.startswith("rand")), {n: 0 for n in names})
            for n in names:
                a[n] += r[n]
    print("\nTOTAL over families/pools (split kind: band-heldout vs random)")
    for (dt, kf, rnd), a in sorted(tot.items()):
        print(f"dt={dt} k/N={kf} {'random' if rnd else 'band  '} " + " ".join(f"{n}={a[n]}" for n in names)
              + f"  gbm/S1e5={a['gbm'] / max(1, a['S1e5']):.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
