"""Train and evaluate learned rank filters vs S(X).  python train.py DATA_DIR OUT_DIR [K ...]

Experiments (all rankings are within one (family, band) pool, so per-family intercepts do not matter):
  A  in-family, held-out parameter band: train on band 1 of every family, test on band 2 (and 2 -> 1)
  B  held-out family: leave-one-family-out over all families (train on both bands of the rest)
Models take only x_p = a_p/sqrt(p) (a_p := 0 at bad primes) for the first K primes.
"""
from __future__ import annotations

import json
import os
import pickle
import sys

import numpy as np

import data as D

KS = [25, 50, 100, 168]
KFRAC = [0.001, 0.005, 0.01, 0.05]


def featurize(ap, bad, K, quad=False):
    x = np.where(bad[:, :K], 0, ap[:, :K]) / D.SQ[:K]
    return np.hstack([x, x * x]) if quad else x


def target(ds):
    """r centred per (family, band); r=-1 (label timeout) rows are dropped by caller."""
    y = ds.r.astype(float).copy()
    for f in np.unique(ds.fam):
        for b in (1, 2):
            m = (ds.fam == f) & (ds.band == b)
            if m.any():
                y[m] -= y[m].mean()
    return y


def fit_models(ds, K, which=("lin", "quad", "gbm", "mlp")):
    from sklearn.linear_model import Ridge
    y = target(ds)
    out = {}
    if "lin" in which:
        out["lin"] = Ridge(alpha=1.0).fit(featurize(ds.ap, ds.bad, K), y)
    if "quad" in which:
        out["quad"] = Ridge(alpha=1.0).fit(featurize(ds.ap, ds.bad, K, True), y)
    if "gbm" in which:
        import lightgbm as lgb
        out["gbm"] = lgb.LGBMRegressor(n_estimators=150, num_leaves=15, learning_rate=0.08, min_child_samples=50,
                                       subsample=0.8, subsample_freq=1, colsample_bytree=0.8, verbose=-1,
                                       n_jobs=6).fit(featurize(ds.ap, ds.bad, K), y)
    if "mlp" in which:
        from sklearn.neural_network import MLPRegressor
        n = min(len(y), 200000)
        idx = np.random.default_rng(0).permutation(len(y))[:n]
        out["mlp"] = MLPRegressor(hidden_layer_sizes=(16,), alpha=1e-3, max_iter=60, early_stopping=True,
                                  random_state=0).fit(featurize(ds.ap[idx], ds.bad[idx], K), y[idx])
    return out


def meta_x(ds, rn):
    x = [ds.S[:, 0], ds.S[:, 1], ds.S[:, 2], ds.logD, ds.nbad, ds.H]
    if rn:
        x.append(ds.rootno)
    return np.stack(x, 1)


def fit_meta(ds, rn):
    import lightgbm as lgb
    return lgb.LGBMRegressor(n_estimators=200, num_leaves=15, learning_rate=0.05, min_child_samples=50, verbose=-1,
                             n_jobs=6).fit(meta_x(ds, rn), target(ds))


def predict(name, model, ap, bad, K):
    return model.predict(featurize(ap, bad, K, name == "quad"))


def topk_hits(score, hit, k):
    idx = np.argsort(-score, kind="stable")[:k]
    return int(hit[idx].sum()), idx


def eval_pool(ds, scores: dict, rows: list, tag: str):
    """scores: name -> score vector over ds. Appends result dicts to rows."""
    g = int(np.median(ds.r))
    N = len(ds)
    for dt in (1, 2, 3, 4):
        hit = ds.r >= g + dt
        nh = int(hit.sum())
        if nh < 5:
            continue
        for kf in KFRAC:
            k = max(1, int(round(kf * N)))
            rec = dict(pool=tag, N=N, generic=g, target=g + dt, base=nh / N, k=k, kf=kf)
            for nm, s in scores.items():
                h, idx = topk_hits(s, hit, k)
                rec[nm] = h
                rec[nm + "_lab_ms"] = float(ds.tl[idx].sum())
            rows.append(rec)


def run(data_dir, out_dir, ks=KS):
    os.makedirs(out_dir, exist_ok=True)
    ds = D.load(data_dir)
    ds = ds.take(ds.r >= 0)
    print("loaded", len(ds), {f: int((ds.fam == f).sum()) for f in np.unique(ds.fam)}, flush=True)
    summ = {}
    for f in np.unique(ds.fam):
        for b in (1, 2):
            m = (ds.fam == f) & (ds.band == b)
            if m.any():
                rr = ds.r[m]
                summ[f"{f}.{b}"] = dict(n=int(m.sum()), rank_hist=np.bincount(rr).tolist(),
                                        frac_r_lt_R=float(((ds.R[m] > rr)).mean()),
                                        lab_ms_mean=float(ds.tl[m].mean()), feat_ms_mean=float(ds.tf[m].mean()),
                                        H_mean=float(ds.H[m].mean()))
    json.dump(summ, open(os.path.join(out_dir, "dataset_summary.json"), "w"), indent=1)
    rows = []
    fams = list(np.unique(ds.fam))
    for K in ks:
        # A: held-out band
        for trb, teb in ((1, 2), (2, 1)):
            tr = ds.take(ds.band == trb)
            models = fit_models(tr, K)
            for f in fams:
                te = ds.take((ds.fam == f) & (ds.band == teb))
                if len(te) < 1000:
                    continue
                sc = {"S_K": D.s_of(te.ap, te.bad, K), "S1e3": te.S[:, 2], "S1e4": te.S[:, 4]}
                for nm, mdl in models.items():
                    sc[nm] = predict(nm, mdl, te.ap, te.bad, K)
                trf = tr.take(tr.fam == f)   # same-family model (what a worker sieving one family would train)
                if len(trf) >= 2000:
                    for nm, mdl in fit_models(trf, K, ("lin", "quad", "gbm")).items():
                        sc[nm + "_f"] = predict(nm, mdl, te.ap, te.bad, K)
                eval_pool(te, sc, rows, f"A|K={K}|{f}|train_b{trb}_test_b{teb}")
            print("A done", K, trb, flush=True)
        # B: leave one family out
        for f in fams:
            tr = ds.take(ds.fam != f)
            models = fit_models(tr, K, ("lin", "quad", "gbm"))
            for b in (1, 2):
                te = ds.take((ds.fam == f) & (ds.band == b))
                if len(te) < 1000:
                    continue
                sc = {"S_K": D.s_of(te.ap, te.bad, K), "S1e3": te.S[:, 2], "S1e4": te.S[:, 4]}
                for nm, mdl in models.items():
                    sc[nm] = predict(nm, mdl, te.ap, te.bad, K)
                eval_pool(te, sc, rows, f"B|K={K}|{f}|heldout_b{b}")
        print("B done", K, flush=True)
        json.dump(rows, open(os.path.join(out_dir, "eval_rows.json"), "w"))
    # meta: learned combination of summary stats (not hot-loop: needs S(1000), disc, root number)
    for exp in ("A", "B"):
        splits = ([(ds.band == 1, (ds.band == 2), "b1->b2"), (ds.band == 2, (ds.band == 1), "b2->b1")] if exp == "A"
                  else [(ds.fam != f, ds.fam == f, f"heldout_{f}") for f in fams])
        for trm, tem, nm in splits:
            tr = ds.take(trm)
            mm = {k: fit_meta(tr, k == "meta_rn") for k in ("meta", "meta_rn")}
            for f in fams:
                for b in (1, 2):
                    te = ds.take(tem & (ds.fam == f) & (ds.band == b))
                    if len(te) < 1000:
                        continue
                    sc = {"S1e3": te.S[:, 2], "S1e4": te.S[:, 4]}
                    for k, mdl in mm.items():
                        sc[k] = mdl.predict(meta_x(te, k == "meta_rn"))
                    eval_pool(te, sc, rows, f"{exp}|meta|{f}|{nm}_b{b}")
        print("meta done", exp, flush=True)
    json.dump(rows, open(os.path.join(out_dir, "eval_rows.json"), "w"))
    # final models on everything (shipped)
    final = {}
    for K in ks:
        final[K] = fit_models(ds, K)
    pickle.dump(final, open(os.path.join(out_dir, "final_models.pkl"), "wb"))
    print("saved", flush=True)


if __name__ == "__main__":
    ks = [int(k) for k in sys.argv[3:]] or KS
    run(sys.argv[1], sys.argv[2], ks)
