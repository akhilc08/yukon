"""y56 learned rank filter -- python scorer (same models as y56_filter.h).

    from scorer import Scorer
    sc = Scorer()                         # loads models.json next to this file
    s = sc.score_ap(ap, "lin50")          # ap: (n, >=K) int array, a_p of 2,3,5,..., 0 at bad primes
    s = sc.score_curves(ainvs, "lin50")   # ainvs: list of [a1,a2,a3,a4,a6]; a_p via gp (ellap)
    sc.S(ap, 50)                          # baseline S over the same primes (msum convention)

Larger score = likelier high rank.  Compare scores only inside one family / parameter pool.
CLI:  python scorer.py MODEL < curves.txt   (one "a1 a2 a3 a4 a6" per line) -> "score a1 a2 a3 a4 a6"
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [p for p in range(2, 1000) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
GP = os.environ.get("Y56_GP", os.path.join(HERE, "../../../../tools/pari/bin/gp"))


class Scorer:
    def __init__(self, path: str | None = None):
        self.m = json.load(open(path or os.path.join(HERE, "model", "models.json")))
        self._lgb = {}

    def models(self):
        return sorted(self.m)

    def score_ap(self, ap, name: str) -> np.ndarray:
        j = self.m[name]
        K = j["K"]
        a = np.asarray(ap, dtype=np.float64)[:, :K]
        if j["kind"] == "linear_ap":
            return j["bias"] + a @ np.array(j["w_ap"])
        if j["kind"] == "quad_ap":
            return j["bias"] + a @ np.array(j["w1_ap"]) + (a * a) @ np.array(j["w2_ap"])
        if j["kind"] == "mlp_relu":
            h = np.maximum(a @ np.array(j["W1_ap"]) + np.array(j["b1"]), 0)
            return h @ np.array(j["W2"]) + j["b2"]
        if j["kind"] == "lgbm":
            import lightgbm as lgb
            if name not in self._lgb:
                self._lgb[name] = lgb.Booster(model_str=j["model_str"])
            return self._lgb[name].predict(a / np.sqrt(np.array(PRIMES[:K])))
        raise ValueError(j["kind"])

    @staticmethod
    def S(ap, K: int) -> np.ndarray:
        w = np.array([0.0] + [math.log(p) / p for p in PRIMES[1:K]])
        return -(np.asarray(ap, dtype=np.float64)[:, :K] @ w)

    def score_curves(self, ainvs, name: str) -> np.ndarray:
        return self.score_ap(ap_from_ainvs(ainvs, self.m[name]["K"]), name)


def ap_from_ainvs(ainvs, K: int, gp: str = GP) -> np.ndarray:
    """a_p for the first K primes via gp's ellap on the given model (0 at primes dividing the discriminant).
    Pass minimal models: bad-prime detection uses the discriminant of the model you give."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for a in ainvs:
            f.write(" ".join(str(int(x)) for x in a) + "\n")
        fn = f.name
    prog = (f'P=primes({K});L=readstr("{fn}");for(i=1,#L,E=ellinit(eval(concat(["[",strjoin(strsplit(L[i]," "),","),"]"])));'
            f'D=E.disc;print(strjoin(vector({K},j,Str(if(D%P[j],ellap(E,P[j]),0)))," ")))')
    out = subprocess.run([gp, "-q", "-s", "256M"], input=prog + "\nquit\n", capture_output=True, text=True, check=True).stdout
    os.unlink(fn)
    return np.array([[int(x) for x in l.split()] for l in out.splitlines() if l.strip()], dtype=np.int64)


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "lin50"
    cur = [l.split() for l in sys.stdin if l.strip()]
    for s, c in zip(Scorer().score_curves(cur, name), cur):
        print(f"{s:.5f} " + " ".join(c))
