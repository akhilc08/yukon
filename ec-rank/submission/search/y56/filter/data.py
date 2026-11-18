"""Load gen_feats.gp output into numpy arrays (y56 learned-filter experiment)."""
from __future__ import annotations

import glob
import math
import os
from dataclasses import dataclass

import numpy as np

PRIMES = [p for p in range(2, 1000) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
assert len(PRIMES) == 168
SQ = np.sqrt(np.array(PRIMES, dtype=np.float64))
LW = np.array([math.log(p) / p for p in PRIMES])


def disc(a1: int, a2: int, a3: int, a4: int, a6: int) -> int:
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


@dataclass
class DS:
    fam: np.ndarray      # str
    band: np.ndarray     # int
    ainv: list           # list of 5-tuples of int
    H: np.ndarray
    logD: np.ndarray
    nbad: np.ndarray
    tors: np.ndarray
    rootno: np.ndarray
    S: np.ndarray        # (n,5) S(100,316,1000,3162,10000)
    r: np.ndarray        # rank lower bound label
    R: np.ndarray        # upper bound (-1 unknown)
    tf: np.ndarray       # feature ms (gp)
    tl: np.ndarray       # label ms
    ap: np.ndarray       # (n,168) int
    bad: np.ndarray      # (n,168) bool  p | disc

    def take(self, m):
        return DS(**{k: ([v[i] for i in np.flatnonzero(m)] if k == "ainv" else v[m])
                     for k, v in self.__dict__.items()})

    def __len__(self):
        return len(self.r)


def load(paths) -> DS:
    if isinstance(paths, str):
        paths = sorted(glob.glob(os.path.join(paths, "*.txt")))
    rows, seen = [], set()
    for fn in paths:
        with open(fn) as f:
            for line in f:
                t = line.split()
                if len(t) != 23:
                    continue
                key = (t[0], t[3], t[4], t[5], t[6], t[7])
                if key in seen:
                    continue
                seen.add(key)
                rows.append(t)
    n = len(rows)
    ap = np.array([[int(x) for x in t[22].split(",")] for t in rows], dtype=np.int16)
    ainv = [tuple(int(x) for x in t[3:8]) for t in rows]
    bad = np.zeros((n, 168), dtype=bool)
    for i, a in enumerate(ainv):
        D = disc(*a)
        bad[i] = [D % p == 0 for p in PRIMES]
    col = lambda j, ty=float: np.array([ty(t[j]) for t in rows])
    return DS(fam=np.array([t[0] for t in rows]), band=col(1, int), ainv=ainv, H=col(8), logD=col(9),
              nbad=col(10, int), tors=col(11, int), rootno=col(12, int),
              S=np.array([[float(x) for x in t[13:18]] for t in rows]), r=col(18, int), R=col(19, int),
              tf=col(20), tl=col(21), ap=ap, bad=bad)


def s_of(ap: np.ndarray, bad: np.ndarray, K: int) -> np.ndarray:
    """Baseline S over the first K primes, msum convention: p>=3, good primes, -a_p log p / p."""
    w = LW[:K].copy()
    w[0] = 0.0
    return -((ap[:, :K] * ~bad[:, :K]) * w).sum(1)
