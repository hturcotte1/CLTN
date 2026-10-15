"""
heur_mixed_cap.py -- HEURISTIC (floats; output is not a certificate).  Question: for sink-class families over an m-node core
that MIX a join class (b = 0) with non-join classes, can the limit ceiling beta_inf = -1/min Re spec(D - bbar 1^T) exceed m?
(Proved caps: no join class -> beta_inf <= m; joins only -> beta_inf <= m - 1.)
Exhaustive over all 0/1 zero-diagonal D for m <= 4 (and random D for m = 5, 6), bbar = w * b with b a 0/1 vector (one non-join
class of weight w, rest join), w on a grid; positivity at beta (u = (D^T + rI)^{-1} 1 > 0, U > 1, b^T u < 1) checked at a beta
just below the ceiling.  Reports the largest ceiling/m found.
"""
import itertools, json, sys
import numpy as np
rng = np.random.default_rng(1)
ws = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
best = {}


def ceiling(D, bbar):
    L = D - np.outer(bbar, np.ones(len(D)))
    mr = min(np.linalg.eigvals(L).real)
    if mr >= -1e-12: return np.inf
    return -1.0 / mr


def positive(D, b, beta):
    r = 1.0 / beta; m = len(D)
    try:
        u = np.linalg.solve(D.T + r * np.eye(m), np.ones(m))
    except np.linalg.LinAlgError:
        return False
    return bool(np.all(u > 1e-12) and u.sum() > 1 + 1e-12 and b @ u < 1 - 1e-12)


def scan(m, Ds):
    top = (0, None)
    for D in Ds:
        if not D.any(): continue
        for bt in itertools.product([0, 1], repeat=m):
            b = np.array(bt, float)
            for w in ws:
                c = ceiling(D, w * b)
                if not np.isfinite(c) or c <= 2: continue
                beta = 2 + 0.999 * (c - 2)
                if positive(D, b if w > 0 else 0 * b, beta):
                    if c / m > top[0]:
                        top = (c / m, {'D': D.astype(int).tolist(), 'b': list(bt), 'w': w, 'ceiling': c})
    return top


for m in (2, 3, 4):
    offs = [(i, j) for i in range(m) for j in range(m) if i != j]
    Ds = []
    for bits in itertools.product([0, 1], repeat=len(offs)):
        D = np.zeros((m, m))
        for (i, j), v in zip(offs, bits): D[i, j] = v
        Ds.append(D)
    best[m] = scan(m, Ds)
    print('m', m, 'max ceiling/m', round(best[m][0], 4), best[m][1], flush=True)
for m in (5, 6):
    Ds = []
    for _ in range(3000 if m == 5 else 1500):
        D = (rng.random((m, m)) < rng.uniform(0.15, 0.6)).astype(float); np.fill_diagonal(D, 0); Ds.append(D)
    best[m] = scan(m, Ds)
    print('m', m, 'max ceiling/m', round(best[m][0], 4), best[m][1], flush=True)
json.dump({str(k): v for k, v in best.items()}, open('heur_mixed_cap.json', 'w'), indent=1, default=str)
