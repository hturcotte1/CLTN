"""EXPLORATION (floats): bouquet D-graphs (hub 0, p petals = directed paths of length K back to the hub)."""
import numpy as np, sys
from hopt import optimize, minre
def petal_D(p, Ks):
    n = 1 + sum(Ks); D = np.zeros((n, n)); idx = 1
    for K in Ks:
        prev = 0
        for k in range(K):
            D[prev, idx] = 1; prev = idx; idx += 1
        D[prev, 0] = 1
    return D
for Ks in [(2,2),(3,3),(4,4),(5,5),(3,3,3),(4,4,4),(6,6),(2,2,2,2),(5,5,5),(8,8)]:
    D = petal_D(len(Ks), Ks); m = len(D)
    beta = 3.0; u = np.linalg.solve(D.T + np.eye(m) / beta, np.ones(m))
    best = (-9, None)
    for s in range(6):
        th0 = np.random.default_rng(s).random(m)
        f, th = optimize(D, list(range(m)), th0 / th0.sum(), iters=300)
        if f > best[0]: best = (f, th)
    print(Ks, 'm', m, 'pos(beta=3)', u.min() > 0, 'minRe', round(best[0], 4), '-1/minRe', round(-1 / best[0], 3))
