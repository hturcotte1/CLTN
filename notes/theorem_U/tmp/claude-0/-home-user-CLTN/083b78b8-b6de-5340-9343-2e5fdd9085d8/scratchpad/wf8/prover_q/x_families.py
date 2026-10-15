"""EXPLORATION (floats): parametric H/Z families; LP-optimized ceiling with positivity at the achieved beta."""
import numpy as np, sys, itertools
from hopt import optimize, minre
def ceiling(D, restarts=5):
    m = len(D); best = (-9, None)
    rng = np.random.default_rng(0)
    for s in range(restarts):
        th0 = np.ones(m) / m if s == 0 else rng.dirichlet(np.ones(m) * 0.3)
        f, th = optimize(D, list(range(m)), th0, iters=250)
        b = -1 / f * 0.999 if f < 0 else 1e6
        try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
        except np.linalg.LinAlgError: continue
        supp = np.flatnonzero(th > 1e-9)
        if u.min() > 0 and u.sum() > 1 and all(u[supp] < 1) and f > best[0]: best = (f, th)
    return best
def fam_ring(c, offs):
    """H cycle h_i -> h_{i+1}; Z child z_i of h_i; Z arcs z_i -> h_{i+o} for o in offs"""
    m = 2 * c; D = np.zeros((m, m))
    for i in range(c):
        D[i, (i + 1) % c] = 1; D[i, c + i] = 1
        for o in offs: D[c + i, (i + o) % c] = 1
    return D
if __name__ == '__main__':
    for c in range(3, 9):
        for offs in [(1,), (2,), (0,), (0, 1), (1, 2), (0, 2)]:
            D = fam_ring(c, offs); f, th = ceiling(D)
            print('ring c=%d offs=%s m=%d minRe=%.4f ceil=%.3f' % (c, offs, len(D), f, -1 / f if f < 0 else 0), flush=True)
