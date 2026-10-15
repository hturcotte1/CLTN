"""EXPLORATION (floats, search only): H/Z designs.  H = c-cycle (D-arcs h_i -> h_{i+1}), Z nodes each with one
H-parent (D-arc parent -> z) and free D-arcs z -> v.  Annealing over parents and Z-arcs; objective = LP-optimized
max_theta min Re spec(D - theta 1^T), counted only if positivity holds at beta = -1/minRe with supp(theta) light."""
import numpy as np, sys, time, json
from hopt import optimize, minre
def build(c, parents, E):
    k = len(parents); m = c + k; D = np.zeros((m, m))
    for i in range(c): D[i, (i + 1) % c] = 1
    for j, p in enumerate(parents): D[p, c + j] = 1
    for (z, v) in E: D[c + z, v] = 1
    return D
def objective(D, th0=None):
    m = len(D)
    f, th = optimize(D, list(range(m)), th0 if th0 is not None else np.ones(m) / m, iters=120)
    if f >= 0: f = -1e-3
    b = -1 / f * 0.999
    try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
    except np.linalg.LinAlgError: return -9, th
    supp = np.flatnonzero(th > 1e-9)
    if u.min() <= 0 or u.sum() <= 1 or any(u[supp] >= 1): return -9, th
    return f, th
def anneal(c, k, rng, iters):
    m = c + k
    parents = list(rng.integers(0, c, k)); E = set()
    for z in range(k):
        for v in rng.choice(m, 2, replace=False):
            if v != c + z: E.add((z, int(v)))
    D = build(c, parents, E); f, th = objective(D); best = (f, parents[:], set(E), th); T = 0.02
    for it in range(iters):
        P2, E2 = parents[:], set(E)
        if rng.random() < 0.2: P2[rng.integers(k)] = int(rng.integers(c))
        else:
            z = int(rng.integers(k)); v = int(rng.integers(m))
            if v == c + z: continue
            E2 ^= {(z, v)}
        D2 = build(c, P2, E2); g, th2 = objective(D2, th)
        if g >= f or rng.random() < np.exp((g - f) / T):
            parents, E, f, th = P2, E2, g, th2
            if f > best[0]: best = (f, parents[:], set(E), th)
        T = max(0.002, T * 0.995)
    return best
if __name__ == '__main__':
    c = int(sys.argv[1]); ks = [int(x) for x in sys.argv[2].split(',')]; iters = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    for k in ks:
        t0 = time.time(); b = anneal(c, k, rng, iters)
        print('c=%d k=%d m=%d minRe=%.4f ceiling~%.3f  %.0fs' % (c, k, c + k, b[0], -1 / b[0] if b[0] < 0 else 0, time.time() - t0), flush=True)
        print('   parents', b[1], 'E', sorted(b[2]), 'theta', np.round(b[3], 3).tolist(), flush=True)
