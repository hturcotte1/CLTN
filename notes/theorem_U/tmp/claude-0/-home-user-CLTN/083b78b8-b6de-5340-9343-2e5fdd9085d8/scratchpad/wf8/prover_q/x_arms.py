"""EXPLORATION (floats, search only): hub + 3-cycle H-core with 'arms' (H node h_i child of a chosen H node,
Z node z_i child of h_i) and Z-arcs z_i -> S_i.  Greedy: add arms one at a time choosing (parent, S_i) maximizing the
LP-optimized ceiling with positivity at the achieved beta."""
import numpy as np, sys, itertools, time
from hopt import optimize
def ceiling(D, th0=None, restarts=3):
    m = len(D); best = (-9.0, None); rng = np.random.default_rng(1)
    starts = [np.ones(m) / m] + ([th0] if th0 is not None else []) + [rng.dirichlet(np.ones(m) * 0.3) for _ in range(restarts - 1)]
    for t0 in starts:
        f, th = optimize(D, list(range(m)), t0, iters=150)
        if f >= 0: f = -1e-6
        b = -1 / f * 0.999
        try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
        except np.linalg.LinAlgError: continue
        supp = np.flatnonzero(th > 1e-9)
        if u.min() > 0 and u.sum() > 1 and all(u[supp] < 1) and f > best[0]: best = (f, th)
    return best
def base():
    # nodes 0,1,2 = cycle 0->1->2->0 (0 = hub)
    D = np.zeros((3, 3)); D[0, 1] = D[1, 2] = D[2, 0] = 1; return D, [0, 1, 2], []
def add_arm(D, Hn, Zn, par, S):
    m = len(D); E = np.zeros((m + 2, m + 2)); E[:m, :m] = D
    h, z = m, m + 1
    E[par, h] = 1; E[h, z] = 1
    for v in S: E[z, v] = 1
    return E, Hn + [h], Zn + [z]
if __name__ == '__main__':
    K = int(sys.argv[1]); maxS = int(sys.argv[2])
    D, Hn, Zn = base(); th = None
    for i in range(K):
        t0 = time.time(); m = len(D); best = None; cnt = 0
        targets = list(range(m)) + [m]           # z may also point to its own H parent h (= m)
        for par in Hn:
            for k in range(1, maxS + 1):
                for S in itertools.combinations(targets, k):
                    E, H2, Z2 = add_arm(D, Hn, Zn, par, S); cnt += 1
                    f, t = ceiling(E, None if th is None else np.append(th, [0, 0]))
                    if best is None or f > best[0]: best = (f, E, H2, Z2, t, par, S)
        f, D, Hn, Zn, th, par, S = best
        print('arm %d: m=%d minRe=%.4f ceiling~%.3f parent=%d S=%s  (%d tried, %.0fs)' % (i + 1, len(D), f, -1 / f, par, S, cnt, time.time() - t0), flush=True)
        print('   theta', np.round(th, 3).tolist(), flush=True)
        print('   D', [''.join(str(int(x)) for x in r) for r in D], flush=True)
