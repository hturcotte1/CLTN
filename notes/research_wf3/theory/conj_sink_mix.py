"""HEURISTIC: do fractional mixtures bbar of light sets (b^T u < 1) make spec(D - bbar 1^T) lie in Re > -1/2?
(the k->infinity limit of many sink classes with different in-sets).  Random cores m = 4..8, floats only."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]))
best = (-9, None); tested = 0
for trial in range(int(sys.argv[2])):
    m = rng.integers(4, 9)
    A = (rng.random((m, m)) < rng.uniform(0.5, 0.85)).astype(int); np.fill_diagonal(A, 0)
    D = 1 - np.eye(m, dtype=int) - A
    for beta in (2.05, 2.3, 2.7):
        r = 1 / beta
        try: u = np.linalg.solve(D.T + r * np.eye(m), np.ones(m))
        except np.linalg.LinAlgError: continue
        if u.min() <= 0 or u.sum() <= 1: continue
        light = [np.array([(bits >> i) & 1 for i in range(m)]) for bits in range(1, 1 << m)]
        light = [b for b in light if b @ u < 1 - 1e-9]
        if not light: continue
        L = np.array(light, float)
        for rep in range(30):
            wts = rng.dirichlet(np.ones(len(L)) * 0.3)
            bb = wts @ L
            tested += 1
            mr = np.linalg.eigvals(D - np.outer(bb, np.ones(m))).real.min()
            if mr > best[0]: best = (mr, (A.tolist(), bb.round(3).tolist(), beta))
print('tested', tested, 'max over tests of min Re spec(D - bbar 1^T):', best[0])
print(best[1])
