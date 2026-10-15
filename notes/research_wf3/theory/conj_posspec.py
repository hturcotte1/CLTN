"""HEURISTIC test of Conjecture P: D a nonzero 0/1 digraph matrix (zero diagonal), 0 < r < 1/2,
u = (D^T + rI)^{-1} 1 > 0  ==>  D has an eigenvalue with Re <= -1/2.   Floats; random D, m = 3..12."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]))
best = (-9, None); tested = 0
for trial in range(int(sys.argv[2])):
    m = int(rng.integers(3, 13))
    D = (rng.random((m, m)) < rng.uniform(0.05, 0.6)).astype(int); np.fill_diagonal(D, 0)
    if D.sum() == 0: continue
    ev = np.linalg.eigvals(D).real.min()
    for r in (0.05, 0.15, 0.3, 0.45, 0.499):
        try: u = np.linalg.solve(D.T + r * np.eye(m), np.ones(m))
        except np.linalg.LinAlgError: continue
        if u.min() <= 1e-12: continue
        tested += 1
        if ev > best[0]: best = (ev, (D.tolist(), r))
print('positive cases tested', tested, ' max min-Re spec(D):', best[0]); print(best[1])
