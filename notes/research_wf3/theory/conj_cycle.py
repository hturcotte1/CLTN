"""HEURISTIC: is min Re spec(D) <= -1/2 for EVERY non-nilpotent 0/1 matrix with zero diagonal (no positivity)?"""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]))
best = (-9, None); cnt = 0
for trial in range(int(sys.argv[2])):
    m = int(rng.integers(2, 11))
    D = (rng.random((m, m)) < rng.uniform(0.02, 0.5)).astype(int); np.fill_diagonal(D, 0)
    ev = np.linalg.eigvals(D)
    if np.abs(ev).max() < 0.5: continue        # (numerically) nilpotent
    cnt += 1
    v = ev.real.min()
    if v > best[0]: best = (v, D.tolist())
print(cnt, best[0]); print(best[1])
