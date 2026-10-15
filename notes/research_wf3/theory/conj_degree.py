"""HEURISTIC: does min in-degree >= 1 and min out-degree >= 1 (0/1, zero diagonal) alone force an eigenvalue with Re <= -1/2?
Adversarial single-flip ascent on min Re spec(D) under the degree constraint."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]))
def ok(D): return D.sum(0).min() >= 1 and D.sum(1).min() >= 1
best = -9
for rs in range(int(sys.argv[2])):
    m = int(rng.integers(3, 12))
    D = (rng.random((m, m)) < 0.35).astype(int); np.fill_diagonal(D, 0)
    s = np.linalg.eigvals(D).real.min() if ok(D) else -9
    for it in range(400):
        i, j = rng.integers(0, m, 2)
        if i == j: continue
        E = D.copy(); E[i, j] ^= 1
        se = np.linalg.eigvals(E).real.min() if ok(E) else -9
        if se >= s: D, s = E, se
    if s > best:
        best = s; print(rs, m, s, D.tolist() if s > -0.49 else '', flush=True)
print('FINAL', best)
