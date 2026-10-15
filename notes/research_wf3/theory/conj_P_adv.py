"""HEURISTIC adversarial search against Conjecture P: maximise min Re spec(D) over 0/1 zero-diagonal D (m = 4..11)
subject to: exists r in a grid of (0, 1/2) with u = (D^T + rI)^{-1} 1 > 0.  Single-arc flips, random restarts."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]))
RS = np.linspace(0.01, 0.499, 60)
def feas(D):
    m = len(D)
    for r in RS:
        try: u = np.linalg.solve(D.T + r * np.eye(m), np.ones(m))
        except np.linalg.LinAlgError: continue
        if u.min() > 1e-10: return True
    return False
def score(D):
    if D.sum() == 0: return -9
    v = np.linalg.eigvals(D).real.min()
    return v if feas(D) else v - 5
best_all = -9
for rs in range(int(sys.argv[2])):
    m = int(rng.integers(4, 12))
    D = (rng.random((m, m)) < 0.3).astype(int); np.fill_diagonal(D, 0)
    s = score(D)
    for it in range(400):
        i, j = rng.integers(0, m, 2)
        if i == j: continue
        E = D.copy(); E[i, j] ^= 1
        se = score(E)
        if se >= s: D, s = E, se
    if s > best_all:
        best_all = s; print('restart', rs, 'm', m, 'best feasible min-Re', s, flush=True)
print('FINAL max min-Re over feasible:', best_all)
