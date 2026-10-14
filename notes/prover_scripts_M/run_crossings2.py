import numpy as np, sys
from crossings import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
npos = 0; ncross = 0; stats = {}; ex = []
for trial in range(40000):
    m = int(rng.integers(3, 9))
    p = rng.uniform(0.0, 1.0)
    A = (rng.random((m, m)) < p).astype(float); np.fill_diagonal(A, 0)
    N = np.ones((m, m)) - np.eye(m) - A.T
    beta = 2 + rng.exponential(2.0)
    a = a_vec(N, beta)
    if not (a > 1e-9).all(): continue
    npos += 1
    cr = crossings(N, beta)
    for (t, y, dl) in cr:
        leg = t * (beta - 2) > 1
        key = (leg, dl > 0)
        stats[key] = stats.get(key, 0) + 1
        if len(ex) < 6 and not dl > 0: ex.append((A.astype(int).tolist(), beta, t, y, dl, leg))
print('positive cases', npos, stats)
for e in ex: print(e)
