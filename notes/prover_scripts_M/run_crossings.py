import numpy as np, sys
from crossings import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
stats = {}
examples = {}
for trial in range(20000):
    m = rng.integers(3, 9)
    p = rng.uniform(0.3, 0.95)
    A = (rng.random((m, m)) < p).astype(float); np.fill_diagonal(A, 0)
    N = np.ones((m, m)) - np.eye(m) - A.T
    beta = 2 + rng.exponential(1.0)
    a = a_vec(N, beta)
    pos = bool((a > 1e-9).all())
    for (t, y, dl) in crossings(N, beta):
        leg = t * (beta - 2) > 1
        key = (pos, leg, dl > 0)
        stats[key] = stats.get(key, 0) + 1
        if pos and not dl > 0 and key not in examples:
            examples[key] = (A.astype(int).tolist(), beta, t, y, dl)
for k in sorted(stats): print('pos=%s legal=%s good=%s' % k, stats[k])
for k, v in examples.items(): print(k, v)
