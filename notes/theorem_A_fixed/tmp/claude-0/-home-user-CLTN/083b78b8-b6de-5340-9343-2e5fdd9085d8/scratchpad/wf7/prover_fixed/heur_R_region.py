"""HEURISTIC (floats steer only): for R = &F^U|m^QTq?, S = all but 2: Psi_k = z^3 g, g = z^5-3z^2-2z-1 - a Q0 - a k Q1.
For a grid of (beta, t) report the set of k in 0..K where min Re roots(g) > -1/beta."""
import numpy as np
Q0 = np.array([7,12,21,13,5.]); Q1 = np.array([1,1,3,1,1.]); B0 = np.array([1,0,0,-3,-2,-1.])
def minre(a, k):
    g = B0.copy(); g[1:] -= a*(Q0 + k*Q1)
    return min(np.roots(g).real)
ks = list(range(0, 60)) + [100, 300, 1000, 10**4, 10**6]
for beta in [2.05, 2.1, 2.15, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]:
    row = []
    for t in [1.01/(beta-2), 1.5/(beta-2), 3/(beta-2), 20, 100, 1e4]:
        if t*(beta-2) <= 1: row.append('  -  '); continue
        a = t/beta; r = 1/beta
        bad = [k for k in ks if minre(a, k) <= -r]
        row.append('t=%.3g bad:%s' % (t, (bad[:3] + ['..', bad[-1]]) if len(bad) > 3 else bad))
    print(beta, ' | '.join(row))
