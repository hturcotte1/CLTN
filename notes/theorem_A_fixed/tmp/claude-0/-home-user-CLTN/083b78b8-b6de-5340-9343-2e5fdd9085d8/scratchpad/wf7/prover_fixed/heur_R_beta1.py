"""HEURISTIC: largest beta1 such that R_k is stable for all k >= 0 at every legal t (a > a1(beta1)), float scan."""
import numpy as np
Q0 = np.array([7,12,21,13,5.]); Q1 = np.array([1,1,3,1,1.]); B0 = np.array([1,0,0,-3,-2,-1.])
def minre(a, k):
    g = B0.copy(); g[1:] -= a*(Q0 + k*Q1); return min(np.roots(g).real)
ks = list(range(0, 40)) + [60, 100, 300, 1000, 10**4, 10**6]
for beta in np.arange(2.10, 2.22, 0.005):
    r = 1/beta; a1 = r*r/(1-2*r)
    worst = min(minre(a, k) + r for a in a1*np.array([1.0001,1.01,1.1,1.5,2,4,10,100,1e4]) for k in ks)
    print('%.3f %.5f' % (beta, worst))
