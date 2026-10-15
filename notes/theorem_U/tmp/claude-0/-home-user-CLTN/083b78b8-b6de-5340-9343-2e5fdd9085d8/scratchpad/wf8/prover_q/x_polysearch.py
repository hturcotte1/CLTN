"""EXPLORATION (floats, search only): integer polynomials p = z^d + z^(d-1) + c_2 z^(d-2) + ... + c_d (c_2 >= 1, c_j >= 1)
for the in-tree / arms-forest design.  Objective: largest beta = 1/r such that (i) every root has Re > -r and
(ii) the closed-form positivity conditions (P2), (P3), (P5) of the arms-forest core hold at r.
Keeps the core small: penalise sum of coefficients."""
import numpy as np, sys, json
def conditions(c, r):
    c = np.asarray(c, float); d = len(c) - 1
    rts = np.roots(c)
    if rts.real.min() <= -r: return False
    js = np.arange(d + 1)
    Pp = float(np.sum(c * (-1.0 / r) ** js)); Pm = float(np.sum(c * (1.0 / r) ** js))
    if Pp <= 0: return False
    c2 = c[2]
    F0 = sum(c[j] * r ** (4 - j) for j in range(4, d + 1, 2))
    F1 = sum(c[j] * r ** (3 - j) for j in range(3, d + 1, 2))
    F2 = c2 - 1
    Den = F0 - r * F1 + r * r * F2 + r ** 2 + r ** 5
    if Den <= 0: return False
    if F0 - r * F1 + r * r <= 0: return False                                    # P2
    if (1 - r + r * r - r ** 3) * F0 - r * F1 + r * r * F2 + r * r - r ** 3 + r ** 4 <= 0: return False   # P3
    lam = r ** 3 * (1 - r + r * r) / Den
    if lam * max(F0, F1, F2, 1.0) >= 1: return False                              # P5 (spine bound f_max <= F_v)
    return True
def best_beta(c):
    rts = np.roots(np.asarray(c, float)); mr = rts.real.min()
    if mr >= 0: return 0.0
    lo = -mr * 1.0001
    for f in (1.0001, 1.001, 1.01, 1.03, 1.1, 1.3):
        if conditions(c, -mr * f): return 1 / (-mr * f)
    return 0.0
if __name__ == '__main__':
    d = int(sys.argv[1]); iters = int(sys.argv[2]); cmax = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    best = (0, None)
    for rs in range(20):
        c = [1, 1] + list(rng.integers(1, cmax + 1, d - 1)); b = best_beta(c)
        T = 0.5
        for it in range(iters):
            c2 = list(c); j = int(rng.integers(2, d + 1)); c2[j] = int(min(cmax, max(1, c2[j] + rng.choice([-2, -1, 1, 2]))))
            b2 = best_beta(c2)
            if b2 >= b or rng.random() < np.exp((b2 - b) / T): c, b = c2, b2
            if b > best[0]: best = (b, list(c))
            T = max(0.01, T * 0.995)
    b, c = best
    print('d=%d cmax=%d best beta=%.4f  sum(c)=%d  c=%s' % (d, cmax, b, sum(c[2:]), json.dumps(c)), flush=True)
