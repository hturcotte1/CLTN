"""HEURISTIC helpers (floats; search steering only, never a decision)."""
import numpy as np
def D_of_A(A):
    A = np.asarray(A, float); m = len(A); return 1 - np.eye(m) - A
def u_vec(D, beta):
    m = len(D)
    return np.linalg.solve(D.T + np.eye(m) / beta, np.ones(m))
def minre(D, theta):
    m = len(D)
    return np.linalg.eigvals(D - np.outer(theta, np.ones(m))).real.min()
def opt_theta(D, light, allow_join=True, starts=6, rng=None, iters=300):
    """maximize min Re spec(D - theta 1^T) over theta >= 0 supported on light, sum <= 1 (join allowed) or == 1."""
    rng = rng or np.random.default_rng(0)
    m = len(D); L = list(light); k = len(L)
    if k == 0: return -9, None
    best = (-9, None)
    def make(x):
        e = np.exp(x[:k] - x[:k].max()); p = e / e.sum()
        lam = 1 / (1 + np.exp(-x[k])) if allow_join else 1.0
        th = np.zeros(m); th[L] = lam * p; return th
    for s in range(starts):
        x = rng.normal(size=k + 1) * 2
        th = make(x); f = minre(D, th)
        step = 1.0
        for it in range(iters):
            y = x + rng.normal(size=k + 1) * step
            g = minre(D, make(y))
            if g > f: x, f = y, g
            else: step *= 0.97
            if step < 1e-3: break
        if f > best[0]: best = (f, make(x))
    return best
def ceiling(D, betas=None, allow_join=True, rng=None):
    """largest beta on grid with positivity and some admissible theta giving minre > -1/beta."""
    betas = betas if betas is not None else np.arange(2.25, 30, 0.25)
    sc, info = 0.0, None
    for beta in betas:
        try: u = u_vec(D, beta)
        except np.linalg.LinAlgError: break
        if u.min() <= 1e-12 or u.sum() <= 1: break
        light = [i for i in range(len(D)) if u[i] < 1 - 1e-12]
        f, th = opt_theta(D, light, allow_join, rng=rng)
        if f > -1 / beta: sc, info = beta, (f, th)
        else: break
    return sc, info
