"""HEURISTIC (floats, search steering only): maximize min Re spec(D - theta 1^T) over the polytope
{theta >= 0, theta_i = 0 off `light`, sum theta <= smax} by an LP trust-region method on eigenvalue gradients."""
import numpy as np
from scipy.optimize import linprog
def eig_data(D, th):
    m = len(D); L = D - np.outer(th, np.ones(m))
    w, VL, VR = None, None, None
    import scipy.linalg as sl
    w, VL, VR = sl.eig(L, left=True, right=True)
    return w, VL, VR
def grads(w, VL, VR, m):
    # d lambda_k / d theta_i = -conj(y_i) (1^T x) / (y^* x)
    G = np.zeros((len(w), m))
    for k in range(len(w)):
        x = VR[:, k]; y = VL[:, k]
        den = np.vdot(y, x)
        if abs(den) < 1e-14: den = 1e-14
        G[k] = (-(np.conj(y)) * x.sum() / den).real
    return G
def minre(D, th):
    return np.linalg.eigvals(D - np.outer(th, np.ones(len(D)))).real.min()
def optimize(D, light, th0=None, smax=1.0, iters=200, step=0.1, tol_act=0.05, rng=None):
    m = len(D); light = list(light)
    if not light: return -9.0, np.zeros(m)
    th = np.zeros(m)
    if th0 is not None: th[light] = np.asarray(th0)[light]
    if th.sum() <= 0: th[light] = smax / len(light)
    if th.sum() > smax: th *= smax / th.sum()
    f = minre(D, th)
    for it in range(iters):
        w, VL, VR = eig_data(D, th)
        re = w.real; act = np.flatnonzero(re <= re.min() + tol_act)
        G = grads(w, VL, VR, m)
        # variables: delta_i (i in light), tau ; maximize tau
        nl = len(light)
        c = np.zeros(nl + 1); c[-1] = -1
        A_ub = []; b_ub = []
        for k in act:
            row = np.zeros(nl + 1); row[:nl] = -G[k, light]; row[-1] = 1
            A_ub.append(row); b_ub.append(re[k])
        row = np.zeros(nl + 1); row[:nl] = 1; A_ub.append(row); b_ub.append(smax - th.sum())
        bounds = [(max(-th[i], -step), step) for i in light] + [(None, None)]
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds, method='highs')
        if res.status != 0: step *= 0.5; continue
        d = np.zeros(m); d[light] = res.x[:nl]
        cand = np.clip(th + d, 0, None)
        if cand.sum() > smax: cand *= smax / cand.sum()
        g = minre(D, cand)
        if g > f + 1e-12: th, f = cand, g; step = min(step * 1.5, 0.5)
        else: step *= 0.5
        if step < 1e-7: break
    return f, th
def best_theta(D, light, smax=1.0, restarts=4, rng=None):
    rng = rng or np.random.default_rng(0); m = len(D); best = (-9.0, None)
    for s in range(restarts):
        th0 = np.zeros(m)
        if s == 0: th0[list(light)] = 1.0
        else: th0[list(light)] = rng.dirichlet(np.ones(len(light)) * 0.5)
        f, th = optimize(D, light, th0 * smax / max(th0.sum(), 1e-12), smax=smax)
        if f > best[0]: best = (f, th)
    return best
