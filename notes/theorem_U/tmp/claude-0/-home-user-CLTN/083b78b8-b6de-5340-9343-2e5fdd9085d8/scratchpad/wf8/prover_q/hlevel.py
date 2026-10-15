"""H-level reduction of H/Z designs (floats here: EXPLORATION ONLY; exact versions live in the certify scripts).

Design: H = {0..h-1} with parent map par (D-arc par[v] -> v), T = h x h nonneg integer matrix realised by Z leaves
(hz_lib.build).  With theta supported on H and sum theta = 1, the limit matrix L = D - theta 1^T satisfies
    chi_L(z) = z^(|Z| - h + 2) * psi_theta(z),   psi_theta(z) = 1_H^T adj(z^2 I - z P - T) theta     (Lemma HZ)
and positivity reads u_H = 1 + r^2 (T^T - r P^T - r^2 I)^{-1} 1,  u_z = (1 - u_par(z))/r.
"""
import numpy as np
from scipy.optimize import linprog
def PH(par):
    h = len(par); P = np.zeros((h, h))
    for v in range(h): P[par[v], v] = 1
    return P
def psi_basis(par, T):
    """coefficient matrix B (h x (2h-1)), B[v] = coefficients (high->low) of [1^T adj(z^2 I - zP - T)]_v"""
    h = len(par); P = PH(par); T = np.asarray(T, float); deg = 2 * h - 2
    xs = np.exp(2j * np.pi * (np.arange(deg + 1) + 0.25) / (deg + 1)) * 1.7
    V = np.vander(xs, deg + 1)              # high->low
    rows = []
    for x in xs:
        M = x * x * np.eye(h) - x * P - T
        rows.append(np.linalg.det(M) * np.linalg.solve(M.T, np.ones(h)))
    vals = np.array(rows)                     # (deg+1) x h
    B = np.linalg.solve(V, vals).T.real       # h x (deg+1)
    return np.round(B, 6)
def roots_of(B, th):
    c = th @ B
    # strip leading zeros
    k = 0
    while k < len(c) - 1 and abs(c[k]) < 1e-12: k += 1
    return np.roots(c[k:])
def minre(B, th):
    r = roots_of(B, th)
    return r.real.min() if len(r) else 0.0
def positivity(par, T, r):
    h = len(par); P = PH(par); T = np.asarray(T, float)
    try: uH = 1 + r * r * np.linalg.solve(T.T - r * P.T - r * r * np.eye(h), np.ones(h))
    except np.linalg.LinAlgError: return None
    uZ = []
    for t in range(h):
        for k in range(int(T[t].max(initial=0))): uZ.append((1 - uH[t]) / r)
    return uH, np.array(uZ)
def optimize_theta(B, light, th0, iters=300, step=0.1):
    h = B.shape[0]; light = list(light)
    th = np.zeros(h); th[light] = np.maximum(th0[light], 0)
    if th.sum() <= 0: th[light] = 1.0
    th /= th.sum()
    f = minre(B, th)
    for it in range(iters):
        c = th @ B; k = 0
        while abs(c[k]) < 1e-12: k += 1
        c = c[k:]; rts = np.roots(c); dc = np.polyder(c)
        re = rts.real; act = np.flatnonzero(re <= re.min() + 0.05)
        G = []
        for a in act:
            z = rts[a]; dp = np.polyval(dc, z)
            if abs(dp) < 1e-12: dp = 1e-12
            G.append([(-np.polyval(B[v][k:] if len(B[v]) > len(c) else B[v], z) / dp).real for v in light])
        nl = len(light)
        cobj = np.zeros(nl + 1); cobj[-1] = -1
        A_ub = [np.append(-np.array(g), 1.0) for g in G]; b_ub = [re[a] for a in act]
        A_eq = [np.append(np.ones(nl), 0.0)]; b_eq = [0.0]
        bounds = [(max(-th[v], -step), step) for v in light] + [(None, None)]
        res = linprog(cobj, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=b_eq, bounds=bounds, method='highs')
        if res.status != 0: step *= 0.5; continue
        cand = th.copy(); cand[light] += res.x[:nl]; cand = np.maximum(cand, 0); cand /= cand.sum()
        g = minre(B, cand)
        if g > f + 1e-12: th, f = cand, g; step = min(0.5, step * 1.5)
        else: step *= 0.5
        if step < 1e-8: break
    return f, th
def best_theta(par, T, starts=8, seed=0):
    """maximize min Re over theta on H (sum 1), with positivity checked at beta = -1/minRe (supp theta light)."""
    h = len(par); B = psi_basis(par, T); rng = np.random.default_rng(seed)
    best = (-9.0, None)
    for s in range(starts):
        th0 = np.ones(h) / h if s == 0 else rng.dirichlet(np.ones(h) * 0.5)
        f, th = optimize_theta(B, range(h), th0)
        if f >= 0: continue
        r = -f * 1.001
        pos = positivity(par, T, r)
        if pos is None: continue
        uH, uZ = pos
        supp = np.flatnonzero(th > 1e-9)
        ok = uH.min() > 0 and (len(uZ) == 0 or uZ.min() > 0) and uH.sum() + uZ.sum() > 1 and all(uH[supp] < 1)
        if not ok:
            # retry restricted to currently light nodes
            light = [v for v in range(h) if uH[v] < 1]
            if not light: continue
            f, th = optimize_theta(B, light, th)
            r = -f * 1.001; pos = positivity(par, T, r)
            if pos is None: continue
            uH, uZ = pos; supp = np.flatnonzero(th > 1e-9)
            ok = f < 0 and uH.min() > 0 and (len(uZ) == 0 or uZ.min() > 0) and uH.sum() + uZ.sum() > 1 and all(uH[supp] < 1)
            if not ok: continue
        if f > best[0]: best = (f, th)
    return best
