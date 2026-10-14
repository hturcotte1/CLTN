"""
HEURISTIC ONLY (floating point). Lemma M restricted to the positivity set.
P := eps' I + N, eps' = 1/beta in (0,1/2); positivity <=> adj(P) 1 > 0.
stable at (beta, tau) <=> g(tau) = min Re spec(N + tau J) > -eps'.
legal <=> tau > eps'^2/(1-2 eps').
"""
import numpy as np
from scan_float import g_curve, g_inf, gP2, N_from_A

def positivity_grid(Ns, epss):
    # returns (B, len(epss)) boolean: adj(eps I + N) 1 > 0
    B, m, _ = Ns.shape
    I = np.eye(m)
    out = np.zeros((B, len(epss)), bool)
    one = np.ones(m)
    for k, e in enumerate(epss):
        P = Ns + e * I
        d = np.linalg.det(P)
        ok = np.abs(d) > 1e-12
        x = np.zeros((B, m))
        x[ok] = np.linalg.solve(P[ok], np.broadcast_to(one, (ok.sum(), m))[..., None])[..., 0]
        a = x * d[:, None]
        out[:, k] = ok & (a > 1e-12).all(axis=1)
    return out

def check_pos(Ns, taus, epss, legal=True, tol=1e-9):
    g = g_curve(Ns, taus)          # (B, T)
    gi = g_inf(Ns)                 # (B,)
    pos = positivity_grid(Ns, epss)  # (B, E)
    L = -epss                      # levels
    B = len(Ns)
    suf = np.minimum.accumulate(g[:, ::-1], axis=1)[:, ::-1]
    sufnext = np.concatenate([suf[:, 1:], np.full((B, 1), np.inf)], axis=1)
    sufnext = np.minimum(sufnext, gi[:, None])
    stable = g[:, :, None] > L[None, None, :] + tol           # (B,T,E)
    if legal:
        leg = (taus[:, None] > (epss**2 / (1 - 2 * epss))[None, :] * (1 + 1e-9))  # (T,E)
    else:
        leg = np.ones((len(taus), len(epss)), bool)
    later_unstable = sufnext[:, :, None] < L[None, None, :] - tol
    inf_unstable = gi[:, None, None] < L[None, None, :] - tol
    base = stable & leg[None] & pos[:, None, :]
    vs = (base & later_unstable).any(axis=(1, 2))
    vw = (base & inf_unstable).any(axis=(1, 2))
    anystable = base.any(axis=(1, 2))
    return vs, vw, anystable
