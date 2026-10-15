"""HEURISTIC (floats; search only).  One-node extensions of base graphs: all 2^(2m) in/out patterns.
For each extension: alpha = max Re spec(A|perp) (batched eigvals of the compressed matrix),
q_lim = sup{beta - 1 : 2 < beta < beta_c, w(beta) > 0} on a beta grid, and flags for alpha ~ -1/2 with an
eigenvalue ~ primitive cube root of unity.  Every hit must be certified exactly afterwards."""
import numpy as np, itertools, json, sys
def ortho_perp(n):
    Q, _ = np.linalg.qr(np.eye(n) - np.ones((n, n)) / n)
    # take n-1 columns spanning 1^perp
    U = np.linalg.svd(np.eye(n) - np.ones((n, n)) / n)[0][:, :n - 1]
    return U
def extensions(A):
    A = np.array(A, dtype=np.int8); m = len(A); n = m + 1
    pats = np.array(list(itertools.product([0, 1], repeat=2 * m)), dtype=np.int8)
    B = np.zeros((len(pats), n, n), dtype=np.int8)
    B[:, :m, :m] = A
    B[:, :m, m] = pats[:, :m]      # S: i -> v
    B[:, m, :m] = pats[:, m:]      # T: v -> i
    return B, pats
def alphas(B):
    n = B.shape[1]; U = ortho_perp(n)
    C = np.einsum('ia,kij,jb->kab', U, B.astype(float), U)
    ev = np.linalg.eigvals(C)
    return ev
def pos_margin(B, betas):
    n = B.shape[1]
    out = np.full((len(B), len(betas)), -np.inf)
    BT = np.transpose(B, (0, 2, 1)).astype(float)
    for j, b in enumerate(betas):
        g = 1 - 1 / b
        w = np.linalg.solve(BT + g * np.eye(n)[None], np.ones((len(B), n, 1)))[..., 0]
        out[:, j] = w.min(axis=1)
    return out
def qlim(B, ev, betas):
    a = ev.real.max(axis=1)
    bc = np.where(a < -0.5, 1 / (1 + np.minimum(a, -0.5 + 1e-15)), 2.0)
    q = np.full(len(B), -1.0)
    idx = np.where(a < -0.5 - 1e-12)[0]
    if len(idx):
        pm = pos_margin(B[idx], betas)
        for r, k in enumerate(idx):
            ok = (pm[r] > 1e-12) & (betas < bc[k])
            if ok.any(): q[k] = betas[ok].max() - 1
    return a, q
