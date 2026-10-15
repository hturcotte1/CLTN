"""HEURISTIC: maximise q_lim over closed-twin blow-ups (integer class sizes mu) of a base graph H, using the
weighted formulas (Lemma W): moving zeros = spec of the compression of M^(1/2) D M^(1/2) onto (M^(1/2) 1)^perp
(D-picture; within-class zero-sum vectors add the eigenvalue 0), positivity: ((A^T+I)M + (gamma-1)I) y = 1, y > 0."""
import numpy as np, json, sys, itertools
import tlib
BETAS = np.concatenate([np.linspace(2.0005, 2.6, 200), np.linspace(2.603, 8.0, 900)])
def q_weighted(A, mu):
    A = np.asarray(A, float); m = len(A); mu = np.asarray(mu, float)
    D = 1 - np.eye(m) - A
    s = np.sqrt(mu); X = s[:, None] * D * s[None, :]
    v = s / np.linalg.norm(s)
    # orthonormal basis of v^perp
    Q, _ = np.linalg.qr(np.column_stack([v, np.eye(m)[:, :m - 1]]))
    U = Q[:, 1:]
    zeta = np.linalg.eigvals(U.T @ X @ U)
    mr = zeta.real.min()
    if mr <= -0.5: return -3 + mr
    bc = 1 / (-mr) if mr < 0 else np.inf
    bs = BETAS[BETAS < bc]
    if len(bs) == 0: return -2.5
    g = 1 - 1 / bs
    Mm = (A.T + np.eye(m)) * mu[None, :]
    best = None
    pmax = -np.inf
    for b, gg in zip(bs, g):
        y = np.linalg.solve(Mm + (gg - 1) * np.eye(m), np.ones(m))
        pm = y.min() / np.abs(y).max()
        pmax = max(pmax, pm)
        if pm > 1e-9: best = b
    return (best - 1) if best is not None else -2 + pmax
def climb(A, mu0, maxw, iters=200, rng=None):
    mu = list(mu0); s = q_weighted(A, mu)
    for it in range(iters):
        improved = False
        cands = []
        for i in range(len(mu)):
            for d in (-1, 1, 2, -2, 5):
                nm = mu[:]; nm[i] += d
                if 1 <= nm[i] <= maxw: cands.append(nm)
        scs = [q_weighted(A, c) for c in cands]
        k = int(np.argmax(scs))
        if scs[k] > s + 1e-9: mu, s = cands[k], scs[k]; improved = True
        if not improved: break
    return s, mu
if __name__ == '__main__':
    d6 = sys.argv[1]; maxw = int(sys.argv[2])
    A = tlib.from_d6(d6); m = len(A)
    print('base', d6, 'q at mu=1:', q_weighted(A, [1] * m))
    rng = np.random.default_rng(0)
    res = []
    for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 20):
        mu0 = [1] * m if trial == 0 else list(rng.integers(1, 4, m))
        s, mu = climb(A, mu0, maxw)
        res.append((s, mu)); print('%.4f' % s, mu, 'n =', sum(mu), flush=True)
    res.sort(key=lambda x: -x[0]); print('BEST', res[0])
