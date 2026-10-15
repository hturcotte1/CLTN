"""EXPLORATION (floats, search only): joint annealing over a core D (non-edge matrix, m nodes) and a
sink-class mixture theta (theta >= 0, sum <= 1, supported on light nodes u_i < 1).
score = min(beta_adm, -1/minRe spec(D - theta 1^T)) where beta_adm = sup beta on grid with core positivity
(u = (D^T + I/beta)^{-1} 1 > 0, sum u > 1) and every supp(theta) node light."""
import numpy as np, sys, json, time
BG = np.concatenate([np.arange(2.05, 10, 0.05), np.arange(10, 60, 0.5)])
def beta_adm(D, supp):
    m = len(D); last = 0.0
    for b in BG:
        try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
        except np.linalg.LinAlgError: return last
        if u.min() <= 1e-12 or u.sum() <= 1 or any(u[i] >= 1 - 1e-12 for i in supp): return last
        last = b
    return last
def score(D, th):
    m = len(D)
    mr = np.linalg.eigvals(D - np.outer(th, np.ones(m))).real.min()
    cap = 1e9 if mr >= 0 else -1 / mr
    supp = np.flatnonzero(th > 1e-9)
    return min(cap, beta_adm(D, supp))
def run(m, rng, iters, D0=None, th0=None, join=True):
    if D0 is None:
        D = (rng.random((m, m)) < 0.25).astype(float); np.fill_diagonal(D, 0)
    else: D = D0.copy()
    th = th0.copy() if th0 is not None else rng.dirichlet(np.ones(m)) * (rng.random() if join else 1)
    s = score(D, th); best = (s, D.copy(), th.copy()); T = 0.3
    for it in range(iters):
        E, ph = D, th
        if rng.random() < 0.5:
            i, j = rng.integers(m, size=2)
            if i == j: continue
            E = D.copy(); E[i, j] = 1 - E[i, j]
        else:
            ph = np.clip(th + rng.normal(size=m) * 0.05 * (rng.random(m) < 0.3), 0, None)
            if rng.random() < 0.1: ph[rng.integers(m)] = 0
            sm = ph.sum()
            if sm <= 1e-9: continue
            if sm > 1 or not join: ph = ph / sm
        t = score(E, ph)
        if t >= s or rng.random() < np.exp((t - s) / T):
            D, th, s = E, ph, t
            if s > best[0]: best = (s, D.copy(), th.copy())
        T = max(0.01, T * 0.9995)
    return best
if __name__ == '__main__':
    seed = int(sys.argv[1]); ms = [int(x) for x in sys.argv[2].split(',')]; iters = int(sys.argv[3]); restarts = int(sys.argv[4])
    join = (sys.argv[5] == 'join') if len(sys.argv) > 5 else True
    rng = np.random.default_rng(seed); out = []
    for m in ms:
        t0 = time.time(); bb = None
        for r in range(restarts):
            b = run(m, rng, iters, join=join)
            if bb is None or b[0] > bb[0]: bb = b
        print(m, round(bb[0], 3), 'theta', np.round(bb[2], 3).tolist(), 'D', [''.join(str(int(x)) for x in row) for row in bb[1]], round(time.time() - t0), 's', flush=True)
