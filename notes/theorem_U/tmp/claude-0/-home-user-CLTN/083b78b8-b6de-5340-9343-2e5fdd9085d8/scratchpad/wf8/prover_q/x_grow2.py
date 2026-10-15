"""EXPLORATION (floats, search only): greedy one-node growth using the LP theta-optimizer (hopt.py).
At each step: for every one-node extension of the current core, require positivity at beta_test; prefilter by
minRe at the parent's theta (padded) after a short optimization; then fully optimize the top candidates."""
import numpy as np, sys, time
from hopt import optimize, minre
def u_at(D, beta):
    try: return np.linalg.solve(D.T + np.eye(len(D)) / beta, np.ones(len(D)))
    except np.linalg.LinAlgError: return None
def ceiling_of(D, th, betas=np.arange(2.5, 60, 0.5)):
    """largest beta on grid where positivity holds and supp(theta) is light and minRe > -1/beta"""
    f = minre(D, th); supp = np.flatnonzero(th > 1e-10); best = 0
    for b in betas:
        u = u_at(D, b)
        if u is None or u.min() <= 0 or u.sum() <= 1 or any(u[supp] >= 1): break
        if f > -1 / b: best = b
        else: break
    return best, f
def step(D, th, beta_test, rng, top=40):
    m = len(D); cands = []
    for pat in range(1 << (2 * m)):
        out = [(pat >> i) & 1 for i in range(m)]; inn = [(pat >> (m + i)) & 1 for i in range(m)]
        if sum(inn) == 0 and sum(out) > 0: continue
        E = np.zeros((m + 1, m + 1)); E[:m, :m] = D; E[m, :m] = out; E[:m, m] = inn
        u = u_at(E, beta_test)
        if u is None or u.min() <= 0 or u.sum() <= 1: continue
        light = np.flatnonzero(u < 1)
        t1 = np.append(th, 0.0) * np.isin(np.arange(m + 1), light)
        cands.append((minre(E, t1), pat, E, light, t1))
    cands.sort(key=lambda c: -c[0])
    res = []
    for f0, pat, E, light, t1 in cands[:top]:
        f, t = optimize(E, light, t1 if t1.sum() > 0 else None, iters=150)
        t2 = np.zeros(m + 1); t2[light] = 1.0
        f2, tt2 = optimize(E, light, t2, iters=150)
        if f2 > f: f, t = f2, tt2
        res.append((f, E, t))
    res.sort(key=lambda c: -c[0])
    return res, len(cands)
if __name__ == '__main__':
    sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory'); import tlib
    A = np.array(tlib.from_d6(sys.argv[1]), float); m = len(A); D = 1 - np.eye(m) - A
    th = np.array([float(x) for x in sys.argv[2].split(',')]); nsteps = int(sys.argv[3]); maxm = 11
    f, th = optimize(D, list(range(m)), th, iters=400)
    print('start m=%d minRe=%.4f ceiling=%s' % (m, f, ceiling_of(D, th)), flush=True)
    rng = np.random.default_rng(0)
    for s in range(nsteps):
        t0 = time.time()
        bt = max(3.0, -1 / f)
        res, nc = step(D, th, bt, rng)
        f, D, th = res[0]; m += 1
        print('m=%d minRe=%.4f (-1/minRe=%.3f) ceiling=%s cands=%d %.0fs theta=%s' % (m, f, -1 / f, ceiling_of(D, th), nc, time.time() - t0, np.round(th, 3).tolist()), flush=True)
        print('   D=', [''.join(str(int(x)) for x in r) for r in D], flush=True)
        print('   eig=', np.round(np.sort_complex(np.linalg.eigvals(D - np.outer(th, np.ones(m)))), 3).tolist(), flush=True)
        if m >= maxm: break
