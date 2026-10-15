"""EXPLORATION (floats, search only): greedy one-node growth, new node's D-out and D-in sets of size <= K.
Positivity required only at beta_test (the limit statement is at a single beta)."""
import numpy as np, sys, time, itertools, json
from hopt import optimize, minre
def u_at(D, beta):
    try: return np.linalg.solve(D.T + np.eye(len(D)) / beta, np.ones(len(D)))
    except np.linalg.LinAlgError: return None
def subsets(m, K):
    for k in range(K + 1):
        for c in itertools.combinations(range(m), k): yield c
def step(D, th, beta_test, K=4, top=30):
    m = len(D); cands = []
    outs = list(subsets(m, K)); ins = [c for c in subsets(m, K) if len(c) >= 1]
    for o in outs:
        for i in ins:
            E = np.zeros((m + 1, m + 1)); E[:m, :m] = D; E[m, list(o)] = 1; E[list(i), m] = 1
            u = u_at(E, beta_test)
            if u is None or u.min() <= 0 or u.sum() <= 1: continue
            light = np.flatnonzero(u < 1)
            t1 = np.append(th, 0.0) * np.isin(np.arange(m + 1), light)
            cands.append((minre(E, t1), E, light, t1))
    cands.sort(key=lambda c: -c[0])
    res = []
    for f0, E, light, t1 in cands[:top]:
        f, t = optimize(E, light, t1 if t1.sum() > 0 else None, iters=150)
        t2 = np.zeros(m + 1); t2[light] = 1.0
        f2, tt2 = optimize(E, light, t2, iters=150)
        if f2 > f: f, t = f2, tt2
        # positivity at -1/f with this support
        b = -1 / f * 0.999
        u = u_at(E, b)
        ok = u is not None and u.min() > 0 and u.sum() > 1 and all(u[np.flatnonzero(t > 1e-10)] < 1)
        if ok: res.append((f, E, t))
    res.sort(key=lambda c: -c[0])
    return res, len(cands)
if __name__ == '__main__':
    D = np.array([[int(c) for c in r] for r in json.loads(sys.argv[1])], float); m = len(D)
    th = np.array([float(x) for x in sys.argv[2].split(',')]); maxm = int(sys.argv[3]); K = int(sys.argv[4])
    f, th = optimize(D, list(range(m)), th, iters=400)
    print('start m=%d minRe=%.4f' % (m, f), flush=True)
    while m < maxm:
        t0 = time.time()
        res, nc = step(D, th, max(3.0, -1 / f), K)
        if not res: print('no extension'); break
        f, D, th = res[0]; m += 1
        print('m=%d minRe=%.4f (-1/minRe=%.3f) cands=%d %.0fs theta=%s' % (m, f, -1 / f, nc, time.time() - t0, np.round(th, 3).tolist()), flush=True)
        print('   D=', json.dumps([''.join(str(int(x)) for x in r) for r in D]), flush=True)
        print('   eig=', np.round(np.sort_complex(np.linalg.eigvals(D - np.outer(th, np.ones(m)))), 3).tolist(), flush=True)
