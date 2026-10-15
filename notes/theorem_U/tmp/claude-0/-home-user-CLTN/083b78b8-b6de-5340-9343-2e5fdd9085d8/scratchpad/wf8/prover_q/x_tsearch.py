"""EXPLORATION (floats, search only): H/Z designs with H = 3-cycle (0->1->2->0) + K children of the hub 2 (nodes 3..K+2).
State = integer matrix T (h x h, entries 0..tmax).  theta restricted to H.
Objective: LP-optimized min Re spec(D - theta 1^T) over theta, with positivity at beta = -1/minRe (supp theta light,
u > 0, sum u > 1).  Annealing over T entries."""
import numpy as np, sys, time, json
from hopt import optimize
from hz_lib import build
def evaluate(par, T, th0=None):
    D = build(par, T).astype(float); m = len(D); h = len(par)
    light = list(range(h))
    t0 = np.zeros(m)
    if th0 is not None: t0[:h] = th0
    else: t0[:h] = 1.0 / h
    f, th = optimize(D, light, t0, iters=120)
    if f >= 0: f = -1e-6
    b = -1 / f * 0.999
    try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
    except np.linalg.LinAlgError: return -9, th[:h]
    supp = np.flatnonzero(th > 1e-9)
    viol = max(0.0, -u.min()) + max(0.0, 1 - u.sum()) + sum(max(0.0, u[i] - 1) for i in supp)
    if viol > 0: return f - 1.0 - viol, th[:h]
    return f, th[:h]
if __name__ == '__main__':
    K = int(sys.argv[1]); iters = int(sys.argv[2]); tmax = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    par = [2, 0, 1] + [2] * K; h = len(par)
    T = np.zeros((h, h), int)
    if len(sys.argv) > 5: T = np.array(json.loads(sys.argv[5]), int)
    else:
        for i in range(3, h): T[i, rng.choice(h, 2, replace=False)] = 1
    f, th = evaluate(par, T); best = (f, T.copy(), th.copy()); Temp = 0.03
    t0 = time.time()
    for it in range(iters):
        T2 = T.copy(); i = rng.integers(h); j = rng.integers(h)
        T2[i, j] = min(tmax, max(0, T2[i, j] + rng.choice([-1, 1])))
        if (T2 == T).all(): continue
        g, th2 = evaluate(par, T2, th)
        if g >= f or rng.random() < np.exp((g - f) / Temp):
            T, f, th = T2, g, th2
            if f > best[0]:
                best = (f, T.copy(), th.copy())
        Temp = max(0.003, Temp * 0.997)
    f, T, th = best
    print('K=%d h=%d minRe=%.4f ceiling~%.3f  m=%d  %.0fs' % (K, h, f, -1 / f, len(build(par, T)), time.time() - t0), flush=True)
    print('  T=', json.dumps(T.tolist()), ' theta=', np.round(th, 4).tolist(), flush=True)
