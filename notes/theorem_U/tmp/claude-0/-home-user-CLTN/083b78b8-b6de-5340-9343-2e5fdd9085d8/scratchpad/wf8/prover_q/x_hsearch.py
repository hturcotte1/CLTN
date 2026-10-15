"""EXPLORATION (floats, search only): annealing at the H level over (parent map, integer T), theta via hlevel.best_theta.
Usage: python3 x_hsearch.py h iters tmax seed [T0json] [parjson]"""
import numpy as np, sys, time, json
from hlevel import best_theta
def score(par, T):
    f, th = best_theta(par, T, starts=4)
    return f, th
if __name__ == '__main__':
    h = int(sys.argv[1]); iters = int(sys.argv[2]); tmax = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    par = [2, 0, 1] + [2] * (h - 3)
    if len(sys.argv) > 6: par = json.loads(sys.argv[6])
    if len(sys.argv) > 5 and sys.argv[5] != '-': T = np.array(json.loads(sys.argv[5]), int)
    else:
        T = np.zeros((h, h), int)
        for i in range(3, h): T[i, rng.choice(h, 2, replace=False)] = 1
    f, th = score(par, T); best = (f, T.copy(), list(par), th); Temp = 0.02; t0 = time.time()
    for it in range(iters):
        T2 = T.copy(); par2 = list(par)
        if rng.random() < 0.1 and h > 3:
            v = int(rng.integers(3, h)); par2[v] = int(rng.choice([u for u in range(h) if u != v]))
        else:
            i = rng.integers(h); j = rng.integers(h)
            T2[i, j] = min(tmax, max(0, T2[i, j] + rng.choice([-1, 1])))
        g, th2 = score(par2, T2)
        if g >= f or rng.random() < np.exp((g - f) / Temp):
            T, par, f, th = T2, par2, g, th2
            if f > best[0]: best = (f, T.copy(), list(par), th)
        Temp = max(0.002, Temp * 0.995)
    f, T, par, th = best
    print('h=%d minRe=%.4f ceiling~%.3f  %.0fs' % (h, f, -1 / f if f < 0 else 0, time.time() - t0), flush=True)
    print('  par=', json.dumps(par), ' T=', json.dumps(T.tolist()), ' theta=', None if th is None else np.round(th, 4).tolist(), flush=True)
