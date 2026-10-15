"""HEURISTIC: for a core D (m nodes) maximise min Re spec(D - bbar 1^T) over bbar in conv(light sets); report the
implied limit beta = -1/minRe and check self-consistency (Pos and lightness at that beta).  Floats only."""
import numpy as np, sys, json, itertools
import tlib
from mixture_limit import light_sets
def opt(D, sets, rng, iters=3000):
    S = np.array(sets); m = len(D)
    f = lambda th: np.linalg.eigvals(D - np.outer(th @ S, np.ones(m))).real.min()
    best_th = None; best = -9
    for t in range(len(S)):
        th = np.zeros(len(S)); th[t] = 1; v = f(th)
        if v > best: best, best_th = v, th
    for t in range(200):
        th = rng.dirichlet(np.ones(len(S)) * 0.3); v = f(th)
        if v > best: best, best_th = v, th
    th = best_th.copy(); step = 0.1
    for it in range(iters):
        th2 = np.abs(th + rng.normal(0, step, len(th))); th2 /= th2.sum(); v = f(th2)
        if v > best: best, th = v, th2
        if it % 500 == 499: step *= 0.6
    return best, th
if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for s in sys.argv[1:]:
        A = tlib.from_d6(s); m = len(A)
        for beta in [3.0, 4.0, 5.0, 6.0]:
            D, sets = light_sets(A, beta)
            if not sets: print(s, beta, 'not positive'); continue
            v, th = opt(D, sets, rng)
            S = np.array(sets)
            supp = [(tuple(int(t) for t in S[i]), round(th[i], 3)) for i in np.argsort(-th)[:4] if th[i] > 1e-3]
            print(s, 'beta', beta, '#light', len(sets), 'max minRe', round(v, 5), 'implied beta_inf', round(-1 / v, 4) if v < 0 else 'inf', 'consistent' if v > -1 / beta else 'NOT', supp, flush=True)
