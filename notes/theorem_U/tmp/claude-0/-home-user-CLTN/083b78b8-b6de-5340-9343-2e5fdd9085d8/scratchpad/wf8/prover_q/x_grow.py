"""EXPLORATION (floats, search only): greedy one-node growth of a core for the sink-class limit.
State: core D (non-edge matrix) + theta (sum <= 1, supported on light nodes).  For each one-node extension
(all 2^(2m) in/out non-edge patterns), check positivity at beta* and locally re-optimize theta from the parent's;
keep the extensions with the best margin minRe(D - theta 1^T) + 1/beta*."""
import numpy as np, sys, itertools, time
def u_at(D, beta):
    m = len(D)
    try: return np.linalg.solve(D.T + np.eye(m) / beta, np.ones(m))
    except np.linalg.LinAlgError: return None
def minre(D, th):
    return np.linalg.eigvals(D - np.outer(th, np.ones(len(D)))).real.min()
def local_theta(D, th, light, rng, iters=60, step=0.05, join=True):
    m = len(D); th = th * np.isin(np.arange(m), light)
    if th.sum() <= 0: th = np.zeros(m); th[light] = 1.0 / len(light)
    f = minre(D, th)
    for it in range(iters):
        ph = th.copy(); k = rng.integers(len(light)); ph[light[k]] = max(0, ph[light[k]] + rng.normal() * step)
        if rng.random() < 0.3:
            k2 = rng.integers(len(light)); ph[light[k2]] = max(0, ph[light[k2]] + rng.normal() * step)
        s = ph.sum()
        if s <= 0: continue
        if s > 1 or not join: ph /= s
        g = minre(D, ph)
        if g > f: th, f = ph, g
        else: step *= 0.98
    return f, th
def extend(D, th, beta, rng, keep=5, join=True):
    m = len(D); res = []
    for pat in range(1 << (2 * m)):
        out = [(pat >> i) & 1 for i in range(m)]            # new node -/-> i
        inn = [(pat >> (m + i)) & 1 for i in range(m)]      # i -/-> new node
        if sum(inn) == 0 and sum(out) > 0: continue          # D-source with out-edges kills positivity
        E = np.zeros((m + 1, m + 1)); E[:m, :m] = D; E[m, :m] = out; E[:m, m] = inn
        u = u_at(E, beta)
        if u is None or u.min() <= 1e-12 or u.sum() <= 1: continue
        light = np.flatnonzero(u < 1 - 1e-12)
        if len(light) == 0: continue
        th1 = np.append(th, 0.0)
        f0 = minre(E, th1 * np.isin(np.arange(m + 1), light))
        res.append((f0, pat, E, th1, light))
    res.sort(key=lambda x: -x[0])
    best = []
    for f0, pat, E, th1, light in res[:200]:
        f, t = local_theta(E, th1, light, rng, join=join)
        best.append((f, E, t))
    best.sort(key=lambda x: -x[0])
    return best[:keep], len(res)
if __name__ == '__main__':
    sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
    import tlib
    A = np.array(tlib.from_d6(sys.argv[1]), float); m = len(A); D = 1 - np.eye(m) - A
    th = np.array([float(x) for x in sys.argv[2].split(',')]); steps = int(sys.argv[3])
    rng = np.random.default_rng(0)
    for s in range(steps):
        lo, hi = 2.0, 40.0
        # find the best achievable beta*: try increasing beta* until no extension has margin > 0
        cur = None
        for beta in np.arange(3.0, 40.0, 0.25):
            t0 = time.time()
            best, nres = extend(D, th, beta, rng, join=True)
            if not best or best[0][0] <= -1 / beta: break
            cur = (beta, best)
            print('  m=%d beta*=%.2f ok: margin %.4f (minRe %.4f) cand %d  %.0fs' % (m + 1, beta, best[0][0] + 1 / beta, best[0][0], nres, time.time() - t0), flush=True)
        if cur is None: print('no extension'); break
        beta, best = cur; f, D, th = best[0]; m += 1
        print('STEP m=%d beta*=%.2f minRe=%.4f theta=%s D=%s' % (m, beta, f, np.round(th, 3).tolist(), [''.join(str(int(x)) for x in r) for r in D]), flush=True)
