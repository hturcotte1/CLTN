"""EXPLORATION (floats): constant D-in-degree d cores. Then u = 1/(d+r) > 0 at every r, and the
limit criterion is: all eigenvalues of D other than the Perron value d have Re > -1/beta.
Search: annealing over column-wise in-neighbour sets, maximize min Re of non-Perron eigenvalues."""
import numpy as np, sys
def score(D, d):
    ev = np.linalg.eigvals(D)
    i = np.argmin(abs(ev - d)); ev = np.delete(ev, i)
    return ev.real.min()
def rand_D(m, d, rng):
    D = np.zeros((m, m))
    for j in range(m):
        ch = rng.choice([i for i in range(m) if i != j], d, replace=False); D[ch, j] = 1
    return D
def anneal(m, d, rng, iters=20000):
    D = rand_D(m, d, rng); s = score(D, d); best = (s, D.copy()); T = 0.05
    for it in range(iters):
        j = rng.integers(m); ins = np.flatnonzero(D[:, j]); outs = [i for i in range(m) if i != j and D[i, j] == 0]
        a = rng.choice(ins); b = rng.choice(outs)
        E = D.copy(); E[a, j] = 0; E[b, j] = 1
        t = score(E, d)
        if t > s or rng.random() < np.exp((t - s) / T):
            D, s = E, t
            if s > best[0]: best = (s, D.copy())
        T *= 0.9997
    return best
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    for d in (1, 2, 3):
        for m in (5, 6, 7, 8, 10, 12):
            if m <= d + 1: continue
            b = max((anneal(m, d, rng, 6000) for _ in range(3)), key=lambda x: x[0])
            print(d, m, round(b[0], 4), flush=True)
