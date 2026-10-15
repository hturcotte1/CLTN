"""
heur_ceiling_search.py -- HEURISTIC (floats steer only; nothing here is a certificate).
Search for cores (D, in-sets) whose sink-class limit Lbar = D - bbar 1^T has a large ceiling beta_inf = -1/min Re spec(Lbar),
with positivity (w_G(beta) > 0, sigma_j(beta) > 1) at some beta close to beta_inf.  By Theorem FP every legal point with
2 < beta < min(beta_inf, sup positivity) then carries an infinite fixed-parameter family.  Current record: Y, 4.9607.
Moves: flip one ordered pair of the core, flip one bit of an in-set, add a node.  Simulated annealing on the effective ceiling.
Usage: python3 heur_ceiling_search.py <seed> <minutes>     writes heur_ceiling_<seed>.json (best found).
"""
import sys, json, time, random
import numpy as np

def eff_ceiling(A, bs, rhos):
    m = A.shape[0]
    D = 1 - np.eye(m) - A
    bbar = sum(r * b for r, b in zip(rhos, bs)) / sum(rhos)
    L = D - np.outer(bbar, np.ones(m))
    ev = np.linalg.eigvals(L)
    mr = min(ev.real)
    if mr >= -1e-12: ceil = 50.0
    else: ceil = min(50.0, -1.0 / mr)
    if ceil <= 2.0: return 0.0, ceil
    # positivity: scan beta downward from ceil
    best = 0.0
    for beta in np.linspace(ceil - 1e-6, 2.0001, 60):
        g = 1 - 1 / beta
        try:
            w = np.linalg.solve(A.T + g * np.eye(m), np.ones(m))
        except np.linalg.LinAlgError:
            continue
        if w.min() <= 1e-12: continue
        if all((1 - b) @ w > 1 + 1e-12 for b in bs):
            best = beta; break
    return best, ceil

def random_start(m):
    A = (np.random.rand(m, m) < 0.75).astype(float); np.fill_diagonal(A, 0)
    bs = [(np.random.rand(m) < 0.2).astype(float) for _ in range(2)]
    for b in bs:
        if b.sum() == 0: b[np.random.randint(m)] = 1
    return A, bs

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = np.zeros((n, n)); bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i, j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A

def main(seed, minutes):
    random.seed(seed); np.random.seed(seed)
    if seed % 2 == 0:
        A = from_d6('&G^ZZZr|KjqTo'); m = 8
        bs = [np.zeros(m), np.zeros(m)]; bs[0][0] = 1; bs[1][2] = 1; rhos = [4, 1]
    else:
        A, bs = random_start(9); rhos = [1, 1]
    cur, _ = eff_ceiling(A, bs, rhos); best = (cur, A.copy(), [b.copy() for b in bs], list(rhos))
    T = 0.05; t0 = time.time(); it = 0
    while time.time() - t0 < 60 * minutes:
        it += 1
        A2 = A.copy(); bs2 = [b.copy() for b in bs]; rh2 = list(rhos)
        mv = random.random(); m = A.shape[0]
        if mv < 0.6:
            i, j = random.sample(range(m), 2); A2[i, j] = 1 - A2[i, j]
        elif mv < 0.85:
            jj = random.randrange(len(bs2)); i = random.randrange(m); bs2[jj][i] = 1 - bs2[jj][i]
            if bs2[jj].sum() == 0: continue
        elif mv < 0.95:
            rh2[random.randrange(len(rh2))] = random.randint(1, 6)
        elif m < 12:
            # add a node: random row/col
            A3 = np.zeros((m + 1, m + 1)); A3[:m, :m] = A2
            A3[m, :m] = (np.random.rand(m) < 0.8); A3[:m, m] = (np.random.rand(m) < 0.8)
            A2 = A3; bs2 = [np.append(b, 0.0) for b in bs2]
        else:
            continue
        val, _ = eff_ceiling(A2, bs2, rh2)
        if val >= cur or random.random() < np.exp((val - cur) / T):
            A, bs, rhos, cur = A2, bs2, rh2, val
            if cur > best[0] + 1e-9:
                best = (cur, A.copy(), [b.copy() for b in bs], list(rhos))
                print('it %d  m %d  eff ceiling %.4f' % (it, A.shape[0], cur), flush=True)
        T = max(0.002, T * 0.9999)
    val, A, bs, rhos = best
    json.dump({'seed': seed, 'eff_ceiling_float': val, 'A': A.astype(int).tolist(), 'bs': [b.astype(int).tolist() for b in bs],
               'rhos': rhos, 'note': 'HEURISTIC search output; not a certificate'}, open('heur_ceiling_%d.json' % seed, 'w'))
    print('best', val)

if __name__ == '__main__':
    main(int(sys.argv[1]), float(sys.argv[2]))
