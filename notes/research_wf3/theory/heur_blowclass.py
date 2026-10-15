# HEURISTIC (floats): for a base graph G and a new node class v (k closed twins, bidirected among themselves)
# with in-set S (who sends to v) and out-set T (whom v sends to), evaluate the limit-stable/positive beta range
# for k = 1..K and the k -> infinity limit spectrum spec(D - b 1^T - 1 c^T).
import numpy as np, itertools, sys, load_hits, tlib
def blow(A, S, T, k):
    m = len(A); n = m + k
    B = np.zeros((n, n), int)
    B[:m, :m] = A
    for i in range(m):
        if S[i]: B[i, m:] = 1
        if T[i]: B[m:, i] = 1
    for a in range(m, n):
        for b2 in range(m, n):
            if a != b2: B[a, b2] = 1
    return B
def alpha_np(B):
    n = len(B)
    c1 = np.poly(B.astype(float)); c2 = np.poly((B + 1).astype(float))
    p = (c1 - c2)[1:]
    return max(np.roots(p).real)
def pos_np(B, beta):
    g = 1 - 1/beta
    w = np.linalg.solve(B.T + g*np.eye(len(B)), np.ones(len(B)))
    return w.min() > 1e-12
def best_beta(B, grid):
    a = alpha_np(B)
    if a >= -0.5: return None
    bc = 1/(1+a)
    ok = [b for b in grid if b < bc and pos_np(B, b)]
    return max(ok) if ok else None
if __name__ == '__main__':
    H6 = load_hits.n6_hits()
    named = {'G1': H6[1000255], 'H3': H6[1387528]}
    grid = np.linspace(2.0005, 3.0, 400)
    for nm, A in named.items():
        A = np.array(A); m = len(A)
        D = 1 - np.eye(m, dtype=int) - A
        for S in itertools.product([0,1], repeat=m):
            for T in itertools.product([0,1], repeat=m):
                bs = [best_beta(blow(A, S, T, k), grid) for k in (1, 2, 3, 5, 8)]
                if sum(x is not None for x in bs) >= 3:
                    b = 1 - np.array(S); c = 1 - np.array(T)
                    L = D - np.outer(b, np.ones(m)) - np.outer(np.ones(m), c)
                    ev = np.linalg.eigvals(L)
                    print(nm, S, T, [None if x is None else round(x, 3) for x in bs], 'limit minRe(spec(D-b1-1c)) = %.4f' % min(ev.real))
