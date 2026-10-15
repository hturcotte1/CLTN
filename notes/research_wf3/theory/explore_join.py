import numpy as np, load_hits, tlib
H6 = load_hits.n6_hits()
named = {'G1': H6[1000255], 'H3': H6[1387528], 'G2': H6[1424629]}
def join(A, k):
    m = len(A); n = m + k
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if i < m and j < m: B[i][j] = A[i][j]
            else: B[i][j] = 1
    return B
for nm, A in named.items():
    m = len(A)
    D = np.array([[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)], float)
    ev = np.linalg.eigvals(D)
    print(nm, 'spec(D) =', np.round(sorted(ev, key=lambda z: z.real), 4), 'min Re =', round(min(ev.real), 5), ' -> beta_inf =', round(-1/min(ev.real), 5) if min(ev.real) < 0 else 'inf')
    for k in range(0, 12):
        B = join(A, k)
        a = tlib.heur_alpha(B)
        print('   k=%d alpha=%.5f beta_c=%.5f' % (k, a, 1/(1+a)))
