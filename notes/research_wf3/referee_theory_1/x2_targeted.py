"""Targeted X2/X2' test: bias toward positive base graphs (so the criterion is non-trivially exercised on both sides).
Compares the k-free criterion against ctln.analyze_support positivity on explicit graphs. Exact only."""
import sys, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp, ctln
random.seed(int(sys.argv[1]))
def w_vec(A, gamma):
    m = len(A); P = sp.Matrix(A).T + sp.Rational(gamma) * sp.eye(m)
    return None if P.det() == 0 else P.LUsolve(sp.ones(m, 1))
def ext_multi(A, classes):
    m = len(A); n = m + sum(k for S, k in classes); B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    owner = [c for c, (S, k) in enumerate(classes) for _ in range(k)]
    for v in range(m, n):
        S = classes[owner[v - m]][0]
        for i in range(m): B[i][v] = S[i]; B[v][i] = 1
        for u in range(m, n):
            if u != v: B[v][u] = 1
    return B
stats = {'cases': 0, 'pos': 0, 'crit_true': 0, 'wG_pos_but_sigma_le_1': 0, 'fail': 0}
for trial in range(400):
    m = random.randint(3, 6)
    A = [[0 if i == j else (1 if random.random() < 0.8 else 0) for j in range(m)] for i in range(m)]
    eps = F(1, random.randint(2, 30)); beta = 1 + F(random.randint(1, 40), 20); delta = (beta - 1) * eps
    gamma = 1 - 1 / beta
    w = w_vec(A, gamma)
    if w is None or not all(v > 0 for v in w): continue
    for _ in range(3):
        nc = random.randint(1, 3)
        classes = [([random.choice([0, 1, 1, 1]) for _ in range(m)], random.randint(1, 3)) for _ in range(nc)]
        sig = [sum((w[i] for i in range(m) if S[i]), sp.Integer(0)) for S, k in classes]
        crit = all(s > 1 for s in sig)
        if not crit: stats['wG_pos_but_sigma_le_1'] += 1
        B = ext_multi(A, classes)
        r = ctln.analyze_support(B, list(range(len(B))), eps, delta)
        if r['x'] is None: continue
        stats['cases'] += 1; stats['pos'] += r['positive']; stats['crit_true'] += crit
        if r['positive'] != crit:
            stats['fail'] += 1; print('FAIL', A, classes, eps, delta, r['positive'], crit)
print(stats)
