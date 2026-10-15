"""Independent derivation (own code, sympy) of the FNFL Hurwitz curve Delta_6 and of the
corner value q at beta = beta_pos.  Exact."""
import sys, json, time, pickle
import sympy as sp
from fractions import Fraction as F
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
rows = ['0011110', '0011100', '1101111', '1110111', '0011011', '0100101', '1100000']
A = [[int(c) for c in r] for r in rows]
n = 7
b, t, lam, x, q = sp.symbols('b t lam x q')
# N[i][j] = 1 iff i != j and NOT (j -> i)
N = sp.Matrix(n, n, lambda i, j: 0 if i == j else (0 if A[j][i] == 1 else 1))
# sanity vs CLAUDE.md: I - W = eps*(I + b N + t J) with eps = 1/(t+1), delta = (b-1) eps
epsS, delS = sp.symbols('eps delta')
W = sp.Matrix(n, n, lambda i, j: 0 if i == j else (-1 + epsS if A[j][i] == 1 else -1 - delS))
lhs = sp.eye(n) - W
rhs = (sp.eye(n) + b * N + t * sp.ones(n, n)) / (t + 1)
diff = (lhs - rhs).subs({epsS: 1 / (t + 1), delS: (b - 1) / (t + 1)})
assert sp.simplify(diff) == sp.zeros(n, n), 'convention mismatch'
M = (lam + 1) * sp.eye(n) + b * N
t0 = time.time()
detM = sp.expand(M.det(method='berkowitz'))
adjsum = sp.expand(sum(M.adjugate(method='berkowitz')))
p = sp.Poly(sp.expand(detM + t * adjsum), lam)
print('charpoly built', time.time() - t0, 'deg', p.degree())
a = p.all_coeffs()   # a0 .. a7, a0 = 1
assert a[0] == 1
def ent(i, j):
    k = 2 * i - j
    return a[k] if 0 <= k <= 7 else 0
H6 = sp.Matrix(6, 6, lambda i, j: ent(i + 1, j + 1))
t0 = time.time()
D6 = sp.Poly(sp.expand(H6.det(method='berkowitz')), b, t)
print('Delta6 built', time.time() - t0, 'degrees', D6.degree(b), D6.degree(t))
c, facs = sp.factor_list(D6.as_expr())
print('factors:', [(sp.Poly(f, b, t).degree(b), sp.Poly(f, b, t).degree(t), m) for f, m in facs])
pickle.dump({'p': [str(ai) for ai in a], 'D6facs': [str(f) for f, m in facs]}, open(OUT + 'fnfl_D6.pkl', 'wb'))
