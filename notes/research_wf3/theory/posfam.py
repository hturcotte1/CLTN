"""Exact positivity data as polynomials in beta.
For integer beta:  Mb := beta A^T + (beta - 1) I = beta (A^T + gamma I).  w = beta * adj(Mb) 1 / det(Mb).
adj(Mb)1 and det(Mb) are integer polynomials in beta of degree <= m-1 and m; obtained by exact interpolation
with check points (fastregion.idet / cramer_ones)."""
import sys
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
from fractions import Fraction as F
from fastregion import idet, cramer_ones, interp, to_int, peval
bsym = sp.symbols('beta')
def w_polys(A):
    m = len(A)
    xs = list(range(2, m + 5))
    dets, adjs = [], []
    for x in xs:
        Mb = [[x * A[j][i] + ((x - 1) if i == j else 0) for j in range(m)] for i in range(m)]
        dets.append(idet(Mb)); adjs.append(cramer_ones(Mb))
    def fit(ys, deg):
        p = to_int(interp(xs[:deg + 1], ys[:deg + 1]))
        assert all(peval(p, xx) == yy for xx, yy in zip(xs, ys)), 'interpolation check failed'
        return sp.Poly(list(reversed(p)), bsym)
    d = fit(dets, m)
    a = [fit([v[i] for v in adjs], m - 1) for i in range(m)]
    return d, a     # w_i = beta * a_i / d
def sigma_minus_1(A, S):
    d, a = w_polys(A)
    num = sp.Poly(bsym, bsym) * sum((a[i] for i in range(len(A)) if S[i]), sp.Poly(0, bsym)) - d
    return num, d   # sigma_S - 1 = num / d
