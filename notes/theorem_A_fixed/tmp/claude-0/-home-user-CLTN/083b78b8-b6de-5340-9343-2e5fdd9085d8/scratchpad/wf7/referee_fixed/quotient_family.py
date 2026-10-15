"""Referee's own multi-class sink-family checker (equitable-partition quotient, not Lemma P).
Core A0 (m nodes), classes j with in-sets S_j and sizes rho_j*k; classes are bidirected cliques, mutually bidirected,
send to every core node. charpoly(M) = chiQ(k) * (lam+eps)^(sum sizes - #classes). Hurwitz minors of chiQ as
polynomials in k; reports the exact threshold k* (all minors > 0 for every real k >= k*) and positivity of the
quotient solve as rational functions in k (sign on [k*, oo)).
Usage: python3 quotient_family.py FAMILY eps delta out.json"""
import sys, json
import sympy as sp
from fractions import Fraction as Fr
from reflib import *
import ctln
FAM = {'R': ('&F^U|m^QTq?', [(1,1,0,1,1,1,1)], [1]),
       'Y': ('&G^ZZZr|KjqTo', [(0,1,1,1,1,1,1,1), (1,1,0,1,1,1,1,1)], [4, 1]),
       'Q2': ('&F^U|m^QTq?', [(1,1,0,1,1,1,1), (1,1,1,0,1,1,1)], [1, 3])}
name = sys.argv[1]; eps = sp.Rational(sys.argv[2]); delta = sp.Rational(sys.argv[3]); out = sys.argv[4]
d6, Ss, rhos = FAM[name]
A0 = d6_decode(d6); m = len(A0); s = len(Ss)
k, lam = sp.symbols('k lam')
sizes = [r*k for r in rhos]
N = m + s
Q = sp.zeros(N, N)
w = -1 + eps; sgl = -1 - delta
for i in range(m):
    for j in range(m):
        Q[i, j] = -1 if i == j else (w if A0[j][i] else sgl)
    for c in range(s): Q[i, m+c] = sizes[c]*w
for c in range(s):
    for j in range(m): Q[m+c, j] = w if Ss[c][j] else sgl
    for c2 in range(s):
        Q[m+c, m+c2] = (-1 + (sizes[c]-1)*w) if c == c2 else sizes[c2]*w
chi = sp.Poly(Q.charpoly(lam).as_expr(), lam)
res = dict(family=name, eps=str(eps), delta=str(delta), legal=bool(delta*(1-eps)/eps > 1))
# identity check vs ctln at k = 1, 2
def build(kk):
    sz = [r*kk for r in rhos]; n = m + sum(sz)
    A = [[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): A[i][j] = A0[i][j]
    cls = []
    for c in range(s): cls += [c]*sz[c]
    for a in range(m, n):
        for j in range(n):
            if j != a and (j >= m): A[a][j] = 1
        for i in range(m):
            A[a][i] = 1
            if Ss[cls[a-m]][i]: A[i][a] = 1
    return A
chk = {}
for kk in [1, 2]:
    A = build(kk); n = len(A)
    W = ctln.weight_matrix(A, Fr(str(eps)), Fr(str(delta)))
    M = [[(Fr(-1) if i == j else Fr(0)) + W[i][j] for j in range(n)] for i in range(n)]
    pred = sp.Poly(chi.as_expr().subs(k, kk)*(lam+eps)**(sum(r*kk for r in rhos) - s), lam).all_coeffs()
    chk[kk] = [sp.Rational(str(c)) for c in ctln.charpoly(M)] == [sp.Rational(c) for c in pred]
res['identity_check'] = chk; print('identity', chk, flush=True)
co = [sp.expand(c) for c in chi.all_coeffs()]; d = len(co) - 1
a_ = lambda i: co[i] if 0 <= i <= d else 0
H = sp.Matrix(d, d, lambda i, j: a_(2*(i+1) - (j+1)))
kstar = 0; eventually = True; mins = []
for j in range(1, d+1):
    num, den = sp.fraction(sp.together(sp.expand(H[:j, :j].det(method='berkowitz'))))
    pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
    assert pd.degree() <= 0
    sgn = sp.sign(pd.LC())
    lc = pn.LC()*sgn
    rts = pn.real_roots() if pn.degree() > 0 else []
    big = max(rts) if rts else None
    if lc <= 0: eventually = False
    if big is not None: kstar = max(kstar, int(sp.floor(big)) + 1)
    mins.append(dict(j=j, deg=pn.degree(), LC_sign=int(sp.sign(lc)), largest_root=(float(big) if big is not None else None)))
res['minors'] = mins
res['eventually_stable'] = eventually
res['k_threshold'] = kstar if eventually else None
print('eventually stable:', eventually, 'k*', kstar, flush=True)
if eventually:
    def hz(kk):
        return ctln.hurwitz_stable([Fr(str(sp.Rational(c.subs(k, kk)))) for c in co])
    res['at_kstar'] = str(hz(kstar)); res['at_kstar_minus_1'] = str(hz(kstar-1)) if kstar >= 2 else None
    print('k*', hz(kstar), 'k*-1', hz(kstar-1) if kstar >= 2 else None, flush=True)
    # positivity: (-Q) y = 1 at symbolic k
    y = (-Q).LUsolve(sp.ones(N, 1))
    pos = True
    k0 = max(kstar, 1)
    for i in range(N):
        num, den = sp.fraction(sp.factor(sp.together(y[i])))
        pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
        rn = int(pn.count_roots(k0, None)) if pn.degree() > 0 else 0
        rd = int(pd.count_roots(k0, None)) if pd.degree() > 0 else 0
        pos &= bool(rn == 0 and rd == 0 and y[i].subs(k, k0) > 0)
    res['positive_k_ge_threshold'] = pos
    print('positive', pos)
json.dump(res, open(out, 'w'), indent=1)
