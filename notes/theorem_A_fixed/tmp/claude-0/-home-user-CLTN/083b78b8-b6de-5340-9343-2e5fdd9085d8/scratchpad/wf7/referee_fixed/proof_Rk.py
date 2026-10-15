"""Independent all-k proof for R_k at a fixed (eps, delta), via the equitable-partition
quotient of M = -I + W (NOT via the prover's Lemma P / Psi / D-picture).
  M P = P Q, P = indicator of cells {0},...,{6},{class}; complement (class vectors summing to 0)
  is M-invariant with eigenvalue -eps. Hence charpoly(M) = charpoly(Q) (lam+eps)^(k-1), k >= 1.
Then Hurwitz minors of charpoly(Q) are polynomials in k; show > 0 on [1, oo) exactly (Sturm).
Positivity: x = P y with (-Q) y = 1; check sign of y as rational functions of k on [1, oo).
Usage: python3 proof_Rk.py eps delta out.json"""
import sys, json
import sympy as sp
from fractions import Fraction as Fr
from reflib import *
import ctln

eps = sp.Rational(sys.argv[1]); delta = sp.Rational(sys.argv[2]); out = sys.argv[3]
k, lam = sp.symbols('k lam')
A0 = d6_decode(R_D6); S = S_R; m = 7

def Wc(i, j):  # influence of j on i in the core
    if i == j: return 0
    return -1 + eps if A0[j][i] == 1 else -1 - delta

Q = sp.zeros(m+1, m+1)
for i in range(m):
    for j in range(m):
        Q[i, j] = (-1 if i == j else 0) + Wc(i, j)
    Q[i, m] = k*(-1 + eps)
for j in range(m):
    Q[m, j] = (-1 + eps) if S[j] else (-1 - delta)
Q[m, m] = -1 + (k-1)*(-1 + eps)

cp = sp.Poly(Q.charpoly(lam).as_expr(), lam)
coeffs = [sp.expand(c) for c in cp.all_coeffs()]
res = {'eps': str(eps), 'delta': str(delta), 'legal': bool(delta*(1-eps)/eps > 1)}

# cross-check identity against ctln.charpoly on the explicit graph, k = 1..6 and k = 30
chk = {}
for kk in [1, 2, 3, 4, 5, 6, 30]:
    A = build_Rk(kk); n = len(A)
    W = ctln.weight_matrix(A, Fr(str(eps)), Fr(str(delta)))
    M = [[(Fr(-1) if i == j else Fr(0)) + W[i][j] for j in range(n)] for i in range(n)]
    cpM = ctln.charpoly(M)
    pred = sp.Poly(cp.as_expr().subs(k, kk) * (lam + eps)**(kk-1), lam).all_coeffs()
    chk[kk] = [sp.Rational(str(c)) for c in cpM] == [sp.Rational(c) for c in pred]
res['identity_check'] = chk
print('identity', chk, flush=True)

# Hurwitz minors in k
d = len(coeffs) - 1
def a(i): return coeffs[i] if 0 <= i <= d else 0
H = sp.Matrix(d, d, lambda i, j: a(2*(i+1) - (j+1)))
minors = []
allpos = True
for j in range(1, d+1):
    mj = sp.factor(sp.expand(H[:j, :j].det(method='berkowitz')))
    num, den = sp.fraction(sp.together(mj))
    pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
    # sign on [1, oo): no real root in [1, oo) and positive at k=1
    rn = int(pn.count_roots(1, None)) if pn.degree() > 0 else 0
    rd = int(pd.count_roots(1, None)) if pd.degree() > 0 else 0
    v1 = mj.subs(k, 1)
    # leading coefficient sign (for k -> oo)
    pos = (rn == 0 and rd == 0 and v1 > 0)
    # also check real k in [0,1)
    v0 = mj.subs(k, 0)
    rn0 = pn.count_roots(0, 1) if pn.degree() > 0 else 0
    allpos &= pos
    coefpos = all(c > 0 for c in pn.all_coeffs()) and all(c > 0 for c in pd.all_coeffs()) or \
              all(c < 0 for c in pn.all_coeffs()) and all(c < 0 for c in pd.all_coeffs())
    minors.append(dict(j=j, deg_k=pn.degree(), roots_k_ge_1=rn, den_roots=rd,
                       value_k1=str(v1), positive_on_k_ge_1=bool(pos), all_coeffs_same_sign=bool(coefpos),
                       num=str(pn.as_expr()), den=str(pd.as_expr())))
    print('minor', j, 'deg', pn.degree(), 'roots>=1:', rn, 'den roots:', rd, 'v1>0', v1 > 0, 'coefpos', coefpos, flush=True)
res['charpoly_Q_coeffs'] = [str(c) for c in coeffs]
res['minors'] = minors
res['stable_all_k_ge_1'] = bool(allpos and coeffs[0] == 1)

# positivity
y = (-Q).LUsolve(sp.ones(m+1, 1))
ypos = True; yinfo = []
for i in range(m+1):
    yi = sp.factor(sp.together(y[i]))
    num, den = sp.fraction(yi)
    pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
    rn = int(pn.count_roots(1, None)) if pn.degree() > 0 else 0
    rd = int(pd.count_roots(1, None)) if pd.degree() > 0 else 0
    ok = rn == 0 and rd == 0 and yi.subs(k, 1) > 0
    ypos &= ok
    yinfo.append(dict(i=i, y=str(yi), positive_on_k_ge_1=bool(ok)))
    print('y', i, yi, ok, flush=True)
res['y'] = yinfo
res['positive_all_k_ge_1'] = bool(ypos)
res['VERDICT_SFP_all_k_ge_1'] = bool(allpos and ypos and res['legal'] and all(chk.values()))
print('VERDICT', res['VERDICT_SFP_all_k_ge_1'])
json.dump(res, open(out, 'w'), indent=1)
