"""Independent all-k proof for the self-join R^{vk} (n = 7k) at fixed (eps, delta).
W = I_k (x) W_R + (J_k - I_k) (x) (eps-1) J_7. Diagonalising J_k:
  charpoly(-I+W) = charpoly(B0 + (k-1)(eps-1)J) * charpoly(B0 + (1-eps)J)^(k-1),  B0 = -I + W_R.
Positivity: x = 1 (x) y, (I - W_R + (k-1)(1-eps)J) y = 1.
Usage: python3 proof_SJ.py eps delta out.json"""
import sys, json
import sympy as sp
from fractions import Fraction as Fr
from reflib import *
import ctln
eps = sp.Rational(sys.argv[1]); delta = sp.Rational(sys.argv[2]); out = sys.argv[3]
k, lam = sp.symbols('k lam')
A0 = d6_decode(R_D6); m = 7
def Wc(i, j):
    if i == j: return 0
    return -1 + eps if A0[j][i] == 1 else -1 - delta
B0 = sp.Matrix(m, m, lambda i, j: (-1 if i == j else 0) + Wc(i, j))
J = sp.ones(m, m)
P1 = sp.Poly((B0 + (k-1)*(eps-1)*J).charpoly(lam).as_expr(), lam)
P2 = sp.Poly((B0 + (1-eps)*J).charpoly(lam).as_expr(), lam)
res = {'eps': str(eps), 'delta': str(delta), 'legal': bool(delta*(1-eps)/eps > 1)}
chk = {}
for kk in [1, 2, 3, 4]:
    A = selfjoin(A0, kk); n = len(A)
    W = ctln.weight_matrix(A, Fr(str(eps)), Fr(str(delta)))
    M = [[(Fr(-1) if i == j else Fr(0)) + W[i][j] for j in range(n)] for i in range(n)]
    cpM = ctln.charpoly(M)
    pred = sp.Poly(P1.as_expr().subs(k, kk) * P2.as_expr()**(kk-1), lam).all_coeffs()
    chk[kk] = [sp.Rational(str(c)) for c in cpM] == [sp.Rational(c) for c in pred]
res['identity_check'] = chk; print('identity', chk)
def hurwitz_minors(coeffs):
    d = len(coeffs) - 1
    a = lambda i: coeffs[i] if 0 <= i <= d else 0
    H = sp.Matrix(d, d, lambda i, j: a(2*(i+1) - (j+1)))
    return [sp.factor(sp.expand(H[:j, :j].det(method='berkowitz'))) for j in range(1, d+1)]
# fixed factor
st2 = ctln.hurwitz_stable([Fr(str(c)) for c in P2.all_coeffs()])
res['fixed_factor_hurwitz'] = str(st2); print('fixed factor', st2)
ok = st2[0] is True
mins = []
for j, mj in enumerate(hurwitz_minors([sp.expand(c) for c in P1.all_coeffs()]), 1):
    num, den = sp.fraction(sp.together(mj))
    pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
    rn = int(pn.count_roots(1, None)) if pn.degree() > 0 else 0
    rd = int(pd.count_roots(1, None)) if pd.degree() > 0 else 0
    good = rn == 0 and rd == 0 and mj.subs(k, 1) > 0
    ok &= bool(good); mins.append(dict(j=j, num=str(pn.as_expr()), den=str(pd.as_expr()), positive_k_ge_1=bool(good)))
    print('minor', j, good)
res['minors_k_factor'] = mins
# Sherman-Morrison: (B + c 1 1^T)^{-1} 1 = B^{-1} 1 / (1 + c 1^T B^{-1} 1), B = I - W_R, c = (k-1)(1-eps) >= 0
B = sp.eye(m) - (B0 + sp.eye(m))
v = B.LUsolve(sp.ones(m, 1))
ypos = all(vi > 0 for vi in v)          # then denominator 1 + c*sum(v) > 0 and y > 0 for all k >= 1
yk = (B + 2*(1-eps)*J).LUsolve(sp.ones(m, 1))
assert list(yk) == [vi/(1 + 2*(1-eps)*sum(v)) for vi in v]
res['y'] = [str(vi) for vi in v]
res['positive_all_k'] = bool(ypos)
res['VERDICT_SFP_all_k_ge_1'] = bool(ok and ypos and res['legal'] and all(chk.values()))
print('VERDICT', res['VERDICT_SFP_all_k_ge_1'])
json.dump(res, open(out, 'w'), indent=1)
