"""Independent check of the (1/2, 2) instance (prover's Theorem Lambda, E4): the 65-node arms-forest core C (stored d6,
decoded by the referee's own decoder), single sink class with in-set = all but node 2 (c2), G_k = C[S, all; k].
Method (referee's own, not Lemma P): equitable-partition quotient Q of M = -I + W (66x66); its last column is affine in k
and the rest constant, so det(lam I - Q) is affine in k:  chiQ(k) = chiQ(0) + k (chiQ(1) - chiQ(0)) (formal k=0 quotient).
charpoly(M) = chiQ(k) (lam + eps)^(k-1). Positivity via Sherman-Morrison on the affine column.
Writes check_lambda_instance.json."""
import json, sys, time
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import ctln

def d6_decode_any(s):
    assert s[0] == '&'
    body = [ord(c) - 63 for c in s[1:]]
    if body[0] == 63:
        assert body[1] != 63
        n = (body[1] << 12) | (body[2] << 6) | body[3]; rest = body[4:]
    else:
        n = body[0]; rest = body[1:]
    bits = []
    for v in rest:
        for b in range(5, -1, -1): bits.append((v >> b) & 1)
    A = [[bits[i*n+j] for j in range(n)] for i in range(n)]
    assert all(A[i][i] == 0 for i in range(n))
    return A

pf = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf7/prover_fixed/extend_lambda.json'
A0 = d6_decode_any(json.load(open(pf))['core_d6'])
m = len(A0); print('m', m)
S = [0 if i == 2 else 1 for i in range(m)]
eps, delta = Fr(1, 2), Fr(2)
assert ctln.legal(eps, delta)
t0 = time.time()

def Q_of(kk):
    Q = [[Fr(0)]*(m+1) for _ in range(m+1)]
    for i in range(m):
        for j in range(m):
            Q[i][j] = Fr(-1) if i == j else ((-1 + eps) if A0[j][i] else (-1 - delta))
        Q[i][m] = kk*(-1 + eps)
    for j in range(m):
        Q[m][j] = (-1 + eps) if S[j] else (-1 - delta)
    Q[m][m] = -1 + (kk - 1)*(-1 + eps)
    return Q

c0 = ctln.charpoly(Q_of(0)); c1 = ctln.charpoly(Q_of(1)); c2 = ctln.charpoly(Q_of(2))
assert [2*b - a for a, b in zip(c0, c1)] == c2, 'affine-in-k check failed'
print('charpolys', round(time.time()-t0, 1), flush=True)
k, lam = sp.symbols('k lam')
chi = sp.Poly(sum((sp.Rational(a.numerator, a.denominator) + k*sp.Rational((b-a).numerator, (b-a).denominator))*lam**(m+1-i)
                  for i, (a, b) in enumerate(zip(c0, c1))), lam)
e = sp.Rational(1, 2)
# strip (lam + eps) factors
j0 = 0; g = chi
while True:
    qq, rr = sp.div(g, sp.Poly(lam + e, lam))
    if rr.is_zero: g = qq; j0 += 1
    else: break
print('factor (lam+eps)^', j0, 'remaining degree', g.degree(), flush=True)
co = [sp.expand(c) for c in g.all_coeffs()]
d = len(co) - 1
a_ = lambda i: co[i] if 0 <= i <= d else 0
H = sp.Matrix(d, d, lambda i, j: a_(2*(i+1) - (j+1)))
mins = []; kstar = 0
for j in range(1, d+1):
    mj = sp.together(sp.expand(H[:j, :j].det(method='berkowitz')))
    num, den = sp.fraction(mj)
    pn = sp.Poly(sp.expand(num), k); pd = sp.Poly(sp.expand(den), k)
    assert pd.degree() <= 0 and pd.LC() > 0
    rts = [r for r in pn.real_roots()]
    big = max([r for r in rts], default=None)
    lc = pn.LC()
    mins.append(dict(j=j, deg=pn.degree(), LC_positive=bool(lc > 0), largest_real_root=str(big),
                     largest_real_root_float_display=float(big) if big is not None else None))
    if big is not None: kstar = max(kstar, int(sp.floor(big)) + 1)
    assert lc > 0
    print('minor', j, 'deg', pn.degree(), 'LC>0', lc > 0, 'largest root ~', float(big) if big is not None else None, flush=True)
# threshold: all minors > 0 for every real k > max root; check integers kstar-1 and kstar exactly
def stable_at(kk):
    cs = [Fr(str(sp.Rational(c.subs(k, kk)))) for c in co]
    return ctln.hurwitz_stable(cs)
res = dict(eps=str(eps), delta=str(delta), m=m, lam_plus_eps_power=j0, g_degree=d, minors=mins,
           k_threshold=kstar, at_kstar=str(stable_at(kstar)), at_kstar_minus_1=str(stable_at(kstar-1)))
print('kstar', kstar, stable_at(kstar), stable_at(kstar-1))
# positivity: (-Q(k)) y = 1, -Q(k) = B + k u e_m^T
B = [[-x for x in row] for row in Q_of(0)]
u = [-(Q_of(1)[i][m] - Q_of(0)[i][m]) for i in range(m+1)]
Binv1 = ctln.solve_exact(B, [Fr(1)]*(m+1)); Binvu = ctln.solve_exact(B, u)
alpha, beta_ = Binv1[m], Binvu[m]
# y_i(k) = [Binv1_i (1 + k beta) - k Binvu_i alpha] / (1 + k beta)
def sign_ok_from(k0, p0, p1):   # p0 + p1 k > 0 for all real k >= k0
    return (p0 + p1*k0 > 0) and (p1 >= 0)
den_ok = sign_ok_from(1, Fr(1), beta_)
num_ok = all(sign_ok_from(1, Binv1[i], Binv1[i]*beta_ - Binvu[i]*alpha) for i in range(m+1))
# sanity: compare with a direct solve at k=3
y3 = ctln.solve_exact([[-x for x in row] for row in Q_of(3)], [Fr(1)]*(m+1))
pred3 = [(Binv1[i]*(1+3*beta_) - 3*Binvu[i]*alpha)/(1+3*beta_) for i in range(m+1)]
res.update(positivity_den_ok=den_ok, positivity_num_ok=num_ok, SM_check_k3=(y3 == pred3))
print('positivity k>=1:', den_ok, num_ok, 'SM check', y3 == pred3)
res['VERDICT'] = (f'SFP for every integer k >= {kstar}' if (den_ok and num_ok and y3 == pred3) else 'FAIL')
res['seconds'] = round(time.time() - t0, 1)
json.dump(res, open('check_lambda_instance.json', 'w'), indent=1)
print(res['VERDICT'])
