"""Independent check of the FNFL supremum: own Delta_6 factor, own positivity boundary,
own resultant; exact only (floats only in printed labels)."""
import pickle, time
import sympy as sp
from fractions import Fraction as F
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
b, t, lam, x, q = sp.symbols('b t lam x q')
d = pickle.load(open(OUT + 'fnfl_D6.pkl', 'rb'))
facs = [sp.Poly(sp.sympify(s), b, t) for s in d['D6facs']]
f0, f1 = sorted(facs, key=lambda P: P.degree(t))
rows = ['0011110', '0011100', '1101111', '1110111', '0011011', '0100101', '1100000']
A = [[int(c) for c in r] for r in rows]
n = 7
N = sp.Matrix(n, n, lambda i, j: 0 if i == j else (0 if A[j][i] == 1 else 1))
C = sp.eye(n) + b * N
adjC = C.adjugate(method='berkowitz')
a_num = [sp.factor(sp.expand(sum(adjC[i, j] for j in range(n)))) for i in range(n)]
detC = sp.factor(C.det(method='berkowitz'))
s1 = sp.factor(sp.expand(sum(adjC)))
print('adj(C)1 entries:', a_num)
print('det C =', detC, ' 1^T adj C 1 =', s1)
m = sp.Poly(x**5 - 3*x**4 + 3*x**3 - 3*x**2 + 2*x - 1, x)
ivm = m.intervals(eps=sp.Rational(1, 10**40), inf=2, sup=3)
print('m real roots in [2,3]:', len(ivm))
bpl, bph = ivm[0][0]
# which entries / D change sign at beta_pos ?
for k, e in enumerate(a_num + [detC, s1]):
    Pe = sp.Poly(e, b)
    print(k, 'divisible by m:', sp.rem(Pe, sp.Poly(m.as_expr().subs(x, b), b)).is_zero, 'sign just below/above beta_pos:',
          sp.sign(Pe.eval(bpl)), sp.sign(Pe.eval(bph)), 'sign at 2.01:', sp.sign(Pe.eval(sp.Rational(201, 100))))
# t-roots of f1 and f0 at b0 = 53/26
b0 = sp.Rational(53, 26)
for nm, f in (('f0', f0), ('f1', f1)):
    u = sp.Poly(f.as_expr().subs(b, b0), t)
    print(nm, 'roots t>26 at b0:', [(str(i[0][0]), str(i[0][1])) for i in u.intervals(eps=sp.Rational(1, 100), inf=26)])
cert_tpoly = sp.Poly([11692481028755823439445120, -609495408176817336149678976, -944295302658237141597537344, 257170292685528203722425880, 766315152811392722530354986, 133579293167248823762723449], x)
u1 = sp.Poly(f1.as_expr().subs(b, b0).subs(t, x), x)
print('certificate t_critical poly divides f1(b0,.):', sp.rem(u1, cert_tpoly).is_zero)
# G(beta, q) from f1
dd = f1.degree(t)
fe = f1.as_expr()
coeffs = [sp.expand(fe).coeff(t, k) for k in range(dd + 1)]
G = sp.expand(sum(coeffs[k] * q**k * (b - 1 - q)**(dd - k) for k in range(dd + 1)))
t0 = time.time()
R = sp.Poly(sp.resultant(G, m.as_expr().subs(x, b), b), q)
print('Res_beta(G,m) degree', R.degree(), time.time() - t0)
claimed = sp.Poly([40826503094343665477600174971, -753366531827304304481324176101, 6920794832000172401966938985993, -39605390504173099485395318553433, 155984378007350311862036281188332, -450035732605417009056688726739125, 992591101455307379568685457222638, -1724565669721312192494727158366384, 2411178845096551677587105674831221, -2753753811325699214759131879769065, 2594983763217058248472189972075213, -2029630423563061369032505400144017, 1320138591083831496294796327256552, -712390162299956463405343776194497, 316405081054021495028562333521936, -113713138681834205966504859913148, 31919860709284149283091464529960, -6415526052269403098847283704038, 645832494412960041944283180166, 105975949944847104525916618240, -65845669129385166875179235187, 16329682753834842966110031267, -2560426304908442501475975573, 263227361372335465557317601, -16413519949660025315455046, 475439792660147674691197], q)
print('claimed irreducible:', claimed.is_irreducible, 'degree', claimed.degree())
print('claimed divides Res_beta(G,m):', sp.rem(R, claimed).is_zero)
qa, qb = sp.Rational(685003, 657781), sp.Rational(2170535, 2084278)
print('claimed roots in [qa,qb]:', claimed.count_roots(qa, qb), ' width', sp.N(abs(qb - qa)), 'qa<qb', qa < qb)
# independent: t-roots of f1(beta_pos, t) via Res_beta(f1, m) in t, and their q values
Rt = sp.Poly(sp.resultant(f1.as_expr(), m.as_expr().subs(x, b), b), t)
ivt = Rt.intervals(eps=sp.Rational(1, 10**12), inf=0)
print('t-roots (positive) of f1(beta_pos, t):', [sp.N((i[0][0] + i[0][1]) / 2, 12) for i in ivt])
cq = [sp.N((bpl - 1) * ((i[0][0] + i[0][1]) / 2) / ((i[0][0] + i[0][1]) / 2 + 1), 15) for i in ivt]
print('q at those corners (label only):', cq)
pickle.dump({'f0': str(f0.as_expr()), 'f1': str(f1.as_expr()), 'bpl': str(bpl), 'bph': str(bph)}, open(OUT + 'fnfl_f.pkl', 'wb'))
