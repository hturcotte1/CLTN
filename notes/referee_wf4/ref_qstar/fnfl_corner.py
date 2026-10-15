"""Exact enclosure of the corner value q(beta_pos, s(beta_pos)) via interval arithmetic on
rational enclosures of beta_pos and of the t-root of Res_beta(f1, m) near 25."""
import pickle
import sympy as sp
from fractions import Fraction as F
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
b, t, x = sp.symbols('b t x')
d = pickle.load(open(OUT + 'fnfl_f.pkl', 'rb'))
f1 = sp.Poly(sp.sympify(d['f1']), b, t)
m = b**5 - 3*b**4 + 3*b**3 - 3*b**2 + 2*b - 1
bi = sp.Poly(m, b).intervals(eps=sp.Rational(1, 10**40), inf=2, sup=3)[0][0]
Rt = sp.Poly(sp.resultant(f1.as_expr(), m, b), t)
ti = [i[0] for i in Rt.intervals(eps=sp.Rational(1, 10**40), inf=13, sup=40)]
assert len(ti) == 1, ti   # unique t-root of Res in [13,40] (legal side at beta_pos: t > 12.05)
# confirm this t-root is a root of f1(beta_pos, .): f1(b,t) changes sign on the box? use exact check:
# g(t) = Res_beta(f1, m) root -> some conjugate beta of m; beta_pos is m's only real root, so the common
# root with real t must have real beta = beta_pos (f1 real coefficients, t real => conjugate pair issue):
# verify directly that f1(beta, t) has a sign change across the t-interval at both beta-enclosure ends.
bl, bh = bi; tl, th = ti[0]
for bb in (bl, bh):
    v1 = f1.as_expr().subs({b: bb, t: tl}); v2 = f1.as_expr().subs({b: bb, t: th})
    print('sign change in t at beta end', sp.sign(v1), sp.sign(v2))
Fr = lambda r: F(int(r.p), int(r.q))
bl, bh, tl, th = map(Fr, (bl, bh, tl, th))
ql = (bl - 1) * tl / (tl + 1); qh = (bh - 1) * th / (th + 1)    # q increasing in both
QA, QB = F(685003, 657781), F(2170535, 2084278)
print('corner q enclosure inside (QA,QB):', QA < ql and qh < QB, ' width', float(qh - ql))
# --- exact: s(beta_pos) is the unique root of f1(beta_pos, .) in [12, 40] and lies in [tl, th]
print('Res roots in [12,40]:', Rt.count_roots(12, 40), ' t_min(beta_pos) < 12.1:', 1 / (bl - 2) < F(121, 10))
mpoly = sp.Poly(m, b)
def sign_at_beta_pos(tt):
    g = sp.Poly(f1.as_expr().subs(t, sp.Rational(tt.numerator, tt.denominator)), b)
    lo, hi = sp.Rational(bl.numerator, bl.denominator), sp.Rational(bh.numerator, bh.denominator)
    # shrink the beta_pos enclosure until g has no root in it (g(beta_pos) != 0 since t-root is irrational-free check)
    ivb = mpoly.intervals(eps=sp.Rational(1, 10**200), inf=2, sup=3)[0][0]
    assert g.count_roots(ivb[0], ivb[1]) == 0
    return sp.sign(g.eval(ivb[0]))
print('sign f1(beta_pos, tl), f1(beta_pos, th):', sign_at_beta_pos(tl), sign_at_beta_pos(th))
