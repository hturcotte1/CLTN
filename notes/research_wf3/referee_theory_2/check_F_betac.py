"""Referee: beta_c(k) table, plastic plateau k = 2..8, and monotone decrease for k >= 8.
Exact: limit-stable at beta  iff  R_k(-1/beta - x) Hurwitz (Lemma D1 + I.3(b); zero roots of z^(k+1) have Re 0 > -1/beta).
Floats only for display."""
import sympy as sp
from fractions import Fraction as Fr
z, k, x, r = sp.symbols('z k x r')
Rk = (6 + k)*z**4 + (11 + k)*z**3 + (16 + k)*z**2 + 10*z + 3
co = sp.Poly(sp.expand(Rk.subs(z, -r - x)), x).all_coeffs()
a0, a1, a2, a3, a4 = co
D2 = sp.expand(a1*a2 - a0*a3); D3 = sp.expand(a3*D2 - a1**2*a4)
minors = [a0, a1, D2, D3, a4]
def hurwitz(kk, rr):
    vals = [m.subs({k: kk, r: rr}) for m in minors]
    if any(v == 0 for v in vals): return None
    return all(v > 0 for v in vals)
table = {0: '2.16635502785', 1: '2.27474704987', 2: '2.33914467333', 3: '2.37197201365', 4: '2.38331632949', 5: '2.38085323090',
         6: '2.37010799260', 7: '2.35485699604', 8: '2.33755471291', 9: '2.31971710694', 10: '2.30222853360', 11: '2.28556715757', 12: '2.26996011736'}
for kk, s in table.items():
    bc = Fr(s); lo = bc - Fr(1, 10**9); hi = bc + Fr(1, 10**9)
    # stable on a grid between 2 and lo (sanity), at lo; not stable at hi
    assert hurwitz(kk, sp.Rational(lo.denominator, lo.numerator)) is True, kk
    assert hurwitz(kk, sp.Rational(hi.denominator, hi.numerator)) is False, kk
    for j in range(1, 20):
        bb = 2 + (lo - 2) * Fr(j, 20)
        assert hurwitz(kk, sp.Rational(bb.denominator, bb.numerator)) is True
print('beta_c(k) brackets [bc - 1e-9, bc + 1e-9] confirmed exactly for k = 0..12 (stable below on a grid, unstable above)')
# check minimal polynomials from certificate vanish near the bracket (sign change)
import json
cert = json.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/family_F_certificate.json'))
bsym = sp.symbols('beta')
for row in cert['F5_beta_c_table']:
    p = sp.sympify(row['beta_c_minpoly'], locals={'beta': bsym}); bc = Fr(table[row['k']])
    lo = sp.Rational(str(bc - Fr(1, 10**9))); hi = sp.Rational(str(bc + Fr(1, 10**9)))
    assert p.subs(bsym, lo) * p.subs(bsym, hi) < 0, row['k']
print('certificate minimal polynomials change sign across each bracket: OK')
# k = 0: Geneson's polynomial 797 b^6 - 6468 b^5 + ...
assert cert['F5_beta_c_table'][0]['beta_c_minpoly'].startswith('797*beta**6 - 6468*beta**5')
# plastic: b_p in (23247/10000, 23248/10000)
pl = lambda bb: bb**3 - 3*bb**2 + 2*bb - 1
assert pl(Fr(23247, 10000)) < 0 < pl(Fr(23248, 10000))
for kk in range(2, 9):
    assert hurwitz(kk, sp.Rational(10000, 23248)) is True, kk     # limit-stable just above b_p => beta_c(k) > b_p
for kk in [0, 1] + list(range(9, 13)):
    assert hurwitz(kk, sp.Rational(10000, 23247)) is False, kk    # beta_c(k) < b_p
print('plateau: beta_c(k) > b_p exactly for k = 2..8; beta_c(k) < b_p for k = 0, 1, 9..12: OK')
# Which minor binds?  find, per k, which minor vanishes at beta_c
for kk in table:
    bc = sp.Rational(table[kk])
    vals = {nm: sp.N(m.subs({k: kk, r: 1/bc}), 6) for nm, m in zip(['a0', 'a1', 'D2', 'D3', 'a4'], minors)}
    print(kk, vals)
print('D3(k, r) =', sp.factor(D3)); print('a4 =', sp.factor(a4))
open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/referee_theory_2/minors.txt', 'w').write('D3 = %s\na4 = %s\nD2 = %s\na1 = %s\n' % (D3, a4, D2, a1))
