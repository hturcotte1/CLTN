"""
conjP_refutation.py -- EXACT refutation of 'Conjecture P' (notes/research_wf3/theory/notes_theory.md section 5.4):
  'If D != 0 is a 0/1 matrix with zero diagonal and (D^T + rI)^{-1} 1 > 0 for some r in (0, 1/2), then D has an eigenvalue
   with Re <= -1/2.'
Counterexample: D = J - I - A for A = &F^U|m^QTq? (one of the signed-off n=7 counterexamples), r = 10/21:
  (P1) chi_D(z) = z^2 (z^5 - 3z^2 - 2z - 1); (-1)^7 chi_D(-1/2 - x) and (-1)^7 chi_D(-10/21 - x) are Hurwitz (exact minors):
       every eigenvalue of D has Re > -10/21 > -1/2.
  (P2) u = (D^T + (10/21) I)^{-1} 1 > 0 entrywise (exact rationals); also U = 1^T u > 1.
Consequence noted in the notes as 'if P is true, joins G v K_k fail for large k' -- the opposite holds here:
  (P3) the join family R v K_k (S = all nodes) is a stable positive full support at (eps, delta) = (1/21, 11/210) for EVERY k >= 0
       (Lemma P with b = 0: Hurwitz minors in k positive on [0, inf), exact), and the members k = 1, 5, 20 are certified by
       ctln.analyze_support.
Writes conjP_refutation.json.
"""
import json
import sympy as sp
from fractions import Fraction as F
import fixlib as fl

A = fl.from_d6('&F^U|m^QTq?'); m = len(A)
D = sp.Matrix(fl.nonedge(A))
z, x = fl.z, fl.x
chi = sp.expand((z * sp.eye(m) - D).det(method='berkowitz'))
assert sp.expand(chi - z ** 2 * (z ** 5 - 3 * z ** 2 - 2 * z - 1)) == 0
out = {'A_rows': [''.join(map(str, r)) for r in A], 'chi_D': str(sp.factor(chi))}
for name, rr in (('-1/2', sp.Rational(1, 2)), ('-10/21', sp.Rational(10, 21))):
    h = sp.Poly(sp.expand((-1) ** m * chi.subs(z, -rr - x)), x)
    v = fl.ctln.hurwitz_stable([F(str(sp.Rational(c))) for c in h.all_coeffs()])
    out['P1_all_eigs_Re_>_' + name] = v
    assert v[0] is True
r = sp.Rational(10, 21)
u = (D.T + r * sp.eye(m)).LUsolve(sp.ones(m, 1))
assert all(val > 0 for val in u) and sum(u) > 1
out['P2_u'] = [str(val) for val in u]; out['P2_U'] = str(sum(u))
# (P3) join family at the fixed point
S = (1,) * m
beta, t = F(21, 10), F(20)
k0, mins = fl.k_threshold(A, [S], [1], beta, t, kmin=0)
assert k0 == 0
out['P3_join_k0_at_(1/21,11/210)'] = k0
pos, w, sig = fl.positivity_core(A, [S], beta)
assert pos
inst = []
for kk in (1, 5, 20):
    rec = fl.certify_member(A, [S], [kk], beta, t)
    ok = rec['legal'] and rec['positive'] and rec['stable'] is True and not rec['clique']
    assert ok
    inst.append({'k': kk, 'n': rec['n'], 'd6': rec['d6'], 'legal': True, 'positive': True, 'stable': True,
                 'detail': rec['stability_detail'], 'clique': False})
out['P3_instances'] = inst
out['conclusion'] = ('Conjecture P is FALSE (exact).  The join R v K_k is an SFP at (1/21, 11/210) for every k >= 0 '
                     '(members are CANDIDATE counterexamples pending hand-verification).')
json.dump(out, open('conjP_refutation.json', 'w'), indent=1)
print(json.dumps(out, indent=1)[:2500])
