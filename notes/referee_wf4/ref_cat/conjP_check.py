"""Exact check: catalogue entries with S = all (joins G v K_k) whose D = J - I - A has every eigenvalue Re > -1/2
while u = (D^T + rI)^{-1} 1 > 0 at r = 1/2 and r = 49/100.  Such a D contradicts Conjecture P (notes_theory.md §5.4)."""
import json, sympy as sp
from fractions import Fraction as F
import reprove as R
z, x = sp.symbols('z x')
c = json.load(open(R.CAT + 'families_catalogue_n7.json'))
out = []
for i, e in enumerate(c):
    if not all(e['S']): continue
    A = R.d6_decode(e['base']); m = len(A)
    D = [[0 if a == b else 1 - A[a][b] for b in range(m)] for a in range(m)]
    chiD = sp.Matrix(D).charpoly(z).as_expr()
    co = sp.Poly(sp.expand((-1) ** m * chiD.subs(z, -sp.Rational(1, 2) - x)), x).all_coeffs()
    hur = co[0] > 0 and all(cc > 0 for cc in co) and all(mn > 0 for mn in R.hurwitz_minors(co))
    pos = {}
    for r in (F(1, 2), F(49, 100)):
        uu = R.solve_q([[F(D[b][a]) + (r if a == b else 0) for b in range(m)] for a in range(m)], [F(1)] * m)
        pos[str(r)] = uu is not None and all(v > 0 for v in uu)
    out.append({'i': i, 'base': e['base'], 'chi_D': str(sp.factor(chiD)), 'all_eigs_Re_gt_-1/2': bool(hur), 'u_pos': pos})
viol = [o for o in out if o['all_eigs_Re_gt_-1/2'] and all(o['u_pos'].values())]
print('S=all entries', len(out), ' D with spec Re > -1/2 and u > 0 at r=49/100:', len(viol))
for o in viol[:5]: print(o)
json.dump({'all': out, 'violations': viol}, open('conjP_check.json', 'w'), indent=1)
