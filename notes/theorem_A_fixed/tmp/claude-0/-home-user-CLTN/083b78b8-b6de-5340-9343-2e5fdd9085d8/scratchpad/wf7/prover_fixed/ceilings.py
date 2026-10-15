"""
ceilings.py -- EXACT limit polynomials chi_{Lbar} (Lbar = D - bbar 1^T, bbar = sum rho_j b_j / sum rho_j) of the families used,
the exact algebraic ceiling beta_inf = -1/rho (rho = min real part of a nonzero root of chi_{Lbar}; minimal polynomial via the
resultant Res_w(p(w), p(s - w)), whose roots include 2 Re z for every root z), and exact Hurwitz checks:
 - chi_{Lbar}(-1/beta - x) is Hurwitz (after removing the z^j factor) at the rational beta1 of every region theorem, and
 - NOT Hurwitz at a rational beta just above beta_inf (so by Theorem FP(ii) the family fails there for all large k, every t).
For F: chi_L = z^4 (z^2 + z + 1), rho = -1/2 exactly, so the family fails at every legal (beta > 2) fixed point.
Writes ceilings.json.  Decimal values are for display only.
"""
import json
import sympy as sp
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

z, s, w = sp.symbols('z s w')
G1 = [[0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0]]
FAM = dict(FAMILIES)
FAM['F'] = (None, [(1, 1, 1, 1, 1, 0)], [1])


def lbar_poly(A, classes, rhos):
    m = len(A)
    D = sp.Matrix(fl.nonedge(A))
    bbar = sp.zeros(m, 1)
    for S, rho in zip(classes, rhos):
        for i in range(m):
            if not S[i]: bbar[i] += sp.Rational(rho, sum(rhos))
    return sp.expand((z * sp.eye(m) - D + bbar * sp.ones(1, m)).det(method='berkowitz'))


def hurwitz_at(p, beta):
    """p with zero roots removed; exact Hurwitz test of all roots Re > -1/beta."""
    P = sp.Poly(p, z); deg = P.degree()
    h = sp.Poly(sp.expand((-1) ** deg * p.subs(z, -1 / sp.Rational(beta) - fl.x)), fl.x)
    cs = [F(str(sp.Rational(c))) for c in h.all_coeffs()]
    lc = cs[0]
    cs = [c / lc for c in cs]          # make monic (dividing by a nonzero rational does not move roots)
    return fl.ctln.hurwitz_stable(cs)


out = {}
REGION_BETA1 = {'R': ['43/20', '11/5', '23/10', '12/5', '5/2', '13/5', '27/10', '14/5'], 'Q': ['21/10', '5/2'],
                'Q2': ['5/2', '3', '7/2', '19/5'], 'Y': ['3', '7/2', '4', '9/2', '24/5'], 'F': []}
for nm, (base, classes, rhos) in FAM.items():
    A = G1 if nm == 'F' else fl.from_d6(base)
    p = lbar_poly(A, classes, rhos)
    P = sp.Poly(p, z)
    j = 0
    while P.eval(0) == 0 and P.degree() > 0:
        P = sp.Poly(sp.cancel(P.as_expr() / z), z); j += 1
    q = P.as_expr()
    # exact minimal polynomial of 2*Re(root) candidates
    qs = sp.sqf_part(q)
    res = sp.Poly(sp.resultant(qs.subs(z, w), qs.subs(z, s - w), w), s)
    facs = [f for f, e in sp.factor_list(res.as_expr())[1]]
    rts = sp.Poly(sp.sqf_part(q), z).nroots(n=40, maxsteps=500)   # display/selection only
    rho = min(sp.re(r) for r in rts)
    # identify the factor of res that vanishes at 2*rho (display decimals select the factor; the factor itself is exact)
    best = min(facs, key=lambda f: abs(sp.N(f.subs(s, 2 * rho), 50)))
    beta_inf = -1 / rho
    rec = {'chi_Lbar': str(sp.factor(p)), 'zero_root_multiplicity': j, 'min_Re_display': str(sp.N(rho, 15)),
           'minpoly_of_2rho_factor': str(best), 'beta_inf_display': str(sp.N(beta_inf, 15))}
    checks = []
    for b1 in REGION_BETA1[nm]:
        v = hurwitz_at(q, sp.Rational(b1))
        checks.append({'beta1': b1, 'limit_Hurwitz': v[0]})
        assert v[0] is True, (nm, b1, v)
    # a rational just above the ceiling: limit not Hurwitz
    if nm != 'F':
        above = sp.Rational(int(sp.ceiling(beta_inf * 1000)) + 1, 1000)
        v = hurwitz_at(q, above)
        assert v[0] is False, (nm, above, v)
        checks.append({'beta_above_ceiling': str(above), 'limit_Hurwitz': v[0], 'detail': v[1]})
    else:
        # rho = -1/2 exactly: z^2 + z + 1 divides chi_L
        assert sp.rem(sp.Poly(q, z), sp.Poly(z ** 2 + z + 1, z)).is_zero
        checks.append({'z^2+z+1 divides chi_L': True, 'consequence': 'root with Re = -1/2 < -1/beta for every beta > 2'})
    rec['checks'] = checks
    out[nm] = rec
    print(nm, rec['chi_Lbar'], 'min Re', rec['min_Re_display'], 'beta_inf', rec['beta_inf_display'], flush=True)
json.dump(out, open('ceilings.json', 'w'), indent=1)
