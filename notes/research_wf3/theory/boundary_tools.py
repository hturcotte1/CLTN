"""Exact tools for alpha = -1/2 (cube-root) boundary graphs: finite-eps-only counterexample mechanism."""
import sympy as sp
from fractions import Fraction as F
import tlib, ctln, rstar
l = sp.symbols('l')
PHI3 = sp.Poly(l**2 + l + 1, l)
def perp_poly(A):
    return sp.Poly(tlib.chi_perp_times_m(A), l)          # n chi_{A|perp}(l), A-picture
def boundary_split(A):
    """returns (mult, quotient_ok): multiplicity of l^2+l+1 in n chi_perp, and whether every OTHER root has Re < -1/2."""
    P = perp_poly(A); mult = 0
    while True:
        q, r = sp.div(P, PHI3)
        if r.is_zero: P = q; mult += 1
        else: break
    co = [F(int(c)) for c in P.all_coeffs()]
    v, _ = tlib.alpha_less_than(None, F(-1, 2), P=[int(c) for c in P.all_coeffs()]) if P.degree() > 0 else (True, '')
    return mult, v
def witness_search(A, R=None, ts=(30, 48, 60, 100, 200, 500, 1000, 3000), xs=(F(21, 20), F(11, 10), F(5, 4), F(3, 2), F(7, 4))):
    n = len(A)
    for t in ts:
        for x in xs:
            beta = 2 + x / t; t_ = F(t)
            eps = 1 / (1 + t_); delta = (beta - 1) / (1 + t_)
            if not ctln.legal(eps, delta): continue
            rec = ctln.analyze_support(A, list(range(n)), eps, delta)
            if rec['stable_fixed_point'] and rec['stable'] is True and not rec['target_free']:
                return {'t': str(t_), 'beta': str(beta), 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)),
                        'legal': True, 'positive': rec['positive'], 'stable': rec['stable'], 'detail': rec['stability_detail'],
                        'clique': rec['clique'], 'x': rec['x']}
    return None
def analyse(A):
    mult, other_ok = boundary_split(A)
    out = {'d6': tlib.to_d6(A), 'n': len(A), 'phi3_mult': mult, 'others_strict': other_ok}
    if mult == 1 and other_ok:
        R, val = rstar.R_at(A)
        out['R'] = None if R is None else str(R)
        out['R_gt_1/8'] = (R is not None and R > sp.Rational(1, 8))
        out['pos_near_2'] = tlib.positive_at(A, F(2) + F(1, 10**4))
    return out
