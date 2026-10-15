"""
theorem_B.py -- EXACT region theorems: a sink-class family is an SFP at EVERY legal (eps, delta) with beta <= beta1, for
every k >= k0 (one k0 for the whole region, i.e. uniformly in t).

For a family (core A, in-sets S_j, proportions rho_j; class sizes rho_j k) and a rational beta1 > 2:
 (B1) positivity on [2, beta1]: det(Mb), every entry of adj(Mb) 1 (Mb = beta A^T + (beta-1) I) and every numerator of
      sigma_{S_j} - 1 have no real root in [2, beta1] and the sign pattern 'w > 0, sigma_j > 1' holds at one point
      (exact Sturm counts)  =>  w_G(beta) > 0 and sigma_{S_j}(beta) > 1 on [2, beta1]  =>  (Lemma X2') every member positive.
 (B2) r1 = 1/beta1, a1 = r1^2/(1 - 2 r1)  (legal with beta <= beta1  =>  r >= r1 and a = t/beta > a1).
      g = Psi / z^j0 (Lemma P), h(x) = (-1)^deg g(-r1 - x); every Hurwitz minor of h, a polynomial in (a, k), is shown > 0 on
      {a >= a1, k >= k0} by the sign-pattern certificate: after a = a1 + u, k = k0 + v all coefficients are >= 0 and the
      constant term is > 0.
 (B3) monotonicity in r: spec(aJ + D') does not depend on r, so 'every eigenvalue has Re > -r1' gives Re > -r for r >= r1.
 Conclusion: for every legal (eps, delta) with 1 + delta/eps <= beta1 and every integer k >= k0, the full support of G_k is
 a stable positive fixed point.
Writes theorem_B_certificate.json.  Floats are used ONLY in heur_k0 to guess where to start k0; every decision is exact.
"""
import json, sys, time
import numpy as np
import sympy as sp
from fractions import Fraction as F
import fixlib as fl

bsym = sp.symbols('beta')

FAMILIES = {
    'R': ('&F^U|m^QTq?', [(1, 1, 0, 1, 1, 1, 1)], [1]),
    'Q': ('&F^S{^ReL|?', [(1, 1, 1, 1, 0, 1, 1)], [1]),
    'Q2': ('&F^U|m^QTq?', [(1, 1, 0, 1, 1, 1, 1), (1, 1, 1, 0, 1, 1, 1)], [1, 3]),
    'Y': ('&G^ZZZr|KjqTo', [(0, 1, 1, 1, 1, 1, 1, 1), (1, 1, 0, 1, 1, 1, 1, 1)], [4, 1]),
}


def positivity_interval(A, classes, beta1):
    m = len(A)
    Mb = sp.Matrix(m, m, lambda i, j: bsym * A[j][i] + ((bsym - 1) if i == j else 0))
    d = sp.expand(Mb.det(method='berkowitz'))
    adjv = Mb.adjugate() * sp.ones(m, 1)
    adj = [sp.expand(adjv[i]) for i in range(m)]
    nums = [sp.expand(bsym * sum(adj[i] for i in range(m) if S[i]) - d) for S in classes]
    polys = [('det', d)] + [('adj%d' % i, p) for i, p in enumerate(adj)] + [('sigma%d-1_num' % j, p) for j, p in enumerate(nums)]
    rec = {}
    bm = (2 + sp.Rational(beta1)) / 2
    sd = sp.sign(d.subs(bsym, bm))
    for name, p in polys:
        P = sp.Poly(p, bsym)
        nr = P.count_roots(2, sp.Rational(beta1)) if P.degree() > 0 else 0
        val = P.eval(bm)
        ok = (nr == 0) and (sp.sign(val) == sd) and val != 0
        rec[name] = {'roots_in_[2,beta1]': int(nr), 'sign_matches_det_at_mid': bool(ok)}
        if not ok:
            return False, rec
    # w_i = beta adj_i / det > 0 and sigma_j - 1 = num_j / det > 0 on [2, beta1] since no sign change and same sign as det
    return True, rec


def heur_k0(g, r1, a1, kmax=5000):
    """float steer only: smallest k (on a grid) beyond which min Re roots(g) > -r1 for sampled a >= a1."""
    P = sp.Poly(g, fl.z)
    cs = [sp.lambdify((fl.a_s, fl.k), c, 'numpy') for c in P.all_coeffs()]
    r1 = float(r1); a1 = float(a1)
    def ok(a, kk):
        return min(np.roots([float(c(a, kk)) for c in cs]).real) > -r1
    asamp = a1 * np.array([1.0001, 1.01, 1.1, 1.5, 2, 4, 10, 100, 1e4, 1e6])
    ks = list(range(0, 200)) + list(range(200, kmax, 10))
    bad = [kk for kk in ks if not all(ok(a, kk) for a in asamp)]
    return (max(bad) + 1) if bad else 0


def region(name, beta1, k0_hint=None, max_tries=12):
    base, classes, rhos = FAMILIES[name]
    A = fl.from_d6(base)
    beta1 = sp.Rational(beta1); r1 = 1 / beta1; a1 = fl.legal_a_bound(r1)
    t0 = time.time()
    pos_ok, pos_rec = positivity_interval(A, classes, beta1)
    g, j0 = fl.reduced_psi(A, classes, rhos)
    co, mins = fl.region_minors(g, r1)
    k0 = k0_hint if k0_hint is not None else heur_k0(g, r1, a1)
    tries = []
    for _ in range(max_tries):
        res = [fl.sign_pattern_certificate(mn, a1, k0) for mn in mins]
        tries.append({'k0': int(k0), 'per_minor_ok': [bool(x[0]) for x in res], 'neg_coeffs': [int(x[1]) for x in res]})
        if all(x[0] for x in res):
            break
        k0 = max(k0 + 1, int(k0 * 1.25) + 1)
    proved = all(tries[-1]['per_minor_ok']) and pos_ok
    out = {'family': name, 'base': base, 'classes': [list(S) for S in classes], 'rhos': rhos,
           'beta1': str(beta1), 'r1': str(r1), 'a1': str(a1), 't1_legal_boundary': str(a1 * beta1),
           'positivity_on_[2,beta1]': pos_ok, 'positivity_detail': pos_rec, 'z_power_removed': j0,
           'g': str(g), 'tries': tries, 'k0': int(tries[-1]['k0']) if proved else None, 'PROVED': proved,
           'statement': ('for every legal (eps, delta) with 1 + delta/eps <= %s and every integer k >= %s, the full support of '
                         '%s_k is a stable positive fixed point' % (beta1, tries[-1]['k0'], name)) if proved else 'not proved',
           'seconds': round(time.time() - t0, 1)}
    print(name, beta1, 'pos', pos_ok, 'k0', out['k0'], 'tries', [(x['k0'], sum(x['neg_coeffs'])) for x in tries], out['seconds'], 's', flush=True)
    return out


JOBS = [('R', '43/20', 0), ('R', '11/5', None), ('R', '23/10', None), ('R', '12/5', None), ('R', '5/2', None),
        ('R', '13/5', None), ('R', '27/10', None), ('R', '14/5', None),
        ('Q', '21/10', None), ('Q', '5/2', None),
        ('Q2', '5/2', None), ('Q2', '3', None), ('Q2', '7/2', None), ('Q2', '19/5', None),
        ('Y', '3', None), ('Y', '7/2', None), ('Y', '4', None), ('Y', '9/2', None), ('Y', '24/5', None)]

if __name__ == '__main__':
    sel = sys.argv[1:]
    results = []
    for name, b1, hint in JOBS:
        if sel and name not in sel: continue
        results.append(region(name, b1, hint))
        json.dump(results, open('theorem_B_certificate%s.json' % ('_' + '_'.join(sel) if sel else ''), 'w'), indent=1)
