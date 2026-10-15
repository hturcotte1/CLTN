"""
family_Q_proof.py -- a clique-class family whose limit q does NOT tend to 1 (uniform lower bounds on q_n*).
Run: python3 family_Q_proof.py [--instances]     (exact; writes family_Q_certificate.json)

B := &F^S{^ReL|? (one of the n=7 limit candidates);  S = (1,1,1,1,0,1,1);  Q_k := B[S, all; k]  (n = 7 + k):
k bidirected clique nodes, each sending to all of B and receiving from all of B except node 4.
 (Q1) Lemma X1 checked on explicit graphs k = 1..8;  Phi_k = z psi_D + k chi_L  with
      z psi_D = z^3 (7z^4+12z^3+19z^2+11z+4),  chi_L = z^3 (z^4+z^3+3z^2+z+1).
 (Q2) positivity of B and sigma_S for EVERY beta > 2 (no root of det / adj / sigma-1 numerator in (2, inf), signs at 3):
      by Lemma X2, Q_k is positive at every beta > 2, for every k.
 (Q3) for beta* in {5/2, 27/10, 14/5}: all Hurwitz minors of (-1)^7 Phi_k(-1/beta* - x) are polynomials in k with positive
      leading coefficient and no real root >= k0(beta*) = 25, 79, 286 respectively  =>  Q_k limit-stable at beta* for all k >= k0.
 (Q4) as k -> inf the moving spectrum tends to spec(L) (root continuity), whose minimum real part rho is the real part of the
      roots of z^4+z^3+3z^2+z+1 of smaller real part; beta_inf = -1/rho ~ 2.8443 (display only).
 (Q5) [--instances] exact witness k = 25 (n = 32) at beta = 5/2, t = 10^4 via ctln.analyze_support.
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import tlib, classext, posfam, family_prover as fp, ctln
B = tlib.from_d6('&F^S{^ReL|?'); S = (1, 1, 1, 1, 0, 1, 1); T = (1,) * 7
def Q_k(k): return classext.blow(B, S, T, k)
z, k, x, b = sp.symbols('z k x beta')
def main(instances=False):
    cert = {'family': 'Q_k = B[S, all; k], B = &F^S{^ReL|?, S = (1,1,1,1,0,1,1)', 'B_rows': [''.join(map(str, r)) for r in B], 'S': S}
    ok, _ = classext.check_X1(B, S, T, range(1, 9)); assert ok
    Q, P = classext.Phi(B, S, T)
    assert sp.expand(Q.as_expr() - z**3 * (7*z**4 + 12*z**3 + 19*z**2 + 11*z + 4)) == 0
    assert sp.expand(P.as_expr() - z**3 * (z**4 + z**3 + 3*z**2 + z + 1)) == 0
    cert['Q1'] = {'z_psi_D': str(sp.factor(Q.as_expr())), 'chi_L': str(sp.factor(P.as_expr()))}
    d, adj = posfam.w_polys(B); num, _ = posfam.sigma_minus_1(B, S)
    polys = {'det': d.as_expr(), 'sigma_minus_1_num': num.as_expr()}
    for i, p in enumerate(adj): polys['adj%d' % i] = p.as_expr()
    q2 = {}
    for name, p in polys.items():
        Pp = sp.Poly(p, b)
        nroots = Pp.count_roots(2, None) if Pp.degree() > 0 else 0
        # count_roots(2, None) counts roots in [2, inf)
        val3 = Pp.eval(3)
        q2[name] = {'poly': str(sp.factor(p)), 'roots_in_[2,inf)': int(nroots), 'value_at_3': int(val3)}
        assert nroots == 0 and val3 > 0
    cert['Q2_positivity_all_beta_gt_2'] = q2
    q3 = {}
    for bstar, k0 in [(sp.Rational(5, 2), 25), (sp.Rational(27, 10), 79), (sp.Rational(14, 5), 286)]:
        Phi = sp.expand(Q.as_expr() + k * P.as_expr())
        h = sp.Poly(sp.expand(-Phi.subs(z, -1 / bstar - x)), x)
        co = [sp.expand(c) for c in h.all_coeffs()]
        mins = fp.hurwitz_minors(co)
        rec = []
        for mn in [co[0]] + mins:
            Pk = sp.Poly(mn, k)
            assert Pk.LC() > 0
            roots_ge = Pk.count_roots(k0, None) if Pk.degree() > 0 else 0
            assert roots_ge == 0
            rec.append({'degree_in_k': Pk.degree(), 'LC_positive': True, 'real_roots_>=k0': int(roots_ge)})
        q3[str(bstar)] = {'k0': k0, 'n0': 7 + k0, 'minors': rec,
                          'conclusion': 'for every integer k >= %d, Q_k (n = 7 + k >= %d) is limit-stable at beta = %s; with Q2 it is an SFP for all small eps, q -> %s' % (k0, 7 + k0, bstar, bstar - 1)}
        # sanity: explicit graph at k0 exact
        if k0 <= 25: assert tlib.limit_ok_at(Q_k(k0), F(int(bstar.p), int(bstar.q)))   # explicit-graph sanity check (small n only)
    cert['Q3'] = q3
    rts = sp.Poly(z**4 + z**3 + 3*z**2 + z + 1, z).nroots(n=30)
    rho = min(sp.re(r) for r in rts)
    cert['Q4_display_only'] = {'min_Re_spec_L': str(rho), 'beta_inf': str(-1 / rho), 'q_inf': str(-1 / rho - 1)}
    if instances:
        A = Q_k(25); n = len(A); eps = F(1, 10001); delta = F(3, 2) / 10001
        r = ctln.analyze_support(A, list(range(n)), eps, delta)
        assert ctln.legal(eps, delta) and r['positive'] and r['stable'] is True and not r['clique']
        cert['Q5_witness'] = {'k': 25, 'n': n, 'd6': tlib.to_d6(A), 'beta': '5/2', 't': '10000', 'eps': str(eps), 'delta': str(delta),
                              'q': str(ctln.q_of(eps, delta)), 'legal': True, 'positive': True, 'stable': True, 'detail': r['stability_detail'], 'x': r['x']}
    json.dump(cert, open('family_Q_certificate.json', 'w'), indent=1)
    print(json.dumps({kk: cert[kk] for kk in ('Q1', 'Q3', 'Q4_display_only')}, indent=1)[:3000])
if __name__ == '__main__':
    main('--instances' in sys.argv)
