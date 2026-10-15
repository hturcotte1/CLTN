"""
theorem_A.py -- EXACT proof that ONE legal parameter point carries an infinite family of non-clique stable full supports.

Family R_k := R[S, all; k],  R = &F^U|m^QTq? (one of the 124 signed-off n=7 counterexamples),  S = all nodes but node 2:
k new nodes forming a bidirected clique; each new node sends to every node of R; node i of R sends to the new nodes iff i != 2.
n = 7 + k.  R_0 = R.

Fixed point: (eps, delta) = (1/21, 11/210)  <=>  t = 20, beta = 21/10;  q = 22/21 > 1 (legal).  This is the same point at which
Geneson's clique-class family F_k is stable only for k <= 40 (see neg_F.py).

Proof steps (each asserted exactly, no floats):
 (A1) Lemma P identity det(zI - aJ_n - D') = z^(k-1) Psi_k(z) on the explicit graphs k = 1..4 at a in {3/7, 5, 200/21}.
 (A2) Psi_k(z) = z det(zI - D - aJ) - a k chi_L(z), chi_L(z) = det(zI - D + b 1^T) = z^3 (z^4 + z^3 + 3z^2 + z + 1).
 (A3) at a = t/beta = 200/21, r = 1/beta = 10/21: h_k(x) = (-1)^8 Psi_k(-r - x) is monic of degree 8 in x; its 8 Hurwitz
      minors are polynomials in k; each is > 0 for every real k >= 0 (exact: positive leading coefficient, positive at 0,
      no real root in [0, inf) by Sturm count; in fact every coefficient is positive).
 (A4) positivity: w_R(21/10) > 0 entrywise and sigma_S(21/10) = sum_{i in S} w_i > 1 (exact rationals)  =>  (Lemma X2) every R_k
      has x > 0 at (eps, delta).
 (A5) R is not a clique (so no R_k is).
 Conclusion: for EVERY integer k >= 0, the full support of R_k (n = 7 + k nodes) is a stable positive fixed point at
 (eps, delta) = (1/21, 11/210).  Writes theorem_A_certificate.json.
"""
import json, time, sys
import sympy as sp
from fractions import Fraction as F
import fixlib as fl

BASE = '&F^U|m^QTq?'
S = (1, 1, 0, 1, 1, 1, 1)
BETA = F(21, 10); T = F(20)


def main(base=BASE, S=S, beta=BETA, t=T, kmin=0, out='theorem_A_certificate.json', check_ks=(1, 2, 3, 4)):
    t0 = time.time()
    A = fl.from_d6(base); m = len(A)
    cert = {'family': '%s[S, all; k], S = %s, n = %d + k' % (base, list(S), m), 'base_rows': [''.join(map(str, r)) for r in A],
            'S': list(S), 'beta': str(beta), 't': str(t)}
    eps, delta = fl.eps_delta(beta, t)
    assert fl.ctln.legal(eps, delta)
    cert.update(eps=str(eps), delta=str(delta), q=str(fl.ctln.q_of(eps, delta)), legal=True)
    av = sp.Rational(str(t / beta)); rv = sp.Rational(str(1 / beta))
    # (A1)
    for kk in check_ks:
        assert fl.check_identity(A, [S], [kk], avals=(sp.Rational(3, 7), sp.Integer(5), av)), kk
    cert['A1_identity_checked_on_explicit_graphs_k'] = list(check_ks)
    # (A2)
    P0, Ps = fl.psi_parts(A, [S])
    chiL = Ps[0]
    cert['A2'] = {'z_det(zI-D-aJ)': str(sp.factor(P0)), 'chi_L': str(sp.factor(chiL))}
    Psi = sp.expand(P0 - fl.a_s * fl.k * chiL)
    # (A3)
    co = fl.shifted_coeffs(Psi.subs(fl.a_s, av), rv)
    mins = fl.hurwitz_minors(co)
    rec = []
    for j, mn in enumerate(mins):
        ok, info = fl.positive_for_all_k_ge(mn, kmin)
        assert ok, (j, info)
        Pk = sp.Poly(mn, fl.k)
        allpos = all(c > 0 for c in Pk.all_coeffs())
        rec.append({'minor': j + 1, 'degree_in_k': Pk.degree(), 'positive_on_[%d,inf)' % kmin: True, 'info': info,
                    'all_coefficients_positive': allpos, 'poly': str(Pk.as_expr())})
    cert['A3_hurwitz_minors_in_k'] = rec
    cert['A3_coefficients_of_h_k'] = [str(c) for c in co]
    # (A4)
    pos, w, sig = fl.positivity_core(A, [S], beta)
    assert pos
    cert['A4_positivity'] = {'w_G': [str(v) for v in w], 'sigma_S': [str(s) for s in sig], 'ok': True}
    # (A5)
    assert not all(A[i][j] for i in range(m) for j in range(m) if i != j)
    cert['A5_core_not_clique'] = True
    # independent cross-check against ctln on explicit graphs (small k): Psi-based verdict == ctln verdict == stable
    xc = []
    for kk in range(kmin, kmin + 6):
        G = fl.build(A, [S], [kk]); n = len(G)
        r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
        v = fl.stable_at_k(A, [S], [1], beta, t, kk)
        assert r['positive'] and r['stable'] is True and v[0] is True and not r['clique']
        xc.append({'k': kk, 'n': n, 'ctln_stable': True, 'ctln_positive': True, 'psi_hurwitz': True})
    cert['crosscheck_ctln_vs_Psi'] = xc
    cert['conclusion'] = ('for every integer k >= %d the full support of %s[S,all;k] (n = %d + k) is a stable positive fixed point '
                          'at (eps, delta) = (%s, %s), q = %s; none of these graphs is a clique.  CANDIDATE counterexamples pending '
                          'hand-verification (CLAUDE.md section 3).' % (kmin, base, m, eps, delta, cert['q']))
    cert['seconds'] = round(time.time() - t0, 1)
    json.dump(cert, open(out, 'w'), indent=1)
    print(json.dumps({kk: cert[kk] for kk in ('A2', 'conclusion')}, indent=1))
    print('minor degrees', [x['degree_in_k'] for x in rec], 'all coeffs positive:', [x['all_coefficients_positive'] for x in rec])
    return cert


if __name__ == '__main__':
    main()
