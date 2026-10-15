"""
family_Q2_proof.py -- two sink classes: q of full-support counterexamples beyond 2 (exact).
Run: python3 family_Q2_proof.py [--instances]   (writes family_Q2_certificate.json)

B' := &F^U|m^QTq? (an n=7 candidate).  W_k := B' + class V1 (size k, receives from all of B' except node 2)
                                          + class V2 (size 3k, receives from all of B' except node 3);
every class node sends to all other nodes; V1, V2 bidirected cliques, mutually bidirected.  n = 7 + 4k.
 (W1) Lemma X1' on explicit graphs k = 1, 2:  n chi_{D'|perp} = z^(4k-1) Phi_k,  Phi_k = z psi_D + 4k chi_D + 1^T adj(zI-D)(k e_2 + 3k e_3).
 (W2) limit polynomial  Phi_k / (4k) -> det(zI - D + bbar 1^T),  bbar = (e_2 + 3 e_3)/4:  = z^3 (2z^2+z+2)^2 / 4,  min Re = -1/4.
 (W3) positivity: det, adj(Mb)1, and the sigma - 1 numerators for S = all-but-2 and all-but-3 have no root in [2, inf) and
      are positive at 3 => (Lemma X2') W_k is positive at EVERY beta > 2, every k.
 (W4) k0(beta*): all Hurwitz minors of -Phi_k(-1/beta* - x) positive for every k >= k0:  5/2 -> 1, 3 -> 3, 7/2 -> 10, 19/5 -> 33.
 (W5) [--instances] exact witnesses n = 19 (k=3, beta=3) and n = 47 (k=10, beta=7/2), t = 10^4.
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import tlib, ctln, multiclass as mc, posfam, classext
A = tlib.from_d6('&F^U|m^QTq?')
CLASSES = [((1, 1, 0, 1, 1, 1, 1), 1), ((1, 1, 1, 0, 1, 1, 1), 3)]
z, k, b = sp.symbols('z k beta')
def main(instances):
    cert = {'family': 'W_k = &F^U|m^QTq? + sink class V1 (size k, in-set all but node 2) + sink class V2 (size 3k, in-set all but node 3); n = 7 + 4k',
            'B_rows': [''.join(map(str, r)) for r in A], 'classes': [[list(S), rho] for S, rho in CLASSES]}
    assert mc.check(A, CLASSES, ks=(1, 2)); cert['W1_formula_checked_k_1_2'] = True
    m = 7; D = sp.Matrix(classext.nonedge_digraph(A))
    bbar = sp.zeros(m, 1); bbar[2] = sp.Rational(1, 4); bbar[3] = sp.Rational(3, 4)
    lim = sp.factor((z * sp.eye(m) - D + bbar * sp.ones(1, m)).det())
    assert sp.expand(lim - z**3 * (2*z**2 + z + 2)**2 / 4) == 0
    cert['W2_limit_polynomial'] = str(lim)
    d, adj = posfam.w_polys(A)
    polys = {'det': d.as_expr()}
    for i, p in enumerate(adj): polys['adj%d' % i] = p.as_expr()
    for S, rho in CLASSES:
        num, _ = posfam.sigma_minus_1(A, S); polys['sigma-1_num_' + ''.join(map(str, S))] = num.as_expr()
    w3 = {}
    for name, p in polys.items():
        Pp = sp.Poly(p, b); nr = Pp.count_roots(2, None) if Pp.degree() > 0 else 0
        assert nr == 0 and Pp.eval(3) > 0
        w3[name] = {'poly': str(sp.factor(p)), 'roots_in_[2,inf)': 0, 'value_at_3': int(Pp.eval(3))}
    cert['W3_positivity_all_beta_gt_2'] = w3
    w4 = {}
    for bstar, k0exp in [(sp.Rational(5, 2), 1), (sp.Integer(3), 3), (sp.Rational(7, 2), 10), (sp.Rational(19, 5), 33)]:
        k0 = mc.k0_for(A, CLASSES, bstar)
        assert k0 == k0exp, (bstar, k0)
        w4[str(bstar)] = {'k0': k0, 'n0': 7 + 4 * k0, 'q_limit': str(bstar - 1)}
    cert['W4_thresholds'] = w4
    if instances:
        wits = []
        for kv, beta in [(3, F(3)), (10, F(7, 2))]:
            G = mc.build(A, CLASSES, kv); n = len(G)
            assert tlib.limit_ok_at(G, beta)
            eps = F(1, 10001); delta = (beta - 1) / 10001
            r = ctln.analyze_support(G, list(range(n)), eps, delta)
            assert ctln.legal(eps, delta) and r['positive'] and r['stable'] is True and not r['clique']
            wits.append({'k': kv, 'n': n, 'd6': tlib.to_d6(G), 'beta': str(beta), 't': '10000', 'eps': str(eps), 'delta': str(delta),
                         'q': str(ctln.q_of(eps, delta)), 'legal': True, 'positive': True, 'stable': True, 'detail': r['stability_detail']})
        cert['W5_witnesses'] = wits
    json.dump(cert, open('family_Q2_certificate.json', 'w'), indent=1)
    print(json.dumps({kk: cert[kk] for kk in ('W2_limit_polynomial', 'W4_thresholds')}, indent=1))
    if instances: print([(w['n'], w['q']) for w in cert['W5_witnesses']])
if __name__ == '__main__':
    main('--instances' in sys.argv)
