"""
family_Y_proof.py -- 8-node core with two sink classes: q -> 3.96 along n = 8 + 5k (exact).
Run: python3 family_Y_proof.py [--instances]   (writes family_Y_certificate.json)

C8 := &G^ZZZr|KjqTo (found by a HEURISTIC core-design search; any 8-node digraph would do for the proof).
Y_k := C8 + class V1 (size 4k, receives from all of C8 except node 0) + class V2 (size k, receives from all except node 2);
class nodes send to every other node; V1, V2 cliques, mutually bidirected.  n = 8 + 5k.
 (Y1) Lemma X1' on the explicit graph k = 1.
 (Y2) limit polynomial det(zI - D + bbar 1^T), bbar = (4 e_0 + e_2)/5 = z^2 (5z^6+5z^5+12z^4+7z^3+7z^2+2z+1)/5 (min Re ~ -0.2016).
 (Y3) positivity at every beta > 2 (det, adj(Mb)1, sigma-1 numerators: no root in [2, inf), signs at 3) => Lemma X2'.
 (Y4) k0(beta*): 3 -> 13, 7/2 -> 24, 4 -> 46, 9/2 -> 117, 24/5 -> 363.
 (Y5) [--instances] exact witness k = 13 (n = 73), beta = 3, t = 10^4 (slow: minutes).
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import tlib, ctln, multiclass as mc, posfam, classext
A = tlib.from_d6('&G^ZZZr|KjqTo')
CLASSES = [((0, 1, 1, 1, 1, 1, 1, 1), 4), ((1, 1, 0, 1, 1, 1, 1, 1), 1)]
z, b = sp.symbols('z beta')
def main(instances):
    cert = {'family': 'Y_k = &G^ZZZr|KjqTo + sink class V1 (size 4k, in-set all but node 0) + sink class V2 (size k, in-set all but node 2); n = 8 + 5k',
            'C8_rows': [''.join(map(str, r)) for r in A], 'classes': [[list(S), rho] for S, rho in CLASSES]}
    assert mc.check(A, CLASSES, ks=(1,)); cert['Y1_formula_checked_k_1'] = True
    m = 8; D = sp.Matrix(classext.nonedge_digraph(A))
    bbar = sp.zeros(m, 1); bbar[0] = sp.Rational(4, 5); bbar[2] = sp.Rational(1, 5)
    lim = sp.factor((z * sp.eye(m) - D + bbar * sp.ones(1, m)).det())
    assert sp.expand(lim - z**2 * (5*z**6 + 5*z**5 + 12*z**4 + 7*z**3 + 7*z**2 + 2*z + 1) / 5) == 0
    cert['Y2_limit_polynomial'] = str(lim)
    d, adj = posfam.w_polys(A)
    polys = {'det': d.as_expr()}
    for i, p in enumerate(adj): polys['adj%d' % i] = p.as_expr()
    for S, rho in CLASSES:
        num, _ = posfam.sigma_minus_1(A, S); polys['sigma-1_num_' + ''.join(map(str, S))] = num.as_expr()
    y3 = {}
    sgn_det = sp.sign(sp.Poly(d.as_expr(), b).eval(3))
    for name, p in polys.items():
        Pp = sp.Poly(p, b); nr = Pp.count_roots(2, None) if Pp.degree() > 0 else 0
        v3 = Pp.eval(3)
        assert nr == 0 and sp.sign(v3) == sgn_det      # every numerator has the sign of det on (2, inf): w > 0, sigma > 1
        y3[name] = {'poly': str(sp.factor(p)), 'roots_in_[2,inf)': 0, 'value_at_3': int(v3)}
    cert['Y3_positivity_all_beta_gt_2'] = y3
    y4 = {}
    for bstar, k0exp in [(sp.Integer(3), 13), (sp.Rational(7, 2), 24), (sp.Integer(4), 46), (sp.Rational(9, 2), 117), (sp.Rational(24, 5), 363)]:
        k0 = mc.k0_for(A, CLASSES, bstar); assert k0 == k0exp, (bstar, k0)
        y4[str(bstar)] = {'k0': k0, 'n0': 8 + 5 * k0, 'q_limit': str(bstar - 1)}
    cert['Y4_thresholds'] = y4
    if instances:
        kv, beta = 13, F(3)
        G = mc.build(A, CLASSES, kv); n = len(G)
        assert tlib.limit_ok_at(G, beta)
        eps = F(1, 10001); delta = (beta - 1) / 10001
        r = ctln.analyze_support(G, list(range(n)), eps, delta)
        assert ctln.legal(eps, delta) and r['positive'] and r['stable'] is True and not r['clique']
        cert['Y5_witness'] = {'k': kv, 'n': n, 'd6': tlib.to_d6(G), 'beta': str(beta), 't': '10000', 'eps': str(eps), 'delta': str(delta),
                              'q': str(ctln.q_of(eps, delta)), 'legal': True, 'positive': True, 'stable': True, 'detail': r['stability_detail']}
    json.dump(cert, open('family_Y_certificate.json', 'w'), indent=1)
    print(json.dumps({kk: cert[kk] for kk in ('Y2_limit_polynomial', 'Y4_thresholds')}, indent=1))
    if instances: print(cert['Y5_witness']['n'], cert['Y5_witness']['q'])
if __name__ == '__main__':
    main('--instances' in sys.argv)
