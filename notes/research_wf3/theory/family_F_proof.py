"""
family_F_proof.py -- exact proof data for the infinite LIMIT family F_k (n = 6 + k, k >= 0).
Run: python3 family_F_proof.py      (writes family_F_certificate.json; no floating point in any decision)

G1 = Geneson's graph (paper labels 1..6 -> our 0..5, notes/geneson_2607.21396.md):
    A = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
F_k = G1 + a class V of k nodes, bidirected among themselves; every v in V sends to all six G1 nodes and
      receives from G1 nodes 0..4 (paper vertices 1..5), not from node 5 (paper vertex 6).
Checks:
 (F1) Lemma X1 on explicit graphs k = 1..10.
 (F2) closed forms: z psi_D = z^2 (6z^4+11z^3+16z^2+10z+3),  chi_L = z^4 (z^2+z+1),
      so Phi_k = z^2 R_k,  R_k = (6+k) z^4 + (11+k) z^3 + (16+k) z^2 + 10 z + 3.
 (F3) R_k(-1/2 - x) has positive Hurwitz data for every real k > -1:
      coefficients 16(k+6), 16(k+1), 8(2k+17), 12(k+1), 3k+16;  Delta2 = 64(k+1)(k+16), Delta3 = 8192(k+1)^2.
      Hence alpha(F_k|perp) < -1/2 for every k >= 0.
 (F4) positivity polynomials of G1 (beta A^T + (beta-1) I):  det = -2b^5+7b^4-12b^3+10b^2-6b+1,
      adj1 components, sigma_S - 1 = -(b^2-b+1)/det; Sturm counts show: on [2, 93/40] det < 0 and adj_i < 0
      except adj_4 = b^3-3b^2+2b-1 which has exactly one root b_p in (2, 93/40) (b_p - 1 = plastic number).
      So w_G1 > 0 and sigma_S > 1 exactly on (2, b_p) within [2, 93/40], and adj_4 > 0 > adj_1 on (b_p, 93/40].
 (F5) beta_c(k) data: exact isolating intervals for the smallest root > 2 of the binding Hurwitz minor, k = 0..12.
 (F6) explicit witnesses: tlib.limit_ok_at + ctln.analyze_support for k = 1..10.
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import tlib, ctln, classext, posfam
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ

G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
S1 = (1,1,1,1,1,0); T1 = (1,1,1,1,1,1)
def F_k(k): return classext.blow(G1, S1, T1, k)
z, k, x, b = sp.symbols('z k x beta')

def sturm_count(p, lo, hi):
    return sp.Poly(p, b).count_roots(lo, hi)

def main():
    cert = {'family': 'F_k = Geneson G1 + bidirected clique class V (|V| = k); V -> all of G1; G1 nodes 0..4 -> V; node 5 -/-> V',
            'G1_rows': [''.join(map(str, r)) for r in G1], 'S': S1, 'T': T1}
    assert tlib.to_d6(G1) == open('/home/user/CLTN/n7/geneson.d6').read().strip()
    # (F1)
    ok, bad = classext.check_X1(G1, S1, T1, range(1, 11)); assert ok
    cert['F1_lemmaX1_explicit_k_1_10'] = True
    # (F2)
    Q, P = classext.Phi(G1, S1, T1)
    assert sp.expand(Q.as_expr() - z**2 * (6*z**4 + 11*z**3 + 16*z**2 + 10*z + 3)) == 0
    assert sp.expand(P.as_expr() - z**4 * (z**2 + z + 1)) == 0
    Rk = (6 + k)*z**4 + (11 + k)*z**3 + (16 + k)*z**2 + 10*z + 3
    assert sp.expand(Q.as_expr() + k*P.as_expr() - z**2 * Rk) == 0
    cert['F2'] = {'z_psi_D': str(sp.factor(Q.as_expr())), 'chi_L': str(sp.factor(P.as_expr())), 'R_k': str(Rk)}
    # (F3)
    h = sp.Poly(sp.expand(16 * Rk.subs(z, -sp.Rational(1, 2) - x)), x)
    a = [sp.factor(c) for c in h.all_coeffs()]
    D2 = sp.factor(a[1]*a[2] - a[0]*a[3]); D3 = sp.factor(a[3]*(a[1]*a[2] - a[0]*a[3]) - a[1]**2*a[4])
    assert [sp.expand(t) for t in a] == [sp.expand(t) for t in [16*(k+6), 16*(k+1), 8*(2*k+17), 12*(k+1), 3*k+16]]
    assert sp.expand(D2 - 64*(k+1)*(k+16)) == 0 and sp.expand(D3 - 8192*(k+1)**2) == 0
    cert['F3'] = {'shifted_coeffs': [str(t) for t in a], 'Delta2': str(D2), 'Delta3': str(D3),
                  'conclusion': 'all Hurwitz data positive for k > -1: every root zeta of R_k has Re zeta > -1/2, so alpha(F_k|perp) = -1 - min Re zeta < -1/2'}
    # exact cross-check on explicit graphs
    for kk in range(0, 16):
        v, _ = tlib.alpha_less_than(F_k(kk), F(-1, 2)); assert v is True
    cert['F3_explicit_alpha_lt_-1/2_k_0_15'] = True
    # (F4)
    d, adj = posfam.w_polys(G1)
    detp = d.as_expr(); assert sp.expand(detp - (-2*b**5 + 7*b**4 - 12*b**3 + 10*b**2 - 6*b + 1)) == 0
    num, _ = posfam.sigma_minus_1(G1, S1); assert sp.expand(num.as_expr() + (b**2 - b + 1)) == 0
    lo, hi = sp.Integer(2), sp.Rational(93, 40)
    counts = {'det': sturm_count(detp, lo, hi)}
    for i, p in enumerate(adj): counts['adj%d' % i] = sturm_count(p.as_expr(), lo, hi)
    signs_at_2 = {'det': detp.subs(b, 2)}
    for i, p in enumerate(adj): signs_at_2['adj%d' % i] = p.as_expr().subs(b, 2)
    assert counts['det'] == 0 and signs_at_2['det'] < 0
    for i in range(6):
        if i == 4: assert counts['adj4'] == 1
        else: assert counts['adj%d' % i] == 0 and signs_at_2['adj%d' % i] < 0
    assert signs_at_2['adj4'] < 0 and adj[4].as_expr().subs(b, hi) > 0 and adj[1].as_expr().subs(b, hi) < 0
    assert sp.expand(adj[4].as_expr() - (b**3 - 3*b**2 + 2*b - 1)) == 0
    qq = sp.symbols('q'); assert sp.expand((b**3 - 3*b**2 + 2*b - 1).subs(b, qq + 1) - (qq**3 - qq - 1)) == 0
    cert['F4'] = {'det': str(detp), 'adj1': [str(sp.factor(p.as_expr())) for p in adj], 'sigma_minus_1_numerator': str(num.as_expr()),
                  'root_counts_in_[2,93/40]': {kk: int(v) for kk, v in counts.items()}, 'values_at_beta_2': {kk: int(v) for kk, v in signs_at_2.items()},
                  'conclusion': 'Pos(F_k, beta) for all k >= 1 and beta in (2, b_p), b_p = 1 + plastic number (root of b^3-3b^2+2b-1); fails on (b_p, 93/40]'}
    # (F5) beta_c(k): limit-stable at beta iff R_k(-1/beta - x) Hurwitz; smallest boundary root > 2
    r = sp.symbols('r')
    hb = sp.Poly(sp.expand(Rk.subs(z, -r - x)), x); c = hb.all_coeffs()
    minors = [c[1], c[1]*c[2] - c[0]*c[3], c[3]*(c[1]*c[2] - c[0]*c[3]) - c[1]**2*c[4], c[4]]
    table = []
    for kk in range(0, 13):
        roots = []
        for mnr in minors:
            pb = sp.Poly(sp.numer(sp.together(sp.expand(mnr.subs({k: kk, r: 1/b})))), b)
            for iv, mult in pb.intervals():
                pass
            rr = [t for t in pb.real_roots() if t > 2]
            roots += rr
        bc = min(roots, key=lambda t: sp.N(t, 30))
        iv = sp.Poly(sp.minimal_polynomial(bc, b), b)
        table.append({'k': kk, 'n': 6 + kk, 'beta_c_minpoly': str(iv.as_expr()), 'beta_c_approx_display_only': str(sp.N(bc, 12)),
                      'q_lim_sup_display_only': str(sp.N(min(bc, sp.CRootOf(b**3 - 3*b**2 + 2*b - 1, 0)) - 1, 12))})
    cert['F5_beta_c_table'] = table
    # (F6) witnesses
    wits = []
    for kk, beta in [(1, F(21, 10)), (2, F(58, 25)), (3, F(58, 25)), (4, F(58, 25)), (5, F(58, 25)), (6, F(21, 10)), (8, F(21, 10)), (10, F(21, 10))]:
        A = F_k(kk); n = len(A)
        assert tlib.limit_ok_at(A, beta)
        t, eps, delta, rec = tlib.sfp_witness(A, beta, ts=(20, 50, 100, 1000))
        wits.append({'k': kk, 'n': n, 'd6': tlib.to_d6(A), 'beta': str(beta), 't': str(t), 'eps': str(eps), 'delta': str(delta),
                     'q': str(ctln.q_of(eps, delta)), 'legal': ctln.legal(eps, delta), 'positive': rec['positive'], 'stable': rec['stable'],
                     'detail': rec['stability_detail'], 'clique': rec['clique']})
    cert['F6_witnesses'] = wits
    json.dump(cert, open('family_F_certificate.json', 'w'), indent=1)
    print(json.dumps({kk: cert[kk] for kk in ('F2', 'F3')}, indent=1))
    for row in table: print(row)
    for w in wits: print(w['n'], w['d6'], w['beta'], w['t'], w['q'], w['detail'])

if __name__ == '__main__':
    main()
