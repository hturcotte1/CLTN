"""
selfjoin.py -- EXACT: a second, structurally different infinite family at ONE fixed point: the k-fold self-join
R^{vk} := R v R v ... v R  (k disjoint copies of R = &F^U|m^QTq?, every pair of nodes in different copies joined by edges in
BOTH directions).  n = 7k.  This is NOT a blow-up of a fixed core (Lemmas U, U'' do not apply) and contains no clique class
of growing size that sends to everything.

Lemma SJ (proved in notes_fixed.md section 12).  For any digraph G (m nodes, non-edge matrix D) and k >= 1:
  (a) spec(D' + aJ_{km}) = spec(D + akJ_m)  u  spec(D)^(k-1)        (D' = I_k (x) D, J_{km} = J_k (x) J_m)
  (b) positivity of G^{vk} at beta  <=>  w_G(beta) = (A^T + (1-1/beta) I)^{-1} 1 > 0   (k-free, t-free)
So with P(z, c) := det(zI - D - cJ) = chi_D(z) - c psi_D(z):  if every root of P(., c) has Re > -r1 for every real c >= 0,
then (Lemma S + monotonicity in r) G^{vk} is a stable full support at every legal point with beta <= 1/r1, for every k >= 1,
as soon as w_G(beta) > 0 there.

This script (no floats in any decision):
 (J1) P(z, c) for R; removes the common z^j0 factor; h(x) = (-1)^deg P(-r1 - x, c) at r1 = 1/beta1; every Hurwitz minor
      (a polynomial in c) is > 0 for every real c >= 0 (positive at 0, positive leading coefficient, no root in [0, inf) by
      exact Sturm count; also reports whether all coefficients are positive).  beta1 in BETAS.
 (J2) exact positivity of w_R on [2, beta1] (Sturm counts on det and adj(Mb) 1, as in theorem_B.positivity_interval with no
      sigma condition).
 (J3) identity (a) checked on explicit R^{vk}, k = 2, 3, by exact evaluation of det(zI - aJ - D') at n + 2 rational z.
 (J4) ctln.analyze_support instances at (eps, delta) = (1/21, 11/210) for k in KS (argv) -> cert_selfjoin_R_eps1_21.json.
Writes selfjoin_certificate.json.
"""
import json, sys, time, hashlib
import sympy as sp
from fractions import Fraction as F
import fixlib as fl
from theorem_B import positivity_interval

z, c, x = sp.symbols('z c x')
BASE = '&F^U|m^QTq?'
BETAS = ['21/10', '43/20', '11/5']


def selfjoin(A, kk):
    m = len(A); n = m * kk
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            G[i][j] = A[i % m][j % m] if (i // m) == (j // m) else 1
    return G


def locus_minors(A, beta1):
    m = len(A)
    D = sp.Matrix(fl.nonedge(A))
    P = sp.expand((z * sp.eye(m) - D - c * sp.ones(m, m)).det(method='berkowitz'))
    Pz = sp.Poly(P, z)
    co = Pz.all_coeffs()[::-1]
    j0 = 0
    while sp.expand(co[j0]) == 0: j0 += 1
    g = sp.expand(sum(cc * z ** (i - j0) for i, cc in enumerate(co) if i >= j0))
    deg = sp.Poly(g, z).degree()
    r1 = 1 / sp.Rational(beta1)
    h = sp.Poly(sp.expand((-1) ** deg * g.subs(z, -r1 - x)), x)
    hc = [sp.expand(cc) for cc in h.all_coeffs()]
    assert hc[0] == 1
    mins = fl.hurwitz_minors(hc)
    rec = []
    ok_all = True
    for j, mn in enumerate(mins):
        Pc = sp.Poly(mn, c)
        if Pc.degree() <= 0:
            ok = Pc.LC() > 0; info = 'constant'
        else:
            ok = (Pc.LC() > 0) and (Pc.eval(0) > 0) and Pc.count_roots(0, None) == 0
            info = 'roots in [0,inf): %d' % Pc.count_roots(0, None)
        ok_all &= bool(ok)
        rec.append({'minor': j + 1, 'degree_in_c': Pc.degree(), 'positive_on_[0,inf)': bool(ok), 'info': info,
                    'all_coefficients_positive': all(v > 0 for v in Pc.all_coeffs()), 'poly': str(Pc.as_expr())})
    return {'P(z,c)': str(sp.factor(P)), 'z_power_removed': j0, 'g': str(g), 'beta1': str(beta1), 'r1': str(r1),
            'h_coeffs': [str(v) for v in hc], 'minors': rec, 'ALL_MINORS_POSITIVE_FOR_ALL_c>=0': ok_all}


def check_decomposition(A, kk, avals=(F(3, 7), F(200, 21))):
    """exact: det(zI - aJ_n - D') == det(zI - D - akJ_m) * chi_D(z)^(k-1) at n + 2 rational z (both monic deg n)."""
    G = selfjoin(A, kk); n = len(G); m = len(A)
    Dp = fl.nonedge(G); D = fl.nonedge(A)
    for av in avals:
        for zi in range(-(n // 2) - 1, n - n // 2 + 1):
            zz = F(zi) + F(1, 3)
            L = fl.ctln.det_exact([[(zz if i == j else F(0)) - av - Dp[i][j] for j in range(n)] for i in range(n)])
            R1 = fl.ctln.det_exact([[(zz if i == j else F(0)) - av * kk - D[i][j] for j in range(m)] for i in range(m)])
            R2 = fl.ctln.det_exact([[(zz if i == j else F(0)) - D[i][j] for j in range(m)] for i in range(m)])
            if L != R1 * R2 ** (kk - 1):
                return False
    return True


def main(ks):
    t0 = time.time()
    A = fl.from_d6(BASE)
    cert = {'family': 'k-fold self-join of %s (n = 7k)' % BASE, 'base_rows': [''.join(map(str, r)) for r in A]}
    # (J1), (J2)
    cert['J1_J2'] = []
    for b1 in BETAS:
        loc = locus_minors(A, b1)
        pos_ok, pos_rec = positivity_interval(A, [], sp.Rational(b1))
        loc['positivity_w_R_on_[2,beta1]'] = pos_ok
        loc['PROVED'] = bool(loc['ALL_MINORS_POSITIVE_FOR_ALL_c>=0'] and pos_ok)
        cert['J1_J2'].append(loc)
        print('beta1', b1, 'minors>0 for all c>=0:', loc['ALL_MINORS_POSITIVE_FOR_ALL_c>=0'], 'pos', pos_ok,
              [(mm['minor'], mm['info'], mm['all_coefficients_positive']) for mm in loc['minors']], flush=True)
    # (J3)
    cert['J3_decomposition_checked_k'] = []
    for kk in (2, 3):
        assert check_decomposition(A, kk), kk
        cert['J3_decomposition_checked_k'].append(kk)
    print('J3 decomposition identity verified on explicit graphs k = 2, 3', flush=True)
    json.dump(cert, open('selfjoin_certificate.json', 'w'), indent=1)
    # (J4) instance certificates
    eps, delta = F(1, 21), F(11, 210)
    assert fl.ctln.legal(eps, delta)
    out = {'job': 'selfjoin_R_eps1_21', 'family': 'self-join R^{vk}', 'base': BASE, 'eps': str(eps), 'delta': str(delta),
           'q': str(fl.ctln.q_of(eps, delta)), 'legal': True, 'members': []}
    for kk in ks:
        t1 = time.time()
        G = selfjoin(A, kk); n = len(G)
        r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
        cp = r.pop('charpoly')
        h = hashlib.sha256(json.dumps(cp).encode()).hexdigest()
        ok = r['positive'] and r['stable'] is True and not r['clique'] and r['stable_fixed_point']
        rec = {'k': kk, 'n': n, 'd6': fl.to_d6(G), 'positive': r['positive'], 'stable': r['stable'],
               'stability_detail': r['stability_detail'], 'clique': r['clique'], 'stable_fixed_point': r['stable_fixed_point'],
               'charpoly_sha256': h, 'x': r['x'], 'CERTIFIED': bool(ok), 'seconds': round(time.time() - t1, 1)}
        out['members'].append(rec)
        print('selfjoin k', kk, 'n', n, 'certified', ok, r['stability_detail'], rec['seconds'], 's', flush=True)
        json.dump(out, open('cert_selfjoin_R_eps1_21.json', 'w'), indent=1)
    cert['seconds'] = round(time.time() - t0, 1)
    json.dump(cert, open('selfjoin_certificate.json', 'w'), indent=1)


if __name__ == '__main__':
    ks = [int(v) for v in sys.argv[1:]] or [1, 2, 3, 5]
    main(ks)
