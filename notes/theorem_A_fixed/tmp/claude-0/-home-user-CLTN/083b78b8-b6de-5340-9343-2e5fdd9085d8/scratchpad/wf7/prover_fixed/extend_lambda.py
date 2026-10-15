"""
extend_lambda.py -- EXACT check that the arms-forest cores of the parallel Theorem U (wf8/prover_q, notes_prover_q.md section 3)
give fixed-parameter infinite families BEYOND beta_Y = 4.9607, verified here with this directory's own tools (ctln.charpoly,
ctln.solve_exact, ctln.hurwitz_stable, Lemma P / Lemma S / Theorem FP of notes_fixed.md).  The core builder is re-implemented
below from the description in notes_prover_q.md section 3.1 (c0 -> c1 -> c2 -> c0 D-cycle; arms c2 -> a -> z_a -> target).

Core c = [1, 1, 11, 7, 10, 3, 1] (the 65-node core of notes_prover_q.md section 5.2), single sink class S = all nodes but c2
(b = e_c2), G_k := core[S, all; k], n = m + k.
 (E1) Lemma T: chi_L(z) = det(zI - D + e_c2 1^T) = z^(m-6) p(z), p = z^6 + z^5 + 11 z^4 + 7 z^3 + 10 z^2 + 3 z + 1
      (ctln.charpoly on the explicit integer matrix; also ctln._charpoly_direct).
 (E2) for beta in BETAS: (-1)^6 p(-1/beta - x) Hurwitz (ctln.hurwitz_stable, exact)  <=>  every root of p has Re > -1/beta.
 (E3) positivity at each beta: w_core = (A^T + (1 - 1/beta) I)^{-1} 1 > 0 and sigma_S > 1 (ctln.solve_exact).
 => Theorem FP(i): at EVERY legal (eps, delta) with that beta, G_k is a stable positive non-clique full support for all large k.
 (E4) exact pointwise threshold at (beta, t) = POINT: Psi_k = z chi_{D+aJ} - a k chi_L (Lemma P) = z^(m-6) g_k(z), g_k of
      degree 7, affine in k; Hurwitz minors of (-1)^7 g_k(-r - x) interpolated exactly in k (degree <= j); k0 = smallest
      integer with every minor > 0 on [k0, inf) (exact real-root isolation, sympy as an assistant on integer polynomials).
 (E5) ctln.analyze_support on the explicit graph G_k for k in CTLN_KS at POINT: verdict compared with the Psi verdict.
Writes extend_lambda.json.
"""
import sys, json, time
from fractions import Fraction as F
import sympy as sp
import fixlib as fl
ctln = fl.ctln

C = [1, 1, 11, 7, 10, 3, 1]
BETAS = [F(5), F(11, 2), F(28, 5), F(6)]
POINT = (F(5), F(1))            # beta = 5, t = 1  <=>  (eps, delta) = (1/2, 2)
CTLN_KS = [int(v) for v in sys.argv[1:]]


def build_core(c):
    d = len(c) - 1
    names = ['c0', 'c1', 'c2']; arcs = [('c0', 'c1'), ('c1', 'c2'), ('c2', 'c0')]
    first = {}
    spec = [(2, 'c2', 1, c[2] - 1)]
    for j in range(3, d + 1):
        spec.append((j, 'c1', (j - 1) // 2, c[j]) if j % 2 else (j, 'c0', j // 2 - 1, c[j]))
    spec.sort(key=lambda s: s[2])
    for j, ex, dep, cnt in spec:
        for i in range(cnt):
            a, zz = 'a%d_%d' % (j, i), 'z%d_%d' % (j, i)
            tgt = ex if dep == 1 else first[(ex, dep - 1)]
            names += [a, zz]; arcs += [('c2', a), (a, zz), (zz, tgt)]
            first.setdefault((ex, dep), a)
    idx = {v: i for i, v in enumerate(names)}; m = len(names)
    D = [[0] * m for _ in range(m)]
    for u, v in arcs: D[idx[u]][idx[v]] = 1
    return D, names


def padd(p, q):
    n = max(len(p), len(q)); p = [F(0)] * (n - len(p)) + p; q = [F(0)] * (n - len(q)) + q
    return [x + y for x, y in zip(p, q)]


def pmul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): out[i + j] += x * y
    return out


def shift(p, r):
    out = [F(0)]
    for c in p: out = padd(pmul(out, [F(-1), -r]), [c])
    while len(out) > 1 and out[0] == 0: out = out[1:]
    return out


def minors(co):
    d = len(co) - 1
    aa = lambda i: co[i] if 0 <= i <= d else F(0)
    H = [[aa(2 * (i + 1) - (j + 1)) for j in range(d)] for i in range(d)]
    return [ctln.det_exact([row[:j] for row in H[:j]]) for j in range(1, d + 1)]


def main():
    t0 = time.time(); out = {'core_c': C}
    D, names = build_core(C); m = len(D); d = len(C) - 1
    A = [[0 if i == j else 1 - D[i][j] for j in range(m)] for i in range(m)]
    S = tuple(0 if nm == 'c2' else 1 for nm in names)
    ic2 = names.index('c2')
    out['m'] = m; out['core_d6'] = fl.to_d6(A)
    assert any(D[i][j] for i in range(m) for j in range(m))          # core is not a clique
    # (E1)
    L = [[F(D[i][j] - (1 if i == ic2 else 0)) for j in range(m)] for i in range(m)]
    chiL = ctln.charpoly(L)
    assert chiL == [F(x) for x in C] + [F(0)] * (m - d), 'Lemma T fails'
    out['E1_lemma_T'] = 'chi_L = z^%d * p, p = %s' % (m - d, C)
    print('E1 ok, m =', m, round(time.time() - t0, 1), 's', flush=True)
    # (E2), (E3)
    out['E2_E3'] = []
    p = [F(x) for x in C]
    for beta in BETAS:
        r = 1 / beta
        h = shift(p, r)                     # p(-r - x); degree 6, leading (-1)^6 = 1
        hv = ctln.hurwitz_stable(h)
        pos, w, sig = fl.positivity_core(A, [S], beta)
        rec = {'beta': str(beta), 'roots_of_p_Re_>_-1/beta': hv[0], 'detail': hv[1], 'positivity': pos,
               'sigma_S': str(sig[0]) if sig else None, 'min_w': str(min(w)) if w else None,
               'FP_i_applies': bool(hv[0] is True and pos)}
        out['E2_E3'].append(rec)
        print('beta', beta, rec['roots_of_p_Re_>_-1/beta'], hv[1], 'pos', pos, flush=True)
    # (E4)
    beta, t = POINT; r = 1 / beta; a = t / beta
    eps, delta = fl.eps_delta(beta, t)
    assert ctln.legal(eps, delta)
    chia = ctln.charpoly([[F(D[i][j]) + a for j in range(m)] for i in range(m)])      # det(zI - D - aJ), hi->lo, deg m
    print('chi_a done', round(time.time() - t0, 1), 's', flush=True)
    zchia = chia + [F(0)]                                                               # deg m+1
    chiLp = [F(0)] + chiL                                                               # pad to deg m+1

    def gk(kk):
        psi = [x - a * kk * y for x, y in zip(zchia, chiLp)]
        j0 = m - d
        assert all(v == 0 for v in psi[len(psi) - j0:]), 'z^(m-d) does not divide Psi'
        g = psi[:len(psi) - j0]
        return g                                                                        # deg d+1 = 7, monic

    def hk(kk):
        g = gk(kk); dg = len(g) - 1
        h = [((-1) ** dg) * v for v in shift(g, r)]
        assert h[0] == 1
        return h
    kk_sym = sp.symbols('k')
    polys = []
    for j in range(1, d + 2):
        xs = list(range(0, j + 1)); ys = [minors(hk(F(v)))[j - 1] for v in xs]
        P = sp.interpolate(list(zip([sp.Integer(v) for v in xs], [sp.Rational(y.numerator, y.denominator) for y in ys])), kk_sym)
        P = sp.Poly(sp.expand(P), kk_sym)
        for v in range(j + 1, j + 5):
            assert P.eval(v) == sp.Rational(*(lambda q: (q.numerator, q.denominator))(minors(hk(F(v)))[j - 1]))
        polys.append(P)
    k0 = 0; rec_m = []
    for j, P in enumerate(polys):
        if P.degree() <= 0:
            assert P.LC() > 0; kj = 0
        else:
            assert P.LC() > 0, ('minor', j + 1, 'leading coefficient not positive')
            rr = P.real_roots(); kj = (int(sp.floor(max(rr))) + 1) if rr else 0; kj = max(kj, 0)
            assert P.count_roots(kj, None) == 0 and P.eval(kj) > 0
        rec_m.append({'minor': j + 1, 'degree_in_k': P.degree(), 'k_from': kj, 'poly': str(P.as_expr())})
        k0 = max(k0, kj)
    # sanity: verdict just below k0 (if k0 > 0) is not 'stable', verdict at k0 is stable
    v_at = ctln.hurwitz_stable(hk(F(k0)))
    v_below = ctln.hurwitz_stable(hk(F(k0 - 1))) if k0 > 0 else None
    assert v_at[0] is True
    out['E4'] = {'beta': str(beta), 't': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)),
                 'g_0': [str(v) for v in gk(F(0))], 'g_1_minus_g_0': [str(x - y) for x, y in zip(gk(F(1)), gk(F(0)))],
                 'minors': rec_m, 'k0': k0, 'n0': m + k0,
                 'verdict_at_k0': v_at[1], 'verdict_at_k0_minus_1': (v_below[1] if v_below else None),
                 'statement': ('at (eps, delta) = (%s, %s) [beta = %s, t = %s, q = %s] the full support of core[S, all; k] is a '
                               'stable positive fixed point for every integer k >= %d (n >= %d); positivity by E3 + Lemma X2\''
                               % (eps, delta, beta, t, ctln.q_of(eps, delta), k0, m + k0))}
    print('E4 k0 =', k0, 'n0 =', m + k0, round(time.time() - t0, 1), 's', flush=True)
    json.dump(out, open('extend_lambda.json', 'w'), indent=1)
    # (E5)
    out['E5'] = []
    for kk in CTLN_KS:
        G = fl.build(A, [S], [kk]); n = len(G)
        t1 = time.time()
        rr = ctln.analyze_support(G, list(range(n)), eps, delta)
        pv = ctln.hurwitz_stable(hk(F(kk)))
        rec = {'k': kk, 'n': n, 'ctln_positive': rr['positive'], 'ctln_stable': rr['stable'], 'ctln_detail': rr['stability_detail'],
               'psi_verdict': pv[0], 'agree': rr['stable'] == pv[0], 'clique': rr['clique'], 'seconds': round(time.time() - t1, 1)}
        out['E5'].append(rec); print('E5', rec, flush=True)
        assert rec['agree'] and rr['positive']
        json.dump(out, open('extend_lambda.json', 'w'), indent=1)
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open('extend_lambda.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
