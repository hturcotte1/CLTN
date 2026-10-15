"""
negative.py -- EXACT negative results at a fixed parameter point (default (eps, delta) = (1/21, 11/210), beta = 21/10, t = 20).

 (N1) Geneson's clique-class family F_k = G1[S, all; k] (S = paper vertices 1..5): via Lemma P, F_k is stable at the point iff
      the Hurwitz minors Delta_j(k) of (-1)^7 Psi_k(-r - x) are all > 0.  Exact: all minors > 0 on [0, k*) and minor 6 (negative
      leading coefficient) is < 0 on [41, inf); so F_k is stable exactly for k = 0..40 and unstable for every k >= 41
      (k = 39, 40 explicit; the known data 'stable k <= 38, unstable 41, 45, 60' is reproduced).
      Mechanism (Lemma E): chi_L = z^4 (z^2+z+1) has the roots omega with Re = -1/2 < -1/beta at EVERY legal beta.
 (N2) Lemma U (uniform blow-ups): for the core R = &F^U|m^QTq? and the point above, the N-fold closed-twin blow-up R^(N)
      (every node -> bidirected N-clique) is stable iff N <= N*, computed exactly; spec(aJ + D) has an eigenvalue with
      negative real part (Perron-Frobenius argument, notes Lemma U) and Re(N mu) -> -inf.
      Explicit-graph check with ctln.analyze_support at N* (stable) and N* + 1 (unstable).
Writes negative_certificate.json.
"""
import json
import sympy as sp
from fractions import Fraction as F
import fixlib as fl

BETA = F(21, 10); T = F(20)
G1 = [[0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0]]
SF = (1, 1, 1, 1, 1, 0)


def blowup_uniform(A, N):
    m = len(A); n = m * N
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            a, b = i // N, j // N
            G[i][j] = 1 if a == b else A[a][b]
    return G


def main():
    cert = {'point': {'beta': str(BETA), 't': str(T)}}
    eps, delta = fl.eps_delta(BETA, T)
    cert['point'].update(eps=str(eps), delta=str(delta), q=str(fl.ctln.q_of(eps, delta)))
    # ---------------- (N1)
    co, mins = fl.fixed_point_minors(G1, [SF], [1], BETA, T)
    P0, Ps = fl.psi_parts(G1, [SF])
    cert['N1_chi_L'] = str(sp.factor(Ps[0]))
    rec = []
    kstar = None
    for j, mn in enumerate(mins):
        P = sp.Poly(mn, fl.k)
        rr = [r for r in P.real_roots() if r >= 0] if P.degree() > 0 else []
        first = min(rr) if rr else None
        rec.append({'minor': j + 1, 'degree': P.degree(), 'LC_sign': int(sp.sign(P.LC())), 'value_at_0_positive': bool(P.eval(0) > 0),
                    'smallest_root_>=0': (str(first) + ' ~ ' + str(sp.N(first, 12))) if first is not None else None})
        assert P.eval(0) > 0
        if first is not None:
            kstar = first if kstar is None else min(kstar, first)
    cert['N1_minors'] = rec
    cert['N1_kstar'] = str(sp.N(kstar, 15))
    assert 40 < kstar < 41
    # minor 6 is negative on [41, inf):
    P6 = sp.Poly(mins[5], fl.k)
    assert P6.LC() < 0 and P6.eval(41) < 0 and P6.count_roots(41, None) == 0
    cert['N1_minor6_negative_on_[41,inf)'] = True
    expl = []
    for kk in (38, 39, 40, 41, 45, 60):
        v = fl.stable_at_k(G1, [SF], [1], BETA, T, kk)
        expl.append({'k': kk, 'psi_hurwitz': v[0], 'detail': v[1]})
    for kk in (40, 41):
        G = fl.build(G1, [SF], [kk]); n = len(G)
        r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
        expl.append({'k': kk, 'n': n, 'ctln_positive': r['positive'], 'ctln_stable': r['stable'], 'detail': r['stability_detail']})
    cert['N1_explicit'] = expl
    assert expl[-2]['ctln_stable'] is True and expl[-1]['ctln_stable'] is False
    cert['N1_conclusion'] = 'F_k is stable at (1/21, 11/210) for k = 0..40 and unstable for every k >= 41 (all positive)'
    # ---------------- (N2) Lemma U on R
    R = fl.from_d6('&F^U|m^QTq?'); m = len(R)
    D = sp.Matrix(fl.nonedge(R)); av = sp.Rational(str(T / BETA)); rv = sp.Rational(str(1 / BETA))
    chi = sp.expand((fl.z * sp.eye(m) - D - av * sp.ones(m, m)).det(method='berkowitz'))
    # trace / Perron data
    cert['N2_trace_aJ+D'] = str(av * m)
    # stability of R^(N): spec(aJ + D') \ {0} = N spec(aJ + D); Re(N mu) > -r  <=>  roots of chi(z) have Re > -r/N
    def stable_N(N):
        rn = rv / N
        h = sp.Poly(sp.expand((-1) ** m * chi.subs(fl.z, -rn - fl.x)), fl.x)
        cs = [F(str(sp.Rational(c))) for c in h.all_coeffs()]
        return fl.ctln.hurwitz_stable(cs)
    # is R itself (N = 1) stable here: yes (Theorem A, k = 0).  Find N*.
    N = 1; hist = []
    while True:
        v = stable_N(N); hist.append((N, v[0]))
        if v[0] is not True: break
        N += 1
        assert N < 10 ** 4
    Nstar = N - 1
    cert['N2_Nstar'] = Nstar
    cert['N2_history'] = [{'N': a, 'stable': b} for a, b in hist]
    # explicit check with ctln on R^(N*) and R^(N*+1)
    xs = []
    for NN in (Nstar, Nstar + 1):
        G = blowup_uniform(R, NN); n = len(G)
        r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
        xs.append({'N': NN, 'n': n, 'positive': r['positive'], 'stable': r['stable'], 'detail': r['stability_detail']})
    cert['N2_explicit'] = xs
    assert xs[0]['stable'] is True and xs[1]['stable'] is False
    json.dump(cert, open('negative_certificate.json', 'w'), indent=1)
    print(json.dumps({k: cert[k] for k in ('N1_kstar', 'N1_conclusion', 'N2_Nstar', 'N2_explicit')}, indent=1))


if __name__ == '__main__':
    main()
