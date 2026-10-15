"""
certify_forest.py -- EXACT certificates for arms-forest instances (no floating point in any decision).

For an integer coefficient vector c = [1, 1, c_2, ..., c_d] (c_2 >= 1, c_j >= 1) and a rational beta* > 2:
  (C1) every root of p(z) = sum_j c_j z^(d-j) has Re > -1/beta*            [exact Hurwitz, forest_exact.roots_right_of]
  (C2) positivity of the core at beta*: u = (D^T + I/beta*)^{-1} 1 > 0, sum u > 1, u_{c2} < 1
       [closed form of Lemma P, exact; plus a direct exact solve when the core is small enough to build]
  (C3) the moving polynomial of the class family: Phi_k = z^(m-d) (G + k p), G(z) = z^d Q(1/z) (Lemma T/Q);
       k0 := least integer such that every Hurwitz minor of (-1)^d (G + k p)(-1/beta* - x) is positive for all k >= k0.
  (C4) [small cores] Lemma T and Lemma Q identities on the explicit matrix (icharpoly), and Lemma X1 on the explicit
       graph G_k for k = 1, 2 (psi of the explicit (m+k)-node non-edge matrix equals z^(k-1) Phi_k).
Conclusion (by Lemma X1/X2/X3 of notes_theory.md and Lemma I.3/I.4 of lemma_I_spectral.md, all PROVED in the repo):
for every k >= k0 the graph G_k = core + sink class of k nodes (received from every core node except c2; sends to
everyone; bidirected inside) is non-clique, positive at beta*, and alpha(A|perp) < -(1 - 1/beta*); hence its full
support is a stable positive fixed point at legal points with q arbitrarily close to beta* - 1.
Instances are CANDIDATE counterexamples pending human hand-verification (CLAUDE.md section 3).
Usage: python3 certify_forest.py   (writes certificates_forest.json)
"""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
import sympy as sp
import ctln
from fastregion import icharpoly
import family_prover as fp
from forest_exact import build_core, L_matrix, roots_right_of, positivity

def arm_spec(c):
    """list of (len j, exit, depth, count) as in forest_exact.build_core"""
    d = len(c) - 1; spec = [(2, 'c2', 1, c[2] - 1)]
    for j in range(3, d + 1):
        spec.append((j, 'c1', (j - 1) // 2, c[j]) if j % 2 else (j, 'c0', j // 2 - 1, c[j]))
    return spec

def closed_form_positivity(c, beta):
    """Lemma P, exact.  Returns dict with all u-values classes and the verdict (no core needed)."""
    r = 1 / F(beta); spec = arm_spec(c)
    levels = {}                                    # (exit, depth) -> count
    for j, ex, dep, cnt in spec:
        if cnt: levels[(ex, dep)] = levels.get((ex, dep), 0) + cnt
    for (ex, dep) in list(levels):
        for k in range(1, dep):
            if levels.get((ex, k), 0) < 1: return {'ok': False, 'reason': 'forest constraint fails'}
    Fv = {'c0': F(0), 'c1': F(0), 'c2': F(0)}; fspine = {}
    for ex in Fv:
        deps = sorted(k for (e, k) in levels if e == ex)
        if not deps: continue
        K = deps[-1]; f = F(1)                     # first (spine) arm at the deepest level is a leaf
        fspine[(ex, K)] = f
        for k in range(K - 1, 0, -1):
            f = 1 + (f + (levels.get((ex, k + 1), 0) - 1)) / r ** 2
            fspine[(ex, k)] = f
        Fv[ex] = fspine[(ex, 1)] + (levels[(ex, 1)] - 1)
    F0, F1, F2 = Fv['c0'], Fv['c1'], Fv['c2']
    Den = F0 - r * F1 + r * r * F2 + r ** 2 + r ** 5
    if Den <= 0: return {'ok': False, 'reason': 'Den <= 0'}
    lam = r ** 3 * (1 - r + r * r) / Den
    e2 = r * (1 - lam); e1 = r * (1 - e2) + lam * F2 / r; e0 = r * (1 - e1) + lam * F1 / r
    assert e2 == r * (1 - e0) + lam * F0 / r                     # Lemma P cycle system, exact
    fmax = max(fspine.values()) if fspine else F(1)
    u_cyc = {'c0': 1 - e0, 'c1': 1 - e1, 'c2': 1 - e2}
    narms = sum(levels.values())
    # sum u: cycle + arms (1 - lam f) + z's (lam f / r);  sum over arms of f = sum over exits of (spine f's + leaves)
    sum_f = sum(fspine.values()) + sum(levels[(ex, k)] - 1 for (ex, k) in levels)
    U = sum(u_cyc.values()) + narms - lam * sum_f + lam * sum_f / r
    ok = (lam > 0 and lam * fmax < 1 and min(u_cyc.values()) > 0 and u_cyc['c2'] < 1 and U > 1)
    return {'ok': ok, 'lambda': lam, 'lambda_fmax': lam * fmax, 'u_c0': u_cyc['c0'], 'u_c1': u_cyc['c1'],
            'u_c2': u_cyc['c2'], 'sum_u': U, 'F': (F0, F1, F2), 'n_arms': narms, 'm': 3 + 2 * narms}

def G_poly(c):
    """G(z) = z^d Q(1/z), Q(x) = [(m + 1 - N)(1 - Gam) + (N - 1) Lam] / (1 - x)  (Lemma Q); returns hi->lo ints."""
    x = sp.symbols('x'); d = len(c) - 1
    spec = arm_spec(c); narms = sum(s[3] for s in spec); m = 3 + 2 * narms; N = 1 + narms
    Lam = sum(sp.Integer(cj) * x ** j for j, cj in enumerate(c))
    Gam = x ** 3 + sum(sp.Integer(cnt) * x ** (1 + j) for j, ex, dep, cnt in spec)
    Q = sp.Poly(sp.cancel(((m + 1 - N) * (1 - Gam) + (N - 1) * Lam) / (1 - x)), x)
    co = Q.all_coeffs()[::-1]                              # co[j] = coeff of x^j
    co = [int(v) for v in co] + [0] * (d + 1 - len(co))
    assert len(co) == d + 1
    return co, m                                           # G(z) = sum_j co[j] z^(d-j)  (hi->lo)

def k0_threshold(G, c, beta):
    z, k, x = sp.symbols('z k x'); d = len(c) - 1
    Gz = sum(sp.Integer(a) * z ** (d - j) for j, a in enumerate(G)); Pz = sum(sp.Integer(a) * z ** (d - j) for j, a in enumerate(c))
    h = sp.Poly(sp.expand((-1) ** d * (Gz + k * Pz).subs(z, -1 / sp.Rational(beta.numerator, beta.denominator) - x)), x)
    co = [sp.expand(cc) for cc in h.all_coeffs()]
    mins = [co[0]] + fp.hurwitz_minors(co)
    k0 = 1
    for mn in mins:
        Pk = sp.Poly(mn, k)
        if Pk.is_zero or Pk.LC() <= 0: return None
        if Pk.degree() > 0:
            rr = Pk.real_roots()
            if rr: k0 = max(k0, int(sp.floor(max(rr))) + 1)
    for mn in mins:
        Pk = sp.Poly(mn, k)
        assert (Pk.degree() == 0 and Pk.LC() > 0) or (Pk.count_roots(k0, None) == 0 and Pk.eval(k0) > 0)
    return k0

def explicit_checks(c, beta, ks=(1, 2)):
    D, names = build_core(c); m = len(D); d = len(c) - 1
    chiL = icharpoly(L_matrix(D, names))
    okT = chiL == list(c) + [0] * (m - d)
    chiD = icharpoly(D); psi = [a - b for a, b in zip(chiD, icharpoly([[D[i][j] + 1 for j in range(m)] for i in range(m)]))][1:]
    G, m2 = G_poly(c); assert m2 == m
    okQ = psi == G + [0] * (m - 1 - d)                    # psi_D = z^(m-1-d) G
    P = positivity(D, names, beta)
    # Lemma X1 on explicit G_k: psi_{D'} = z^(k-1) Phi_k = z^(k-1) z^(m-d) (G + k p)   (degree m+k-1)
    okX = []
    for kk in ks:
        n = m + kk; Dk = [[0] * n for _ in range(n)]
        for i in range(m):
            for j in range(m): Dk[i][j] = D[i][j]
        ic2 = names.index('c2')
        for v in range(m, n): Dk[ic2][v] = 1            # c2 does not send to class nodes; all else sends; class sends to all
        chiDk = icharpoly(Dk); psik = [a - b for a, b in zip(chiDk, icharpoly([[Dk[i][j] + 1 for j in range(n)] for i in range(n)]))][1:]
        target = [G[j] + kk * c[j] for j in range(d + 1)] + [0] * (n - 1 - d)
        okX.append(psik == target)
    return {'m': m, 'lemma_T_explicit': okT, 'lemma_Q_explicit': okQ, 'X1_explicit_k': dict(zip(ks, okX)),
            'direct_positivity_ok': P['ok'], 'direct_u_min': str(P['u_min']), 'direct_u_c2': str(P['u_root'])}

def certify(c, beta, explicit=False):
    t0 = time.time(); beta = F(beta); rec = {'c': c, 'beta*': str(beta), 'q_limit': str(beta - 1), 'd': len(c) - 1}
    hv = roots_right_of(c, 1 / beta); rec['C1_roots_right_of_-1/beta'] = list(map(str, hv))
    pos = closed_form_positivity(c, beta)
    rec['C2_positivity_closed_form'] = {kk: (str(v) if not isinstance(v, tuple) else [str(w) for w in v]) for kk, v in pos.items()}
    G, m = G_poly(c); rec['m_core'] = m; rec['G'] = [str(g) for g in G]
    k0 = k0_threshold(G, c, beta) if hv[0] is True and pos['ok'] else None
    rec['C3_k0'] = k0; rec['n0'] = (m + k0) if k0 else None
    if explicit: rec['C4_explicit'] = explicit_checks(c, beta)
    rec['certified'] = bool(hv[0] is True and pos['ok'] and k0 is not None and (not explicit or (rec['C4_explicit']['lemma_T_explicit'] and rec['C4_explicit']['lemma_Q_explicit'] and all(rec['C4_explicit']['X1_explicit_k'].values()) and rec['C4_explicit']['direct_positivity_ok'])))
    rec['seconds'] = round(time.time() - t0, 1)
    return rec

if __name__ == '__main__':
    out = []
    small = [1, 1, 11, 7, 10, 3, 1]
    for beta in (F(4), F(5), F(11, 2), F(28, 5)):
        out.append(certify(small, beta, explicit=(beta == 5))); print(json.dumps({k: out[-1][k] for k in ('beta*', 'm_core', 'C3_k0', 'n0', 'certified', 'seconds')}), flush=True)
    json.dump(out, open('certificates_forest_small.json', 'w'), indent=1)
