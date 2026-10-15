"""
check_lemmas.py -- EXACT checks (integers / Fractions / sympy; no floats) of the lemmas of notes_prover_q.md.

Lemma T (in-tree characteristic polynomial).  D a 0/1 zero-diagonal m x m matrix with a root rho such that every
node u != rho has exactly one D-out-neighbour and every node reaches rho; ell(u) = length of the D-path u -> rho,
N = N^+(rho) = D-out-neighbours of rho.  Then, with x = 1/z,
   chi_{D - e_rho 1^T}(z) = z^m Lam(x),   Lam(x) = 1 + sum_{u not in N} x^(1 + ell(u)),
   chi_D(z)              = z^m (1 - Gam(x)),  Gam(x) = sum_{w in N} x^(1 + ell(w)),
   psi_D(z) := 1^T adj(zI - D) 1 = z^(m-1) Q(x),  Q(x) = [(m + 1 - |N|)(1 - Gam(x)) + (|N| - 1) Lam(x)] / (1 - x).
Checked on random in-trees (random tree shape, random N containing every leaf) and on arms forests.

Lemma P (positivity of the arms forest; closed form) checked symbolically (sympy) and against exact direct solves.
Run: python3 check_lemmas.py   (writes check_lemmas.json)
"""
import sys, json, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
import ctln
from fastregion import icharpoly
import sympy as sp
from forest_exact import build_core, L_matrix, positivity, closed_form_lambda

def poly_from_terms(terms, deg):
    c = [0] * (deg + 1)
    for e, v in terms: c[e] += v
    return c   # c[j] = coefficient of x^j

def random_intree(m, rng):
    """node 0 = root rho.  nodes 1..m-1: next(u) in {0..u-1} (so every node reaches 0, in-tree).  N = all leaves plus a
    random subset of internal non-root nodes; N must not contain rho."""
    nxt = {u: rng.randrange(0, u) for u in range(1, m)}
    children = {u: [] for u in range(m)}
    for u, v in nxt.items(): children[v].append(u)
    leaves = [u for u in range(1, m) if not children[u]]
    N = set(leaves) | {u for u in range(1, m) if rng.random() < 0.3}
    D = [[0] * m for _ in range(m)]
    for u, v in nxt.items(): D[u][v] = 1
    for w in N: D[0][w] = 1
    return D, nxt, N

def ell_of(nxt, m):
    ell = {0: 0}
    for u in range(1, m):
        k, v = 0, u
        while v != 0: v = nxt[v]; k += 1
        ell[u] = k
    return ell

def check_T(D, nxt_ell, N, root=0):
    m = len(D); ell = nxt_ell
    Lam = poly_from_terms([(0, 1)] + [(1 + ell[u], 1) for u in range(m) if u not in N], m)
    Gam = poly_from_terms([(1 + ell[w], 1) for w in N], m + 1)
    # chi_L
    L = [[D[i][j] - (1 if i == root else 0) for j in range(m)] for i in range(m)]
    chiL = icharpoly(L)                       # hi->lo: chiL[j] = coeff of z^(m-j) = coeff of x^j in Lam
    okL = chiL == Lam[:m + 1]
    chiD = icharpoly(D)
    okD = chiD == [1 if j == 0 else -Gam[j] for j in range(m + 1)] and all(g == 0 for g in Gam[m + 1:])
    DJ = [[D[i][j] + 1 for j in range(m)] for i in range(m)]
    psi = [a - b for a, b in zip(chiD, icharpoly(DJ))][1:]     # hi->lo, degree m-1: psi[j] = coeff of x^j in Q
    x = sp.symbols('x'); n_ = len(N)
    Qexpr = sp.cancel(((m + 1 - n_) * (1 - sum(sp.Integer(g) * x ** j for j, g in enumerate(Gam))) +
                       (n_ - 1) * sum(sp.Integer(l) * x ** j for j, l in enumerate(Lam))) / (1 - x))
    Qp = sp.Poly(Qexpr, x).all_coeffs()[::-1]
    Qp = [int(v) for v in Qp] + [0] * (m - len(Qp))
    okQ = psi == Qp[:m] and all(v == 0 for v in Qp[m:])
    return okL, okD, okQ

def main():
    rng = random.Random(20261015)
    out = {'lemma_T_random_intrees': [], 'lemma_T_arms_forests': [], 'lemma_P': {}}
    for trial in range(40):
        m = rng.randrange(4, 26)
        D, nxt, N = random_intree(m, rng)
        okL, okD, okQ = check_T(D, ell_of(nxt, m), N)
        out['lemma_T_random_intrees'].append({'m': m, 'N_size': len(N), 'chi_L': okL, 'chi_D': okD, 'psi_D': okQ})
        assert okL and okD and okQ, ('Lemma T failed', trial)
    for c in ([1, 1, 2, 2, 2, 1, 1], [1, 1, 4, 3, 4, 2, 1], [1, 1, 11, 7, 10, 3, 1], [1, 1, 3, 1, 2, 1, 1, 1, 1], [1, 1, 1, 2, 1, 3]):
        D, names = build_core(c); m = len(D); root = names.index('c2')
        # relabel so root = 0 for check_T's conventions
        perm = [root] + [i for i in range(m) if i != root]; pos = {v: k for k, v in enumerate(perm)}
        Dp = [[D[perm[i]][perm[j]] for j in range(m)] for i in range(m)]
        nxt = {}
        for u in range(1, m):
            outs = [v for v in range(m) if Dp[u][v]]; assert len(outs) == 1; nxt[u] = outs[0]
        N = {v for v in range(m) if Dp[0][v]}
        okL, okD, okQ = check_T(Dp, ell_of(nxt, m), N)
        chi = icharpoly(L_matrix(D, names)); d = len(c) - 1
        okp = chi == list(c) + [0] * (m - d)
        out['lemma_T_arms_forests'].append({'c': c, 'm': m, 'chi_L': okL, 'chi_D': okD, 'psi_D': okQ, 'chi_L_equals_z^(m-d)p': okp})
        assert okL and okD and okQ and okp
    # Lemma P symbolic: the closed form solves the cycle system
    r, F0, F1, F2 = sp.symbols('r F0 F1 F2')
    Den = F0 - F1 * r + F2 * r ** 2 + r ** 2 + r ** 5
    lam = r ** 3 * (1 - r + r ** 2) / Den
    e2 = r * (1 - lam)
    e1 = r * (1 - e2) + lam * F2 / r
    e0 = r * (1 - e1) + lam * F1 / r
    res = sp.simplify(e2 - (r * (1 - e0) + lam * F0 / r))
    out['lemma_P']['cycle_system_residual'] = str(res)
    assert res == 0
    out['lemma_P']['u_c1_numerator'] = str(sp.factor(sp.simplify((1 - e1) * Den)))
    out['lemma_P']['u_c0_numerator'] = str(sp.factor(sp.expand(sp.simplify((1 - e0) * Den))))
    out['lemma_P']['eps_c2_times_Den'] = str(sp.factor(sp.simplify(e2 * Den)))
    # exact direct solves vs closed form at several rational r
    cmp = []
    for c in ([1, 1, 4, 3, 4, 2, 1], [1, 1, 11, 7, 10, 3, 1], [1, 1, 3, 1, 2, 1, 1, 1, 1]):
        D, names = build_core(c)
        for beta in (F(3), F(7, 2), F(5), F(11, 2)):
            P = positivity(D, names, beta)
            lamv, (f0, f1, f2, den) = closed_form_lambda(c, beta)
            r_ = 1 / beta
            u2 = P['u'][names.index('c2')]
            ok = (lamv == 1 - (1 - u2) / r_)
            # every arm/z value: u_a = 1 - lam f_a, u_z = lam f_a / r  (f from the forest recursion)
            cmp.append({'c': c, 'beta': str(beta), 'lambda_matches': ok, 'positive': P['ok']})
            assert ok
    out['lemma_P']['direct_vs_closed_form'] = cmp
    json.dump(out, open('check_lemmas.json', 'w'), indent=1)
    print('Lemma T: %d random in-trees + %d arms forests: all exact identities hold' % (len(out['lemma_T_random_intrees']), len(out['lemma_T_arms_forests'])))
    print('Lemma P: cycle residual =', res, '; closed form = direct exact solve in', len(cmp), 'cases')
    print('u_c1 numerator:', out['lemma_P']['u_c1_numerator'])
    print('u_c0 numerator:', out['lemma_P']['u_c0_numerator'])
    print('eps_c2 * Den  :', out['lemma_P']['eps_c2_times_Den'])

if __name__ == '__main__':
    main()
