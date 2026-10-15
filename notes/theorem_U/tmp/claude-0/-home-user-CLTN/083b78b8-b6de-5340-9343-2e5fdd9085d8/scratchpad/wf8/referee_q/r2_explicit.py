"""r2_explicit.py -- referee: rebuild the two explicit graphs independently, check X1 on them, run ctln.analyze_support at
the stated (eps, delta), run the same on the prover's stored adjacency rows, and recompute k0 (reduced, exact) with the
explicit limit test at k0-1, k0, k0+1.  Writes r2_explicit.json."""
import sys, json, time
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q')
sys.path.insert(0, '/home/user/CLTN')
from ref_lib import *
import ctln
from fastregion import icharpoly

PQ = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/'
prov = json.load(open(PQ + 'certificates_explicit.json'))


def sink_ext(Acore, root, k):
    m = len(Acore); n = m + k
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = Acore[i][j]
    for v in range(m, n):
        for i in range(m):
            B[v][i] = 1                       # class sends to every core node
            B[i][v] = 0 if i == root else 1   # core node i sends to class iff i != c2
        for w in range(m, n):
            if w != v: B[v][w] = 1
    return B


def psi_int(D):
    m = len(D)
    a = icharpoly(D); b = icharpoly([[D[i][j] + 1 for j in range(m)] for i in range(m)])
    return sp.Poly(sum((ai - bi) * z ** (m - e) for e, (ai, bi) in enumerate(zip(a, b))), z)


def explicit_limit(B, beta):
    ps = psi_int(D_of_A(B))
    # strip zero roots
    e = ps.as_expr(); P = sp.Poly(e, z)
    while P.eval(0) == 0: P = sp.Poly(sp.cancel(P.as_expr() / z), z)
    return strip_right_of(P, sp.Rational(1) / sp.nsimplify(beta)), ps


out = []
for rec in prov:
    c = rec['c']; k = rec['k']; beta = F(rec['beta']); t = F(rec['t'])
    D, idx = arms_core(c)
    Acore = A_of_D(D); m = len(D)
    B = sink_ext(Acore, 0, k); n = len(B)
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    r = {'c': c, 'k': k, 'n': n, 'beta': str(beta), 'eps': str(eps), 'delta': str(delta),
         'eps_matches_prover': str(eps) == rec['eps'], 'delta_matches_prover': str(delta) == rec['delta']}
    # X1 + Lemmas T/Q on the explicit graph
    p, G, mm = reduced_from_c(c); d = len(c) - 1
    ps_expl = psi_int(D_of_A(B))
    ps_red = sp.Poly(z ** (k - 1) * z ** (m - d) * (G.as_expr() + k * p.as_expr()), z)
    r['X1_psi_identity'] = ps_expl == ps_red
    # sympy berkowitz cross-check of psi for the smaller graph
    if n <= 30: r['psi_sympy_crosscheck'] = (psi_sym(D_of_A(B)) == ps_expl)
    # reduced k0 (exact)
    t0 = time.time()
    k0, passes, det = k_threshold(p, G, sp.Rational(beta.denominator, beta.numerator))
    r['k0_referee'] = k0; r['k0_secs'] = round(time.time() - t0, 1)
    # explicit limit test at k0-1, k0, k0+1
    lim = {}
    for kq in [k0 - 1, k0, k0 + 1]:
        Bq = sink_ext(Acore, 0, kq)
        v, _ = explicit_limit(Bq, sp.Rational(beta.numerator, beta.denominator))
        lim[kq] = v
    r['explicit_limit_test'] = {str(a): b for a, b in lim.items()}
    # ctln on own graph
    t0 = time.time()
    a = ctln.analyze_support(B, list(range(n)), eps, delta)
    r['own_graph_ctln'] = {kq: a[kq] for kq in ['q', 'is_fixed_point', 'positive', 'stable', 'stability_detail', 'clique', 'stable_fixed_point', 'agrees']}
    r['own_graph_legal'] = ctln.legal(eps, delta)
    r['own_graph_x_min'] = str(min(F(v) for v in a['x']))
    r['own_ctln_secs'] = round(time.time() - t0, 1)
    # ctln on prover's stored adjacency
    rows = rec['adjacency_rows']
    rows = eval(rows) if isinstance(rows, str) else rows
    Bp = [[int(ch) for ch in row] for row in rows]
    assert len(Bp) == n
    ap = ctln.analyze_support(Bp, list(range(n)), eps, delta)
    r['prover_graph_ctln'] = {kq: ap[kq] for kq in ['is_fixed_point', 'positive', 'stable', 'stability_detail', 'clique', 'stable_fixed_point']}
    # cheap isomorphism evidence: same sorted (outdeg, indeg) multiset and same charpoly of A
    deg = lambda M: sorted((sum(M[i]), sum(M[j][i] for j in range(len(M)))) for i in range(len(M)))
    r['same_degree_multiset'] = deg(B) == deg(Bp)
    r['same_charpoly_A'] = icharpoly(B) == icharpoly(Bp)
    r['same_x_min'] = str(min(F(v) for v in ap['x'])) == r['own_graph_x_min']
    out.append(r)
    print(json.dumps(r, default=str), flush=True)

json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q/r2_explicit.json', 'w'), indent=1, default=str)
