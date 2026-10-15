"""
certify_explicit.py -- EXPLICIT-GRAPH certificates for arms-forest sink-class graphs G_k (no floating point).
G_k = arms-forest core D(c) (forest_exact.build_core) + k bidirected class nodes that send to every core node and
receive from every core node except c2.  For each instance (c, k, beta, t):
 (E0) identity: explicit psi_{D'} (n x n integer charpolys) == z^(k-1) z^(m-d) (G + k p)  [Lemmas T, Q, X1 on this graph]
 (E1) limit test on the EXPLICIT n-node graph: every root of psi_{D'} has Re > -1/beta  (alpha(A|perp) < -gamma)
 (E2) positivity on the EXPLICIT graph: u = (D'^T + I/beta)^{-1} 1 > 0, 1^T u > 1  (w > 0)
 (E3) ctln.analyze_support(A, full support, eps = 1/(1+t), delta = (beta-1)/(1+t)): legal, fixed point, x > 0,
      stable with all Hurwitz minors positive (not None), graph not complete.
(E1)+(E2) use only Lemma I.1/D1/D2 and the limit rule of notes/lemma_I_spectral.md; (E3) uses nothing but ctln.
Run: python3 certify_explicit.py  (writes certificates_explicit.json)
"""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
sys.path.insert(0, '/home/user/CLTN')
import ctln
from exact_core import psi, roots_right_of, positivity, ctln_witness, is_complete
from forest_exact import build_core
from certify_forest import G_poly

def build_Gk(c, k):
    D, names = build_core(c); m = len(D); n = m + k
    ic2 = names.index('c2')
    Dk = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): Dk[i][j] = D[i][j]
    for v in range(m, n): Dk[ic2][v] = 1
    A = [[0 if i == j else 1 - Dk[i][j] for j in range(n)] for i in range(n)]
    return A, Dk, names, m

def run(c, k, beta, t):
    t0 = time.time(); beta = F(beta)
    A, Dk, names, m = build_Gk(c, k); n = len(A); d = len(c) - 1
    rec = {'c': c, 'k': k, 'n': n, 'm_core': m, 'beta': str(beta), 't': str(t)}
    ps = psi(Dk)
    G, m2 = G_poly(c); assert m2 == m
    target = [G[j] + k * c[j] for j in range(d + 1)] + [0] * (n - 1 - d)
    rec['E0_psi_identity'] = ps == target
    lt = roots_right_of(target[:d + 1], 1 / beta)            # zero roots (Re 0) stripped
    rec['E1_limit_alpha_lt_-gamma'] = [str(lt[0]), lt[1]]
    P = positivity(A, beta)
    rec['E2_positive_w'] = P['ok']; rec['E2_u_min'] = str(min(P['u'])); rec['E2_U'] = str(P['U'])
    eps, delta, r = ctln_witness(A, beta, t)
    rec.update({'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)), 'legal': ctln.legal(eps, delta),
                'is_fixed_point': r['is_fixed_point'], 'positive': r['positive'], 'stable': str(r['stable']),
                'stability_detail': r['stability_detail'], 'clique': r['clique'], 'complete': is_complete(A),
                'x_min': str(min(r['x'])) if r['x'] else None})
    rec['adjacency_rows'] = [''.join(map(str, row)) for row in A]
    rec['certified'] = bool(rec['E0_psi_identity'] and lt[0] is True and P['ok'] and rec['legal'] and r['is_fixed_point']
                            and r['positive'] and r['stable'] is True and not r['clique'])
    rec['seconds'] = round(time.time() - t0, 1)
    return rec

if __name__ == '__main__':
    jobs = [([1, 1, 5, 3, 1], 9, F(3), 10 ** 4), ([1, 1, 4, 2, 1], 26, F(13, 4), 10 ** 4)]
    out = []
    for c, k, beta, t in jobs:
        rec = run(c, k, beta, t); out.append(rec)
        print(json.dumps({x: rec[x] for x in rec if x != 'adjacency_rows'}), flush=True)
        json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/certificates_explicit.json', 'w'), indent=1)
