"""r4_validate.py -- referee: validate referee tools on known small-n numbers, using the repo's Geneson graph
(test_ctln.py GENESON_A) rather than the prover's hard-coded G1.  Writes r4_validate.json."""
import sys, json, itertools
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q')
sys.path.insert(0, '/home/user/CLTN')
from ref_lib import *
import ctln

GN = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]
GA = [[0]*6 for _ in range(6)]
for a in range(6):
    for b in range(6):
        if a != b: GA[b][a] = 1 - GN[a][b]
PG1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
iso = None
for perm in itertools.permutations(range(6)):
    if all(PG1[i][j] == GA[perm[i]][perm[j]] for i in range(6) for j in range(6)): iso = perm; break
out = {'prover_G1_isomorphic_to_repo_Geneson': iso is not None, 'perm': iso}


def lim(A, beta):
    ps = psi_sym(D_of_A(A)); P = ps
    while P.eval(0) == 0: P = sp.Poly(sp.cancel(P.as_expr() / z), z)
    return strip_right_of(P, 1 / sp.Rational(beta))


def posi(A, beta):
    D = D_of_A(A); m = len(D); r = 1 / F(beta)
    u = ctln.solve_exact([[F(D[j][i]) + (r if i == j else 0) for j in range(m)] for i in range(m)], [F(1)] * m)
    return u is not None and min(u) > 0 and sum(u) > 1


bl = 1 + sp.Rational(561, 481); bh = 1 + sp.Rational(7, 6)
out['V1'] = {'lim_at_561/481': lim(GA, bl), 'pos_at_561/481': posi(GA, F(1042, 481)), 'lim_at_7/6': lim(GA, bh)}
# direct ctln at a legal point on the ray beta = 1042/481, t = 10^5
beta = F(1042, 481); t = F(10 ** 5)
eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
rec = ctln.analyze_support(GA, list(range(6)), eps, delta)
out['V1_ctln_point'] = {'q': rec['q'], 'stable_fixed_point': rec['stable_fixed_point'], 'detail': rec['stability_detail']}
# sink-class F_1, F_2 with S = all but the image of prover's node 5
S = [1] * 6; S[iso[5]] = 0
def ext(A, S, k):
    m = len(A); n = m + k
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    for v in range(m, n):
        for i in range(m): B[v][i] = 1; B[i][v] = S[i]
        for w in range(m, n):
            if w != v: B[v][w] = 1
    return B
F1 = ext(GA, S, 1); F2 = ext(GA, S, 2)
out['V2'] = {'lim_2.2746': lim(F1, sp.Rational(22746, 10000)), 'lim_2.2748': lim(F1, sp.Rational(22748, 10000))}
out['V3'] = {'pos_2.3247': posi(F2, F(23247, 10000)), 'pos_2.3248': posi(F2, F(23248, 10000)),
             'lim_2.3248': lim(F2, sp.Rational(23248, 10000))}
print(json.dumps(out, default=str, indent=1))
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q/r4_validate.json', 'w'), indent=1, default=str)
