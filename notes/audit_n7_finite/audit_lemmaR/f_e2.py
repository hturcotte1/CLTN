"""e_2 audit (exact). 2,000 random 7-node digraphs x 1 random legal point each (plus 3 extra points for 500 of them):
e_2(K) is computed independently from ctln.weight_matrix: K = (I - W)/eps, e_2 = (tr(K)^2 - tr(K^2))/2, and also
as the sum of the 2x2 principal minors; compared with the closed form
   B(1+2t) - Nn(beta-1)(beta+1+2t) + O(1 - t(beta-2))   and   t Phi(beta) + B + O - Nn(beta^2-1).
Also: e_2 <= (B - Nn)(1 + 2t) at every legal point (the inequality behind 'B > Nn'), the m=7 identity
edges = 21 + B - Nn, and the edge column of n7/n7_pos_all_PU.txt.gz against from_d6 on 20,000 lines."""
import sys, json, random, gzip, itertools
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from n7.curve_test import from_d6
rng = random.Random(31337)
def rand_t():
    return F(rng.randint(1, 10 ** 6), 10 ** rng.randint(0, 4)) / 10
fails = []; n = 0; ineq_fail = 0
for g in range(2000):
    p = rng.choice([0.2, 0.4, 0.5, 0.6, 0.75, 0.9, 0.97])
    A = [[0 if i == j else int(rng.random() < p) for j in range(7)] for i in range(7)]
    B = sum(1 for i, j in itertools.combinations(range(7), 2) if A[i][j] and A[j][i])
    Nn = sum(1 for i, j in itertools.combinations(range(7), 2) if not A[i][j] and not A[j][i])
    O = 21 - B - Nn
    edges = sum(map(sum, A))
    if edges != 2 * B + O or edges != 21 + B - Nn: fails.append(("edge identity", g))
    for k in range(1 if g >= 500 else 4):
        t = rand_t(); eta = F(rng.randint(1, 10 ** 6), 10 ** rng.randint(3, 12))
        beta = 2 + 1 / t + eta; eps = 1 / (t + 1); delta = (beta - 1) * eps
        assert ctln.legal(eps, delta)
        W = ctln.weight_matrix(A, eps, delta)
        K = [[((1 if i == j else 0) - W[i][j]) / eps for j in range(7)] for i in range(7)]
        tr = sum(K[i][i] for i in range(7)); tr2 = sum(K[i][j] * K[j][i] for i in range(7) for j in range(7))
        e2a = (tr * tr - tr2) / 2
        e2b = sum(K[i][i] * K[j][j] - K[i][j] * K[j][i] for i, j in itertools.combinations(range(7), 2))
        closed = B * (1 + 2 * t) - Nn * (beta - 1) * (beta + 1 + 2 * t) + O * (1 - t * (beta - 2))
        Phi = 42 - beta * (2 * Nn + O)
        closed2 = t * Phi + B + O - Nn * (beta * beta - 1)
        n += 1
        if not (e2a == e2b == closed == closed2): fails.append(("e2 mismatch", g, str(t), str(beta)))
        if not e2a <= (B - Nn) * (1 + 2 * t): ineq_fail += 1
        # sanity: K entries are exactly I + beta N + tJ with N[i][j] = 1 iff j -/-> i
        for i in range(7):
            for j in range(7):
                exp = (1 + t) if i == j else (t + beta * (0 if A[j][i] else 1))
                if K[i][j] != exp: fails.append(("K form", g)); break
# edge column of the positivity file
cnt = 0; colfail = 0
with gzip.open('/home/user/CLTN/n7/n7_pos_all_PU.txt.gz', 'rt') as fh:
    for line in fh:
        cnt += 1
        if cnt % 400 != 1: continue
        sp_ = line.split(); A = from_d6(sp_[0])
        if sum(map(sum, A)) != int(sp_[1]) or any(A[i][i] for i in range(7)): colfail += 1
        if cnt > 400 * 20000: break
out = {"graphs": 2000, "points": n, "fails": fails[:20], "n_fails": len(fails), "ineq_e2_le_(B-Nn)(1+2t)_fail": ineq_fail,
       "pos_file_lines_checked": (min(cnt, 400 * 20000) + 399) // 400, "pos_file_edge_col_fail": colfail}
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/f_e2.json', 'w'), indent=1)
print(json.dumps(out))
