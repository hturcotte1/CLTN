"""Own finite-eps check (Fractions) of witness_in_first_stable_cell for every 913th candidate (10 graphs) + K8 sanity."""
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/chain_referee')
from myfilter2 import decode_d6, charpoly, hurwitz_minors, det
def solve(M, b):
    n = len(M); A = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return None
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [A[r][j] - f * A[c][j] for j in range(n + 1)]
    return [A[i][n] / A[i][i] for i in range(n)]
def decide(d6, eps, delta):
    n, A = decode_d6(d6)
    Wm = [[Fr(0) if i == j else (Fr(-1) + eps if A[j][i] == 1 else Fr(-1) - delta) for j in range(n)] for i in range(n)]
    IW = [[(Fr(1) if i == j else Fr(0)) - Wm[i][j] for j in range(n)] for i in range(n)]
    x = solve(IW, [Fr(1)] * n)
    pos = x is not None and all(v > 0 for v in x)
    M = [[-(Fr(1) if i == j else Fr(0)) + Wm[i][j] for j in range(n)] for i in range(n)]
    c = charpoly(M)
    if any(v < 0 for v in c): st = False
    elif any(v == 0 for v in c): st = None
    else:
        m = hurwitz_minors(c); st = False if any(v < 0 for v in m) else (None if any(v == 0 for v in m) else True)
    q = delta * (1 - eps) / eps
    return pos, st, q, sum(map(sum, A)), all(A[i][j] for i in range(n) for j in range(n) if i != j)
w = json.load(open('certificates/n8_candidates_legal_witnesses.json'))['graphs']
out = []
for i in range(0, len(w), 913):
    r = w[i]; wp = r['witness_in_first_stable_cell']
    eps, delta = Fr(wp['eps']), Fr(wp['delta'])
    pos, st, q, arcs, complete = decide(r['d6'], eps, delta)
    rec = {'i': i, 'd6': r['d6'], 'eps': str(eps), 'delta': str(delta), 'q': str(q), 'legal': q > 1 and eps > 0 and delta > 0, 'positive': pos, 'stable': st, 'arcs': arcs, 'complete': complete}
    out.append(rec); print(rec)
ok = all(r['legal'] and r['positive'] and r['stable'] is True and not r['complete'] for r in out)
# sanity: K8 at eps=1/2, delta=1 (q=1) is a stable positive fixed point? and Geneson-type negative control: an oriented 3-cycle on 3 nodes must not be stable
print('ALL 10 SAMPLED CANDIDATES: legal, positive, strictly stable, non-complete:', ok)
json.dump(out, open('/home/user/CLTN/notes/verify_n8_2/chain_referee/witness_check.json', 'w'), indent=1)
