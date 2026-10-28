"""Own exact check of why blow-ups are absent: domination-freeness (filter8's F1 definition) and the strict
limit-spectral condition alpha(A|perp) < -1/2, decided by Hurwitz minors of chi_perp(w - 1/2)."""
import sys, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from fractions import Fraction as F
from collections import Counter
R = '/home/user/CLTN/'; D = R + 'notes/verify_n8_2/completeness_twins/'

def domination_free(A):
    n = len(A)
    for j in range(n):
        for k in range(n):
            if j == k: continue
            if A[j][k] == 1 and A[k][j] == 0 and all(A[i][k] == 1 for i in range(n) if i not in (j, k) and A[i][j] == 1):
                return False
    return True

def chi_perp_shifted(A):
    """chi_{A|perp}(z) = (1/n) 1^T adj(zI - A) 1 = (1/n)(det(zI - A + J) - det(zI - A)); return coefficients (hi to lo)
    of chi_perp(w - 1/2) as Fractions, degree n-1, via interpolation-free exact route: charpoly of A and of A - J."""
    n = len(A)
    M = [[F(-A[i][j]) for j in range(n)] for i in range(n)]
    # det(zI - X) with X = A gives charpoly(A); det(zI - (A - J)) = det(zI - A + J)
    cA = charpoly([[F(A[i][j]) for j in range(n)] for i in range(n)])
    cAJ = charpoly([[F(A[i][j] - 1) for j in range(n)] for i in range(n)])
    diff = [(b - a) / n for a, b in zip(cA, cAJ)]  # degree n polynomial coefficients hi..lo; leading term cancels
    assert diff[0] == 0
    p = diff[1:]  # degree n-1, hi to lo
    assert p[0] == 1, p
    # shift z = w - 1/2: p(w - 1/2)
    m = len(p) - 1
    # Horner-style expansion
    res = [F(0)] * (m + 1)
    res[-1] = p[0]
    cur = [p[0]]
    for c in p[1:]:
        # cur(w) * (w - 1/2) + c
        new = [F(0)] * (len(cur) + 1)
        for i, a in enumerate(cur):
            new[i] += a; new[i + 1] -= a / 2
        new[-1] += c
        cur = new
    return cur

def classify_limit(A):
    c = chi_perp_shifted(A)
    mins = hurwitz_minors(c)
    if any(m == 0 for m in mins): return 'B'
    if all(m > 0 for m in mins): return 'S'
    return 'X'  # alpha > -1/2 (some root with Re > -1/2)

rows = [json.loads(l) for l in open(D + 'part2_blowups.jsonl')]
G = json.load(open(R + 'certificates/stage0_geneson_2607.21396.json'))['A_ours']
print('Geneson (alpha=-0.538 < -1/2, should be S):', classify_limit(G), 'domfree', domination_free(G))
K = [[int(i != j) for j in range(8)] for i in range(8)]; print('K8:', classify_limit(K))
cnt = Counter(); mism = []
for r in rows:
    B = d6_decode(r['blowup_d6']); df = domination_free(B); cl = classify_limit(B)
    own = 'DROPPED' if (not df or cl == 'X') else cl
    f8 = r.get('filter8_stdout', '') if not r['in_hits_own'] else 'S'
    f8flag = f8.split()[-1] if f8 else 'DROPPED'
    r['own_domination_free'] = df; r['own_limit_class'] = cl; r['own_filter_flag'] = own
    cnt[(own, f8flag)] += 1
    if own != f8flag: mism.append(r['blowup_d6'])
print('own flag vs filter8 flag (own, filter8):', dict(cnt)); print('mismatches', mism)
print('dropped reasons:', dict(Counter(('not_domfree' if not r['own_domination_free'] else 'alpha>-1/2') for r in rows if r['own_filter_flag'] == 'DROPPED')))
with open(D + 'part2_blowups_flags.jsonl', 'w') as f:
    for r in rows: f.write(json.dumps({k: r[k] for k in ('blowup_d6', 'parent_d6', 'v', 'in_hits_own', 'own_domination_free', 'own_limit_class', 'own_filter_flag')}) + '\n')
json.dump({'own_vs_filter8': {str(k): v for k, v in cnt.items()}, 'mismatches': mism}, open(D + 'part2b_summary.json', 'w'), indent=1)
