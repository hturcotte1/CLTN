"""Check that each family witness graph is literally the claimed construction (base rows + sink classes)."""
import sys, json
sys.path.insert(0, '.')
from verify_instances_lib import my_d6
TH = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/'
def check(A, base_rows, classes):
    m = len(base_rows); base = [[int(c) for c in r] for r in base_rows]
    n = len(A); assert n == m + sum(k for S, k in classes), (n, m, [k for S,k in classes])
    for i in range(m):
        for j in range(m): assert A[i][j] == base[i][j]
    v = m
    blocks = []
    for S, k in classes:
        blocks.append((S, list(range(v, v + k)))); v += k
    for S, nodes in blocks:
        for x in nodes:
            for i in range(m):
                assert A[i][x] == S[i] and A[x][i] == 1
            for y in range(m, n):
                if y != x: assert A[x][y] == 1
    return True
c = json.load(open(TH + 'family_Q_certificate.json'))
w = c['Q5_witness']; print('Q5', check(my_d6(w['d6']), c['B_rows'], [([int(x) for x in c['S']], w['k'])]))
c = json.load(open(TH + 'family_Q2_certificate.json'))
print(c['family']); print(c['classes'])
for w in c['W5_witnesses']:
    cl = [([int(x) for x in S], r * w['k']) for S, r in c['classes']]
    print('Q2', w['n'], check(my_d6(w['d6']), c['B_rows'], cl))
c = json.load(open(TH + 'family_Y_certificate.json'))
print(c['family']); print(c['classes'])
w = c['Y5_witness']
cl = [([int(x) for x in S], r * w['k']) for S, r in c['classes']]
print('Y5', w['n'], check(my_d6(w['d6']), c['C8_rows'], cl))
