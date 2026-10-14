"""HEURISTIC. Survey of first-order coefficient R at moving zeros z0 with Re z0 near 1/2."""
import numpy as np, sys, gzip
from zeros_R import zeros_R
from scan_n7pos import from_d6
def a_pos(A, beta):
    m = len(A); N = np.ones((m, m)) - np.eye(m) - A.T
    C = np.eye(m) + beta * N
    a = np.linalg.det(C) * np.linalg.solve(C, np.ones(m)); return bool((a > 1e-9).all())
src = sys.argv[1]
rows = []
if src == 'n7S':
    for line in gzip.open('/home/user/CLTN/n7/n7_survivors_all.txt.gz', 'rt'):
        s, e, f = line.split()
        if f == 'S': rows.append(s)
elif src == 'n7B':
    k = int(sys.argv[2])
    for i, line in enumerate(gzip.open('/home/user/CLTN/n7/n7_survivors_all.txt.gz', 'rt')):
        s, e, f = line.split()
        if f == 'B' and i % k == 0: rows.append(s)
bins = {}
best = []
for s in rows:
    A = from_d6(s)
    for z0, R, kind in zeros_R(A):
        if kind != 'moving' or abs(z0.imag) < 1e-9: continue
        x = z0.real
        if 0.3 < x < 0.7:
            b = round(x, 2)
            bins[b] = max(bins.get(b, -9), R)
            if x < 0.52 and R > 0.02:
                best.append((R, x, s, a_pos(A, 2 + 1e-3)))
print('max R by Re z0 bin (complex moving zeros):')
for b in sorted(bins): print('  ', b, round(bins[b], 5))
best.sort(reverse=True)
for r in best[:15]: print(r)
