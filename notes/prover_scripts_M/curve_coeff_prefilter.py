"""EXACT early exit for the curve test: if no t > 0 has all coefficients c_k(t) > 0 (k = 1..m) of det(mu I + L(t)),
then by (b1) the curve is never stable. Measures how often this decides 'never' and the time saved (sample)."""
import sys, time, random, gzip
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
from curve1d import coeff_polys, curve_test
from fastregion import _hi_to_lo, _cells, peval
from certify_M_violation import from_d6
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ
def coeff_never(A):
    m = len(A); P = coeff_polys(A)
    if any(not P[k] for k in range(1, m + 1)): return True
    iso = dup_isolate_real_roots_list([_hi_to_lo(P[k]) for k in range(1, m + 1) if len(P[k]) > 1], ZZ, inf=0, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == 0)]
    for _, _, lo, hi, s in _cells(iso, F(0)):
        if all(peval(P[k], s) > 0 for k in range(1, m + 1)): return False
    return True
P = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/n7pos/n7_pos_all_PU.txt.gz'
random.seed(7); lines = []
for i, l in enumerate(gzip.open(P, 'rt')):
    sp = l.split()
    if int(sp[1]) >= 22 and random.random() < 0.0005: lines.append(sp[0])
lines = lines[:1500]
t0 = time.time(); nc = 0; bad = 0
for s in lines:
    A = from_d6(s)
    if coeff_never(A): nc += 1
t1 = time.time()
for s in lines[:300]:
    A = from_d6(s)
    if coeff_never(A) and curve_test(A)['verdict'] != 'never': bad += 1
print('sample', len(lines), 'decided never by coefficient positivity alone:', nc, f'({nc/len(lines):.1%})',
      'sec/graph', (t1 - t0) / len(lines), 'consistency failures (should be 0):', bad)
