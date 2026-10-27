"""Calibrate the three root counters on polynomials with KNOWN exact roots (built from the roots)."""
import sys, random, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8/v_roots')
from fractions import Fraction as F
from vroots8 import routh_rhp, cauchy_rhp, sympy_rhp, charpoly, peval, det, refine_root, d6_decode, d6_encode

def mulp(a, b):  # high->low
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r

random.seed(20261027)
rep = {'random_known_roots': {}, 'imag_axis': [], 'repeated': [], 'charpoly': {}, 'refine': {}}
stats = {'n': 0, 'routh_ok': 0, 'cauchy_ok': 0, 'sympy_ok': 0, 'sympy_run': 0, 'fail': [], 'by_degree': {}, 'by_rhp': {}}
t0 = time.time()
for trial in range(600):
    deg = random.choice([3, 4, 5, 6, 7, 8, 8, 8])
    p = [F(1)]; rhp = 0; d = 0
    while d < deg:
        if deg - d >= 2 and random.random() < 0.5:
            a = F(random.choice([-1, 1]) * random.randint(1, 40), random.randint(1, 30))
            b = F(random.randint(1, 50), random.randint(1, 20))
            p = mulp(p, [F(1), -2 * a, a * a + b * b]); d += 2      # roots a +- i b
            if a > 0: rhp += 2
        else:
            a = F(random.choice([-1, 1]) * random.randint(1, 40), random.randint(1, 30))
            p = mulp(p, [F(1), -a]); d += 1                          # real root a
            if a > 0: rhp += 1
    k = F(random.choice([-1, 1]) * random.randint(1, 9), random.randint(1, 9))
    p = [k * x for x in p]
    r1, _ = routh_rhp(p); r2, _ = cauchy_rhp(p)
    run_s = (trial % 3 == 0)
    r3, _ = sympy_rhp(p) if run_s else (None, None)
    stats['n'] += 1
    stats['by_degree'][deg] = stats['by_degree'].get(deg, 0) + 1
    stats['by_rhp'][rhp] = stats['by_rhp'].get(rhp, 0) + 1
    stats['routh_ok'] += (r1 == rhp); stats['cauchy_ok'] += (r2 == rhp)
    if run_s:
        stats['sympy_run'] += 1; stats['sympy_ok'] += (r3 == rhp)
    if r1 != rhp or r2 != rhp or (run_s and r3 != rhp):
        stats['fail'].append({'p': [str(x) for x in p], 'true_rhp': rhp, 'routh': r1, 'cauchy': r2, 'sympy': r3})
stats['seconds'] = round(time.time() - t0, 1)
rep['random_known_roots'] = stats
print(json.dumps(stats, indent=1))

# roots ON the imaginary axis must come out undecided (None) from Cauchy and sympy, and Routh must hit a zero pivot
for p, what in (([F(1), F(0), F(1)], 'mu^2+1'),
                (mulp([F(1), F(0), F(4)], [F(1), F(3), F(2)]), '(mu^2+4)(mu+1)(mu+2)'),
                (mulp([F(1), F(0)], [F(1), F(2), F(5)]), 'mu (mu^2+2mu+5)'),
                (mulp(mulp([F(1), F(0), F(9)], [F(1), F(1)] * 1), mulp([F(1), F(2), F(2)], [F(1), F(1, 2)])), '(mu^2+9)(mu+1)(mu^2+2mu+2)(mu+1/2)')):
    r1, n1 = routh_rhp(p); r2, n2 = cauchy_rhp(p); r3, n3 = sympy_rhp(p)
    rep['imag_axis'].append({'poly': what, 'routh': r1, 'cauchy': r2, 'sympy': r3, 'cauchy_note': n2, 'sympy_note': n3})
    assert r2 is None and r3 is None and r1 is None, (what, r1, r2, r3)
print('imaginary-axis cases: all three undecided, OK')

# repeated roots off the axis
for p, rhp, what in ((mulp(mulp([F(1), F(1)], [F(1), F(1)]), mulp([F(1), F(1), F(1)], [F(1), F(1), F(1)])), 0, '(mu+1)^2 (mu^2+mu+1)^2'),
                     (mulp(mulp([F(1), F(-1)], [F(1), F(-1)]), [F(1), F(2)]), 2, '(mu-1)^2 (mu+2)'),
                     (mulp(mulp([F(1), F(-1), F(3)], [F(1), F(-1), F(3)]), [F(1), F(1, 3)]), 4, '(mu^2-mu+3)^2 (mu+1/3)'),
                     (mulp(mulp(mulp([F(1), F(2)], [F(1), F(2)]), [F(1), F(2)]), mulp([F(1), F(2)], [F(1), F(1), F(5)])), 0, '(mu+2)^4 (mu^2+mu+5)')):
    r = (routh_rhp(p)[0], cauchy_rhp(p)[0], sympy_rhp(p)[0])
    rep['repeated'].append({'poly': what, 'true_rhp': rhp, 'routh_cauchy_sympy': list(r)})
    assert all(v is None or v == rhp for v in r) and r[1] == rhp and r[2] == rhp, (what, r)
print('repeated-root cases OK')

# charpoly against matrices with known spectrum
M = [[F(2), F(1)], [F(0), F(-3)]]
assert charpoly(M) == [1, 1, -6]
import random as _r
_r.seed(5)
n = 8
diag = [F(_r.randint(-9, 9), _r.randint(1, 4)) for _ in range(n)]
T = [[diag[i] if i == j else (F(_r.randint(-5, 5)) if j > i else F(0)) for j in range(n)] for i in range(n)]
# similarity transform by a random unimodular matrix to get a dense matrix with the same spectrum
P = [[F(int(i == j)) for j in range(n)] for i in range(n)]
for _ in range(20):
    i, j = _r.sample(range(n), 2); f = F(_r.randint(-3, 3))
    for k in range(n): P[i][k] += f * P[j][k]
Pinv = [[F(int(i == j)) for j in range(n)] for i in range(n)]
# invert P by solving
from vroots8 import solve
cols = [solve(P, [F(int(i == j)) for i in range(n)])[1] for j in range(n)]
Pinv = [[cols[j][i] for j in range(n)] for i in range(n)]
def mm(X, Y): return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
Mdense = mm(mm(P, T), Pinv)
expected = [F(1)]
for d in diag: expected = mulp(expected, [F(1), -d])
got = charpoly(Mdense)
assert got == expected, 'charpoly mismatch'
rhp_true = sum(1 for d in diag if d > 0)
rc = (routh_rhp(got)[0], cauchy_rhp(got)[0], sympy_rhp(got)[0])
rep['charpoly'] = {'n': n, 'diag': [str(d) for d in diag], 'charpoly_matches_product_of_(mu-d)': True,
                   'true_rhp': rhp_true, 'routh_cauchy_sympy': list(rc)}
assert all(v == rhp_true for v in rc), rc
print('charpoly on dense 8x8 with known spectrum OK', rep['charpoly'])

# refine_root on a known algebraic root: sqrt(2) in [1, 2], and a rational root 3/2 of 2x-3
a, b = refine_root([1, 0, -2], F(1), F(2), lambda a, b: b - a < F(1, 10**12))
assert a * a < 2 < b * b and b - a < F(1, 10**12)
r1, r2 = refine_root([2, -3], F(1), F(2), lambda a, b: b - a < F(1, 10**6))
assert r1 == r2 == F(3, 2)
rep['refine'] = {'sqrt2_interval': [str(a), str(b)], 'rational_root_detected': str(r1)}
# digraph6 round trip on a random 8-node digraph
A = [[0 if i == j else _r.randint(0, 1) for j in range(8)] for i in range(8)]
s = d6_encode(A); assert d6_decode(s) == A and len(s) == 13 and s[:2] == '&G'
rep['d6_roundtrip'] = s
json.dump(rep, open('/home/user/CLTN/notes/verify_n8/v_roots/calibration.json', 'w'), indent=1)
print('CALIBRATION DONE')
