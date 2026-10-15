"""Calibrate the three root-counting methods on polynomials with KNOWN exact roots."""
import sys, random, itertools
sys.path.insert(0, '.')
from fractions import Fraction as F
from vroots import routh_rhp, cauchy_rhp, sympy_rhp, charpoly, peval, det, d6_decode
def mulp(a, b):  # high->low
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r
random.seed(1)
nfail = 0; ntest = 0
for trial in range(3000):
    deg = random.choice([3, 4, 5, 6, 7])
    p = [F(1)]; rhp = 0; d = 0
    while d < deg:
        if deg - d >= 2 and random.random() < 0.5:
            a = F(random.choice([-1, 1]) * random.randint(1, 40), random.randint(1, 30))
            b = F(random.randint(1, 50), random.randint(1, 20))
            p = mulp(p, [F(1), -2 * a, a * a + b * b]); d += 2   # roots a +- ib
            if a > 0: rhp += 2
        else:
            a = F(random.choice([-1, 1]) * random.randint(1, 40), random.randint(1, 30))
            p = mulp(p, [F(1), -a]); d += 1
            if a > 0: rhp += 1
    k = F(random.randint(1, 9), random.randint(1, 9))
    p = [k * x for x in p]
    r1, _ = routh_rhp(p); r2, _ = cauchy_rhp(p)
    r3, _ = sympy_rhp(p) if trial % 10 == 0 else (rhp, None)
    ntest += 1
    for r in (r1, r2, r3):
        if r is not None and r != rhp:
            nfail += 1; print('FAIL', p, rhp, r1, r2, r3); break
# imaginary-axis roots must be undecided
for p in ([F(1), F(0), F(1)], mulp([F(1), F(0), F(4)], [F(1), F(3), F(2)]), mulp([F(1), F(0)], [F(1), F(2), F(5)])):
    r1, n1 = routh_rhp(p); r2, n2 = cauchy_rhp(p); r3, n3 = sympy_rhp(p)
    print('imag-axis case', [str(x) for x in p], r1, r2, r3, n2, n3)
    assert r2 is None and r3 is None
    if r1 is not None: print('  NOTE Routh decided', r1)
# charpoly check against a matrix with known spectrum
M = [[F(2), F(1)], [F(0), F(-3)]]
assert charpoly(M) == [1, 1, -6]
print('tests', ntest, 'fail', nfail)
# repeated roots off the axis
for p, rhp in ((mulp(mulp([F(1), F(1)], [F(1), F(1)]), mulp([F(1), F(1), F(1)], [F(1), F(1), F(1)])), 0),
               (mulp(mulp([F(1), F(-1)], [F(1), F(-1)]), [F(1), F(2)]), 2),
               (mulp(mulp([F(1), F(-1), F(3)], [F(1), F(-1), F(3)]), [F(1), F(1, 3)]), 4)):
    r = (routh_rhp(p)[0], cauchy_rhp(p)[0], sympy_rhp(p)[0])
    print('repeated-root case expect', rhp, 'got', r)
    assert all(v is None or v == rhp for v in r) and r[1] == rhp and r[2] == rhp
print('repeated-root tests OK')
