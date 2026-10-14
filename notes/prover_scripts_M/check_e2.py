"""EXACT checks of the e_2 closed form (sum of 2x2 principal minors of K = I + beta N + t J):
   e_2(K) = B(1+2t) - Nn(beta-1)(beta+1+2t) + O(1 - t(beta-2))
          = t*Phi(beta) + B + O - Nn(beta^2-1),   Phi(beta) = m(m-1) - beta*e,  e = #ordered nonedges = 2Nn + O
and of the consequence 'legal-stable somewhere => B > Nn' on all 95 n=5 legal-stable graphs."""
import sys, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from curve1d import curve_test
from run_scan_small import iso
def counts(A):
    m = len(A); B = Nn = O = 0
    for i in range(m):
        for j in range(i + 1, m):
            s = A[i][j] + A[j][i]
            if s == 2: B += 1
            elif s == 0: Nn += 1
            else: O += 1
    return B, Nn, O
random.seed(3)
for trial in range(300):
    m = random.randint(2, 7)
    A = [[0 if i == j else random.randint(0, 1) for j in range(m)] for i in range(m)]
    beta = F(random.randint(1, 60), random.randint(1, 20)); t = F(random.randint(1, 90), random.randint(1, 20))
    N = [[1 if (i != j and A[j][i] == 0) else 0 for j in range(m)] for i in range(m)]
    K = [[F(int(i == j)) + beta * N[i][j] + t for j in range(m)] for i in range(m)]
    e2 = sum(K[i][i] * K[j][j] - K[i][j] * K[j][i] for i in range(m) for j in range(i + 1, m))
    cp = ctln.charpoly([[-K[i][j] for j in range(m)] for i in range(m)])     # det(mu I + K), high->low
    B, Nn, O = counts(A); e = 2 * Nn + O
    assert e2 == cp[2] == B * (1 + 2 * t) - Nn * (beta - 1) * (beta + 1 + 2 * t) + O * (1 - t * (beta - 2))
    assert e2 == t * (m * (m - 1) - beta * e) + B + O - Nn * (beta * beta - 1)
print('e_2 closed form: 300 random exact checks OK')
bad = 0; n = 0
for A in iso(5):
    A = A.astype(int).tolist()
    if curve_test(A)['verdict'] == 'stable_somewhere':
        n += 1; B, Nn, O = counts(A)
        if not B > Nn: bad += 1
print('n=5 legal-stable graphs:', n, ' violating B > Nn:', bad)
