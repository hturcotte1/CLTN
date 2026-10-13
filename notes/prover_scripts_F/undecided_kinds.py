"""Which quantity vanishes first for the undecided n=6 full-support cases at (eps,delta)=(1/2,3/2),
i.e. (beta,t)=(4,1); and the value set of 3x3 principal minors of K = I + J + 4N there.
Same 3000 random graphs (seed 7) as undecided_rate.py / undecided_factor.py. Exact arithmetic."""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
from collections import Counter
from itertools import product
import ctln
rng = random.Random(7)
n = 6
def first_zero(A, eps, delta):
    W = ctln.weight_matrix(A, eps, delta)
    M = [[(F(-1) if i == j else F(0)) + W[i][j] for j in range(n)] for i in range(n)]
    a = ctln.charpoly(M); m = len(a) - 1
    for i, c in enumerate(a):
        if c == 0: return ("coef", i)
    def aidx(k): return a[k] if 0 <= k <= m else F(0)
    H = [[aidx(2 * (i + 1) - (j + 1)) for j in range(m)] for i in range(m)]
    for k in range(1, m + 1):
        if ctln.det_exact([row[:k] for row in H[:k]]) == 0: return ("minor", k)
cnt = Counter()
for _ in range(3000):
    A = [[0 if i == j else rng.randint(0, 1) for j in range(n)] for i in range(n)]
    r = ctln.analyze_support(A, list(range(n)), F(1, 2), F(3, 2))
    if r["stable"] is None:
        cnt[(first_zero(A, F(1, 2), F(3, 2)), "positive=" + str(r["positive"]))] += 1
print("first vanishing quantity (index into [1,c1,..,c6]) x positivity:", dict(cnt))
vals = Counter()
for bits in product([0, 1], repeat=6):
    N = [[0] * 3 for _ in range(3)]; it = iter(bits)
    for i in range(3):
        for j in range(3):
            if i != j: N[i][j] = next(it)
    K = [[F((1 if i == j else 0) + 1 + 4 * N[i][j]) for j in range(3)] for i in range(3)]
    vals[ctln.det_exact(K)] += 1
print("3x3 principal minors of K at (beta,t)=(4,1):", sorted((int(k), v) for k, v in vals.items()))
