"""For a few undecided n=6 full-support cases at (eps,delta)=(1/2,3/2) i.e. (beta,t)=(4,1),
factor the vanishing quantity as a polynomial in (beta,t)."""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
import sympy as sp
import ctln
rng = random.Random(7)
b, t, mu = sp.symbols("beta t mu")
n = 6; shown = 0
from collections import Counter
facs = Counter()
for _ in range(3000):
    A = [[0 if i == j else rng.randint(0, 1) for j in range(n)] for i in range(n)]
    r = ctln.analyze_support(A, list(range(n)), F(1, 2), F(3, 2))
    if r["stable"] is not None:
        continue
    N = sp.Matrix(n, n, lambda i, j: 0 if i == j else 1 - A[j][i])
    K = sp.eye(n) + b * N + t * sp.ones(n, n)
    P = sp.Poly((mu * sp.eye(n) + K).det(method="berkowitz"), mu)
    a = P.all_coeffs()
    m = n
    def aidx(k): return a[k] if 0 <= k <= m else 0
    H = sp.Matrix(m, m, lambda i, j: aidx(2 * (i + 1) - (j + 1)))
    found = None
    for i, c in enumerate(a):
        if sp.expand(c.subs({b: 4, t: 1})) == 0:
            found = ("coef", i, c); break
    if not found:
        for k in range(1, m + 1):
            v = sp.expand(H[:k, :k].det(method="berkowitz"))
            if v.subs({b: 4, t: 1}) == 0:
                found = ("minor", k, v); break
    f = sp.factor_list(found[2])
    vanishing = [fa for fa, e in f[1] if fa.subs({b: 4, t: 1}) == 0]
    facs[str(vanishing)] += 1
    if shown < 4:
        print(found[0], found[1], "positive:", r["positive"], "vanishing factor(s):", vanishing, " detail:", r["stability_detail"])
        shown += 1
print(facs.most_common(10))
