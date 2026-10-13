"""Probe: at eps=1/2, delta=3/2 (q=3/2), random n=6 graphs, full support.
For each UNDECIDED case (hurwitz_stable -> None), find the first vanishing quantity
(char-poly coefficient or Hurwitz minor) and test whether it vanishes at two other
random rational legal parameter points. If it vanishes identically, Lemma O says the
graph is NEVER stable (decided), not undecided.  Exact arithmetic only."""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
import ctln
rng = random.Random(7)

def minors(A, eps, delta):
    n = len(A)
    W = ctln.weight_matrix(A, eps, delta)
    M = [[(F(-1) if i == j else F(0)) + W[i][j] for j in range(n)] for i in range(n)]
    a = ctln.charpoly(M)
    m = len(a) - 1
    def aidx(k): return a[k] if 0 <= k <= m else F(0)
    H = [[aidx(2 * (i + 1) - (j + 1)) for j in range(m)] for i in range(m)]
    return a, [ctln.det_exact([row[:k] for row in H[:k]]) for k in range(1, m + 1)]

def first_zero(A, eps, delta):
    a, h = minors(A, eps, delta)
    for i, c in enumerate(a):
        if c == 0: return ("coef", i)
    for k, v in enumerate(h, 1):
        if v == 0: return ("minor", k)
    return None

def val(A, eps, delta, tag):
    a, h = minors(A, eps, delta)
    return a[tag[1]] if tag[0] == "coef" else h[tag[1] - 1]

n = 6; tot = und = ident = 0; examples = []
probes = [(F(1, 7), F(5, 3)), (F(2, 9), F(13, 4)), (F(1, 30), F(1, 10))]
for _ in range(3000):
    A = [[0 if i == j else rng.randint(0, 1) for j in range(n)] for i in range(n)]
    tot += 1
    r = ctln.analyze_support(A, list(range(n)), F(1, 2), F(3, 2))
    if r["stable"] is not None:
        continue
    und += 1
    tag = first_zero(A, F(1, 2), F(3, 2))
    if tag is None:
        continue
    if all(val(A, e, d, tag) == 0 for e, d in probes):
        ident += 1
        if len(examples) < 3: examples.append((A, tag, r["positive"]))
print("graphs", tot, "undecided", und, "first-vanishing quantity also zero at 3 other legal points", ident)
for e in examples: print(e)
