"""Sanity checks (NOT proofs) for Lemmas D and E. Exact Fractions throughout.
Enumerates ALL supports of random graphs (n = 3..6) having a dominated node."""
import random, sys
from fractions import Fraction as F
from itertools import combinations
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs

def W_of(A, eps, delta):
    n = len(A)
    return [[F(0) if i == j else (eps - 1 if A[j][i] else -1 - delta) for j in range(n)] for i in range(n)]

def fp_on_support(A, sigma, eps, delta, theta, strict):
    """Return (x or None): x solves (I-W_s)x = theta_s, x>0, and off-support inputs <0 (strict) or <=0."""
    W = W_of(A, eps, delta)
    m = len(sigma)
    IW = [[(F(1) if a == b else F(0)) - W[sigma[a]][sigma[b]] for b in range(m)] for a in range(m)]
    x = ctln.solve_exact(IW, [theta[s] for s in sigma])
    if x is None or not all(v > 0 for v in x):
        return None, None
    for k in range(len(A)):
        if k in sigma: continue
        y = theta[k] + sum(W[k][sigma[b]] * x[b] for b in range(m))
        if (y >= 0) if strict else (y > 0):
            return None, None
    M = [[(F(-1) if a == b else F(0)) + W[sigma[a]][sigma[b]] for b in range(m)] for a in range(m)]
    st, _ = ctln.hurwitz_stable(ctln.charpoly(M))
    return x, st

def all_supports(n):
    for r in range(1, n + 1):
        yield from combinations(range(n), r)

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
nD = nE = nG = 0
singular_seen = 0
while nG < int(sys.argv[2]) if len(sys.argv) > 2 else 300:
    n = random.randint(3, 6)
    A = [[0 if i == j else int(random.random() < random.choice([0.3, 0.5, 0.7])) for j in range(n)] for i in range(n)]
    pairs = [(k, j) for j in range(n) for k in range(n) if k != j and graphs.graphically_dominates(A, k, j)]
    if not pairs: continue
    nG += 1
    k, j = random.choice(pairs)
    for _ in range(3):
        eps = random.choice([F(random.randint(1, 99), 100), F(1), F(random.randint(101, 300), 100)])  # includes eps >= 1
        delta = F(random.randint(1, 400), 100)
        th = F(random.randint(1, 9), random.randint(1, 9))
        theta = [th] * n
        # Lemma D: no fixed point with x_j > 0 (weak off-support), any eps > 0
        for sig in all_supports(n):
            if j in sig:
                x, _ = fp_on_support(A, list(sig), eps, delta, theta, strict=False)
                assert x is None, ("Lemma D violated", A, k, j, sig, eps, delta)
                nD += 1
        # Lemma E: stable supports (strict and weak) agree between G and G-j
        keep = [v for v in range(n) if v != j]
        Ap = [[A[a][b] for b in keep] for a in keep]
        for strict in (True, False):
            for sig in all_supports(n - 1):
                sigG = [keep[v] for v in sig]
                xG, sG = fp_on_support(A, sigG, eps, delta, theta, strict)
                xP, sP = fp_on_support(Ap, list(sig), eps, delta, [th] * (n - 1), strict)
                assert (xG is None) == (xP is None), ("Lemma E FP mismatch", A, k, j, sigG, strict)
                if xG is not None:
                    assert xG == xP and sG == sP
                nE += 1
print(f"graphs with a dominated node: {nG}; Lemma D support checks: {nD}; Lemma E support checks: {nE}; all pass")

# Non-uniform theta: the 2-node example j -> k, k -/-> j, theta_j > theta_k/(1-eps)
A = [[0, 1], [0, 0]]       # node 0 = j sends to node 1 = k
assert graphs.graphically_dominates(A, 1, 0)
eps, delta = F(1, 4), F(1, 2)
x, st = fp_on_support(A, [0], eps, delta, [F(10), F(1)], strict=True)
print("non-uniform theta (theta_j=10, theta_k=1): support {j} fixed point x =", x, "stable:", st)
