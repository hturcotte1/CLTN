"""never_side engine -- OWN exact decision code written from CLAUDE.md section 2 / the task definitions.

digraph6 -> A (A[i][j] = 1 iff arc i -> j).
CTLN at eps = 1/(t+1), delta = (beta - 1) eps  (beta = bn/bd rational > 2 + 1/t for legality):
    W[i][j] = 0 (i == j); -1 + eps if A[j][i] == 1; -1 - delta otherwise.
    M = -I + W.  Full support is a stable fixed point iff
      (1) (I - W) x = 1 has a solution with every x_i > 0, and
      (3) every eigenvalue of M has real part < 0.
Everything is exact: integers and fractions.Fraction only.
Stability: Hurwitz minors of the characteristic polynomial of L = S*M, S = bd*(t+1) > 0 (integer matrix with
the same stability as M).  Three-valued: True / False / None (a vanishing minor or coefficient).
"""
from fractions import Fraction as F

def parse_d6(s):
    s = s.strip()
    assert s[0] == '&', "not digraph6"
    n = ord(s[1]) - 63
    assert 1 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, "bad digraph6 char"
        bits += [(v >> (5 - k)) & 1 for k in range(6)]
    assert len(bits) >= n * n and not any(bits[n * n:]), "digraph6 padding not zero"
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def ones_arcs(A):
    return sum(map(sum, A))

# ---------------------------------------------------------------- integer linear algebra

def charpoly_int(M):
    """Faddeev-LeVerrier.  Returns [1, c1, ..., cn]: det(lam I - M) = lam^n + c1 lam^(n-1) + ... + cn."""
    n = len(M)
    c = [1]
    Mk = [row[:] for row in M]
    for k in range(1, n + 1):
        if k > 1:
            for i in range(n):
                Mk[i][i] += c[-1]
            Mk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(Mk[i][i] for i in range(n))
        assert tr % k == 0, "Faddeev-LeVerrier division not exact"
        c.append(-tr // k)
    return c

def det_int(M):
    """Bareiss fraction-free determinant of a square integer matrix (copied)."""
    k = len(M)
    if k == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(k - 1):
        if M[c][c] == 0:
            sw = next((r for r in range(c + 1, k) if M[r][c] != 0), None)
            if sw is None:
                return 0
            M[c], M[sw] = M[sw], M[c]
            sign = -sign
        pc = M[c][c]
        for r in range(c + 1, k):
            mrc = M[r][c]
            row = M[r]
            crow = M[c]
            for j in range(c + 1, k):
                row[j] = (row[j] * pc - mrc * crow[j]) // prev
        prev = pc
    return sign * M[k - 1][k - 1]

def hurwitz(coeffs):
    """coeffs = [a0, a1, ..., ad] of a0 x^d + ... + ad, integers.  Returns (verdict, detail):
       True  = all roots have real part < 0 (all coefficients > 0 and all leading principal Hurwitz minors > 0)
       False = some coefficient < 0 or some Hurwitz minor < 0 with no vanishing one encountered first
       None  = a coefficient or a minor vanishes (boundary / undecided).
    a0 must be > 0."""
    d = len(coeffs) - 1
    a = coeffs
    assert a[0] > 0
    if any(c == 0 for c in a):
        return None, "zero coefficient"
    if any(c < 0 for c in a):
        return False, "negative coefficient"
    def ent(i, j):          # 1-indexed Hurwitz matrix entry a_{2i-j}
        idx = 2 * i - j
        return a[idx] if 0 <= idx <= d else 0
    for k in range(1, d + 1):
        H = [[ent(i, j) for j in range(1, k + 1)] for i in range(1, k + 1)]
        m = det_int(H)
        if m == 0:
            return None, f"Hurwitz minor {k} vanishes"
        if m < 0:
            return False, f"Hurwitz minor {k} negative"
    return True, "all Hurwitz minors positive"

def solve_fraction(M, rhs):
    """Exact Gauss-Jordan over Fractions.  Returns the solution list or None if M is singular."""
    n = len(M)
    aug = [[F(M[i][j]) for j in range(n)] + [F(rhs[i])] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if aug[r][c] != 0), None)
        if p is None:
            return None
        aug[c], aug[p] = aug[p], aug[c]
        pv = aug[c][c]
        aug[c] = [v / pv for v in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [vr - f * vc for vr, vc in zip(aug[r], aug[c])]
    return [aug[i][n] for i in range(n)]

# ---------------------------------------------------------------- the CTLN point decision

def scaled_M(A, t, bn, bd):
    """L = S * (-I + W) with S = bd*(t+1), as an integer matrix (t integer, beta = bn/bd)."""
    n = len(A)
    S = bd * (t + 1)
    weak = -bd * t                 # S*(-1 + eps)
    strong = -(bd * t + bn)        # S*(-1 - delta)
    return [[(-S if i == j else (weak if A[j][i] == 1 else strong)) for j in range(n)] for i in range(n)], S

def decide_point(A, t, bn, bd):
    """Exact verdict at eps = 1/(t+1), delta = (beta-1) eps, beta = bn/bd.
    Returns dict(legal, stable (True/False/None), stable_detail, positive (True/False/None), x (Fractions))."""
    n = len(A)
    t = int(t); bn = int(bn); bd = int(bd)
    beta = F(bn, bd)
    eps = F(1, t + 1); delta = (beta - 1) * eps
    q = delta * (1 - eps) / eps
    legal = eps > 0 and delta > 0 and q > 1
    L, S = scaled_M(A, t, bn, bd)
    # stability of M = L/S: characteristic polynomial of L (monic) must be Hurwitz
    cp = charpoly_int(L)
    st, det = hurwitz(cp)
    # positivity: (I - W) x = 1  <=>  (-L) x = S * 1
    x = solve_fraction([[-v for v in row] for row in L], [S] * n)
    if x is None:
        pos = None
    else:
        pos = all(v > 0 for v in x)
    return {"legal": legal, "q": q, "eps": eps, "delta": delta, "stable": st, "stable_detail": det,
            "positive": pos, "x": x, "charpoly": cp}

def W_matrix(A, eps, delta):
    """the W matrix itself (Fractions), for cross-checks against the definition"""
    n = len(A)
    return [[F(0) if i == j else (-1 + eps if A[j][i] == 1 else -1 - delta) for j in range(n)] for i in range(n)]

def decide_point_direct(A, eps, delta):
    """Slow reference: build W literally, M = -I + W over Fractions, char poly by Fraction Faddeev-LeVerrier
    (scaled to integers by the lcm), Hurwitz as above; positivity by solving (I - W) x = 1."""
    from math import lcm
    n = len(A)
    W = W_matrix(A, eps, delta)
    M = [[W[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)]
    den = 1
    for row in M:
        for v in row:
            den = lcm(den, v.denominator)
    Lm = [[int(v * den) for v in row] for row in M]
    cp = charpoly_int(Lm)
    st, det = hurwitz(cp)
    IW = [[(1 if i == j else 0) - W[i][j] for j in range(n)] for i in range(n)]
    x = solve_fraction(IW, [1] * n)
    pos = None if x is None else all(v > 0 for v in x)
    return {"stable": st, "stable_detail": det, "positive": pos, "x": x}

# ---------------------------------------------------------------- graph-theoretic tests (own code)

def graphically_dominates(A, k, j):
    """Curto 2510.05098 Def 1.1: k dominates j iff (ii) j -> k and k -/-> j, and (i) every i not in {j,k} with
    i -> j also has i -> k."""
    n = len(A)
    if j == k or not (A[j][k] == 1 and A[k][j] == 0):
        return False
    return all(not (A[i][j] == 1 and A[i][k] == 0) for i in range(n) if i != j and i != k)

def is_domination_free(A):
    n = len(A)
    return not any(graphically_dominates(A, k, j) for j in range(n) for k in range(n))

def is_complete(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)
