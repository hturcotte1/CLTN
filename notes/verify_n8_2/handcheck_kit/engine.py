"""handcheck_kit engine -- own exact implementation written from the definitions
(CLAUDE.md section 2).  Fractions only.  Three-valued stability.

A[i][j] == 1 iff arc i -> j.
W[i][j] = 0 (i==j); -1+eps if A[j][i]==1; -1-delta otherwise.
I - W = -delta*I + (1+delta)*J - (eps+delta)*A^T   (checked entrywise below).
Full support stable iff (I-W)x = 1 has x > 0 and all eigenvalues of -I+W have Re < 0,
decided by Hurwitz minors of charpoly(-I+W) (all positive => stable; some negative => unstable;
zero before any negative => UNDECIDED).
"""
from fractions import Fraction as Fr


def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    assert len(bits) >= n * n
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]


def d6_encode(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    out = ['&', chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return ''.join(out)


def W_matrix(A, eps, delta):
    n = len(A)
    return [[Fr(0) if i == j else (Fr(-1) + eps if A[j][i] == 1 else Fr(-1) - delta)
             for j in range(n)] for i in range(n)]


def I_minus_W(A, eps, delta):
    n = len(A)
    W = W_matrix(A, eps, delta)
    M = [[(Fr(1) if i == j else Fr(0)) - W[i][j] for j in range(n)] for i in range(n)]
    # cross-check closed form
    for i in range(n):
        for j in range(n):
            cf = (-delta if i == j else Fr(0)) + (1 + delta) - (eps + delta) * A[j][i]
            assert cf == M[i][j], ("closed form mismatch", i, j)
    return M


def solve(M, b):
    """Gaussian elimination over Fractions; returns x or None if singular."""
    n = len(M)
    a = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None:
            return None
        a[c], a[p] = a[p], a[c]
        pv = a[c][c]
        a[c] = [v / pv for v in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [a[r][k] - f * a[c][k] for k in range(n + 1)]
    return [a[i][n] for i in range(n)]


def det(M):
    n = len(M)
    a = [row[:] for row in M]
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            d = -d
        d *= a[c][c]
        for r in range(c + 1, n):
            if a[r][c] != 0:
                f = a[r][c] / a[c][c]
                a[r] = [a[r][k] - f * a[c][k] for k in range(n)]
    return d


def matmul(X, Y):
    n = len(X)
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def charpoly(M):
    """Faddeev-LeVerrier: coefficients c[0..n] of det(s I - M) = s^n + c1 s^{n-1} + ... + cn, c0 = 1."""
    n = len(M)
    I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    c = [Fr(1)]
    Mk = [[Fr(0)] * n for _ in range(n)]  # M_0 = 0
    for k in range(1, n + 1):
        # M_k = M * M_{k-1} + c_{k-1} I
        Mk = matmul(M, Mk)
        for i in range(n):
            Mk[i][i] += c[k - 1]
        AMk = matmul(M, Mk)
        tr = sum(AMk[i][i] for i in range(n))
        c.append(-tr / k)
    return c


def charpoly_check(M, c):
    """Independent check: evaluate det(sI - M) at s = 0..n and compare with polynomial values."""
    n = len(M)
    for s in range(n + 1):
        N = [[(Fr(s) if i == j else Fr(0)) - M[i][j] for j in range(n)] for i in range(n)]
        val = sum(c[k] * Fr(s) ** (n - k) for k in range(n + 1))
        assert det(N) == val, ("charpoly check failed", s)
    return True


def hurwitz_minors(c):
    """c = [1, a1, ..., an]. Hurwitz matrix H (n x n): H[i][j] = a_{2j - i} (1-indexed), a_0 = 1, out of range 0."""
    n = len(c) - 1
    def a(k):
        return c[k] if 0 <= k <= n else Fr(0)
    H = [[a(2 * (j + 1) - (i + 1)) for j in range(n)] for i in range(n)]
    minors = [det([row[:k] for row in H[:k]]) for k in range(1, n + 1)]
    return H, minors


def classify(minors):
    """Three-valued Routh-Hurwitz: all > 0 -> 'stable'; first non-positive minor is < 0 -> 'unstable';
    first non-positive minor == 0 -> 'undecided'."""
    for m in minors:
        if m == 0:
            return 'undecided'
        if m < 0:
            return 'unstable'
    return 'stable'


def analyze(A, eps, delta):
    n = len(A)
    q = delta * (1 - eps) / eps
    M = I_minus_W(A, eps, delta)
    x = solve(M, [Fr(1)] * n)
    W = W_matrix(A, eps, delta)
    J = [[W[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)]  # -I + W
    c = charpoly(J)
    charpoly_check(J, c)
    H, minors = hurwitz_minors(c)
    stab = classify(minors)
    positive = x is not None and all(v > 0 for v in x)
    return dict(q=q, M=M, x=x, positive=positive, charpoly=c, H=H, minors=minors, stability=stab,
                legal=(eps > 0 and delta > 0 and q > 1),
                complete=all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j))


def missing_pair(A):
    n = len(A)
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] == 0:
                return (i, j)
    return None
