"""Independent exact verifier for n = 8 (pure fractions.Fraction; no sympy/numpy; no floats in decisions).
Implements CLAUDE.md section 2 from scratch; adapted from notes/verify_124/v_fraction/vf.py (n = 7).
Does not import any project module."""
from fractions import Fraction as F

# ---------- digraph6 ----------
def decode_d6(s):
    s = s.strip()
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    need = n * n
    assert len(bits) >= need and len(s[2:]) == (need + 5) // 6, s
    assert all(b == 0 for b in bits[need:]), "nonzero padding"
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def encode_d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(63 + n) + ''.join(chr(63 + int(''.join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))

# ---------- graph notions ----------
def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def has_loops(A):
    return any(A[i][i] for i in range(len(A)))

# ---------- CTLN matrices ----------
def W_matrix(A, eps, delta):
    n = len(A)
    W = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: W[i][j] = F(0)
            elif A[j][i] == 1: W[i][j] = -1 + eps      # j -> i : weak inhibition
            else: W[i][j] = -1 - delta                 # strong inhibition
    return W

def q_of(eps, delta):
    return delta * (1 - eps) / eps

def legal(eps, delta):
    return 0 < eps < 1 and delta > 0 and q_of(eps, delta) > 1

# ---------- linear algebra ----------
def solve(M, b):
    """Gaussian elimination with exact pivot search. Returns (x, det) or (None, 0) if singular."""
    n = len(M)
    a = [list(map(F, M[i])) + [F(b[i])] for i in range(n)]
    det = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None: return None, F(0)
        if p != c: a[c], a[p] = a[p], a[c]; det = -det
        det *= a[c][c]
        for r in range(c + 1, n):
            if a[r][c] != 0:
                f = a[r][c] / a[c][c]
                for k in range(c, n + 1): a[r][k] -= f * a[c][k]
    x = [F(0)] * n
    for r in range(n - 1, -1, -1):
        s = a[r][n] - sum(a[r][k] * x[k] for k in range(r + 1, n))
        x[r] = s / a[r][r]
    for i in range(n):  # exact residual check
        assert sum(F(M[i][j]) * x[j] for j in range(n)) == b[i], "residual nonzero"
    return x, det

def det(M):
    n = len(M)
    if n == 0: return F(1)
    a = [list(map(F, r)) for r in M]; d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: a[c], a[p] = a[p], a[c]; d = -d
        d *= a[c][c]
        for r in range(c + 1, n):
            if a[r][c] != 0:
                f = a[r][c] / a[c][c]
                for k in range(c, n): a[r][k] -= f * a[c][k]
    return d

def matmul(X, Y):
    n, m, p = len(X), len(Y), len(Y[0])
    return [[sum(X[i][k] * Y[k][j] for k in range(m)) for j in range(p)] for i in range(n)]

def charpoly_faddeev(B):
    """Coefficients [1, c1, ..., cn] of det(lambda I - B) by Faddeev-LeVerrier; Cayley-Hamilton check enforced."""
    n = len(B)
    c = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = [[Mk[i][j] + (c[-1] if i == j else 0) for j in range(n)] for i in range(n)]  # M_k = B M_{k-1} + c_{k-1} I
        BM = matmul(B, Mk)
        ck = -sum(BM[i][i] for i in range(n)) / k
        c.append(ck)
        Mk = BM
    assert all(Mk[i][j] + (c[-1] if i == j else 0) == 0 for i in range(n) for j in range(n)), "CH check failed"
    return c

def charpoly_berkowitz(B):
    """Berkowitz algorithm (division free). Returns [1, c1, ..., cn] of det(lambda I - B)."""
    n = len(B)
    vect = [F(1), -B[0][0]]
    for r in range(1, n):
        Rr = [B[r][j] for j in range(r)]
        Cc = [B[i][r] for i in range(r)]
        Ar = [row[:r] for row in B[:r]]
        a = B[r][r]
        col = [F(1), -a]
        v = Cc[:]
        for k in range(r):
            col.append(-sum(Rr[i] * v[i] for i in range(r)))
            v = [sum(Ar[i][j] * v[j] for j in range(r)) for i in range(r)]
        new = []
        for i in range(r + 2):
            s = F(0)
            for j in range(r + 1):
                if 0 <= i - j < len(col): s += col[i - j] * vect[j]
            new.append(s)
        vect = new
    return vect

def hurwitz_minors(c):
    """c = [1, a1, ..., an] monic. Leading principal minors Delta_1..Delta_n of the Hurwitz matrix
    H[i][j] = a_{2j - i} (1-indexed), a_0 = 1, a_k = 0 outside 0..n."""
    n = len(c) - 1
    def a(k): return c[k] if 0 <= k <= n else F(0)
    H = [[a(2 * (j + 1) - (i + 1)) for j in range(n)] for i in range(n)]
    return [det([row[:k] for row in H[:k]]) for k in range(1, n + 1)]

def routh_first_column(c):
    n = len(c) - 1
    r0 = [c[k] for k in range(0, n + 1, 2)]
    r1 = [c[k] for k in range(1, n + 1, 2)]
    w = len(r0)
    r1 = r1 + [F(0)] * (w - len(r1))
    rows = [r0, r1]
    for _ in range(n - 1):
        p, q = rows[-2], rows[-1]
        if q[0] == 0: return None
        new = [(q[0] * p[k + 1] - p[0] * q[k + 1]) / q[0] if k + 1 < w else F(0) for k in range(w)]
        rows.append(new)
    return [r[0] for r in rows[: n + 1]]

def stability(c):
    """Three-valued: 'stable' iff every Delta_k > 0; 'unstable' iff some Delta_k < 0;
    otherwise 'undecided' (a vanishing minor and no negative one)."""
    D = hurwitz_minors(c)
    if all(d > 0 for d in D): v = 'stable'
    elif any(d < 0 for d in D): v = 'unstable'
    else: v = 'undecided'
    return v, D

def analyze(A, eps, delta):
    n = len(A)
    W = W_matrix(A, eps, delta)
    M = [[(F(1) if i == j else F(0)) - W[i][j] for j in range(n)] for i in range(n)]  # I - W
    x, dM = solve(M, [F(1)] * n)
    positive = x is not None and all(xi > 0 for xi in x)
    B = [[W[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)]  # -I + W
    c1 = charpoly_faddeev(B)
    c2 = charpoly_berkowitz(B)
    assert c1 == c2, "charpoly methods disagree"
    assert c1[-1] == det(M), "constant term != det(I-W)"
    assert c1[1] == n, "trace check"
    v, D = stability(c1)
    rc = routh_first_column(c1)
    if rc is not None and all(d != 0 for d in D):
        Dp = [F(1)] + D
        assert rc == [F(1)] + [Dp[k] / Dp[k - 1] for k in range(1, n + 1)], "Routh/Hurwitz mismatch"
        routh = 'stable' if all(r > 0 for r in rc) else 'unstable'
        assert (routh == 'stable') == (v == 'stable')
    return dict(legal=legal(eps, delta), x=x, positive=positive, charpoly=c1, minors=D, stability=v,
                complete=is_complete_bidirected(A))

# ---------- exact root isolation helpers ----------
def peval(coeffs_hi_to_lo, x):
    r = F(0)
    for c in coeffs_hi_to_lo:
        r = r * x + c
    return r

def sgn(v): return (v > 0) - (v < 0)

def refine_root(coeffs, a, b, width):
    """Given a sign-change isolating interval [a,b] of a root of the integer polynomial `coeffs`,
    bisect exactly until b - a <= width. Returns (a, b, note). If an endpoint is an exact root,
    returns that point twice."""
    a, b = F(a), F(b)
    pa, pb = peval(coeffs, a), peval(coeffs, b)
    if pa == 0: return a, a, 'endpoint_a_is_root'
    if pb == 0: return b, b, 'endpoint_b_is_root'
    if sgn(pa) == sgn(pb): return None, None, 'NO_SIGN_CHANGE'
    while b - a > width:
        m = (a + b) / 2
        pm = peval(coeffs, m)
        if pm == 0: return m, m, 'midpoint_is_root'
        if sgn(pm) == sgn(pa): a, pa = m, pm
        else: b, pb = m, pm
    return a, b, 'bisected'

def simplest_between(a, b):
    """Simplest rational (smallest denominator, then numerator) strictly inside the open interval (a,b),
    0 < a < b. b may be None (= +infinity)."""
    assert b is None or a < b
    fl = a.numerator // a.denominator
    c = F(fl + 1)
    if b is None or c < b: return c
    # a, b both in [fl, fl+1)
    ra = a - fl
    rb = b - fl
    if ra == 0:
        inner = simplest_between(1 / rb, None)
    else:
        inner = simplest_between(1 / rb, 1 / ra)
    return fl + 1 / inner
