from fractions import Fraction as F
def berkowitz(M):
    """charpoly det(zI - M), hi->lo, exact (ints or Fractions)."""
    n = len(M)
    C = [1]
    for k in range(n):        # add row/col k
        # leading principal submatrix of size k, column R = M[0..k-1][k], row S = M[k][0..k-1], a = M[k][k]
        a = M[k][k]
        R = [M[i][k] for i in range(k)]; S = [M[k][j] for j in range(k)]
        # Toeplitz column: 1, -a, -S R, -S A R, -S A^2 R, ...
        col = [1, -a]
        v = R[:]
        for p in range(k):
            col.append(-sum(S[j] * v[j] for j in range(k)))
            v = [sum(M[i][j] * v[j] for j in range(k)) for i in range(k)]
        # new C = T * C, T lower-triangular Toeplitz (k+2) x (k+1)
        newC = []
        for i in range(k + 2):
            newC.append(sum(col[i - j] * C[j] for j in range(max(0, i - (k + 1)), min(i, k) + 1) if 0 <= i - j < len(col)))
        C = newC
    return C

def polysub(a, b):
    L = max(len(a), len(b)); a = [0] * (L - len(a)) + a; b = [0] * (L - len(b)) + b
    r = [x - y for x, y in zip(a, b)]
    while r and r[0] == 0: r = r[1:]
    return r

def psi_A(A):
    n = len(A)
    AmJ = [[A[i][j] - 1 for j in range(n)] for i in range(n)]
    return polysub(berkowitz(AmJ), berkowitz(A))   # det(zI-A+J) - det(zI-A) = 1^T adj(zI-A) 1

def shift(p, c):
    """coefficients of p(x + c), hi->lo."""
    p = [F(v) for v in p]; n = len(p) - 1
    out = p[:]
    for i in range(n):
        for j in range(1, n - i + 1):
            out[j] += c * out[j - 1]
    return out

def hurwitz_all_positive(p):
    """True iff leading coeff > 0 and all leading Hurwitz minors > 0 (exact). False/None otherwise (None: zero minor)."""
    p = [F(v) for v in p]
    if p[0] < 0: p = [-v for v in p]
    d = len(p) - 1
    if d == 0: return True
    a = lambda i: p[i] if 0 <= i <= d else F(0)
    H = [[a(2 * (i + 1) - (j + 1)) for j in range(d)] for i in range(d)]
    for k in range(d):
        piv = H[k][k]
        if piv == 0: return None
        if piv < 0: return False     # leading minors' ratio negative -> some minor negative (product of pivots)
        for i in range(k + 1, d):
            f = H[i][k] / piv
            if f:
                for j in range(k, d): H[i][j] -= f * H[k][j]
    return True

def polydiv(a, b):
    a = [F(v) for v in a]; q = []
    while len(a) >= len(b):
        c = a[0] / b[0]; q.append(c)
        a = [x - c * y for x, y in zip(a, b + [0] * (len(a) - len(b)))][1:]
    while a and a[0] == 0: a = a[1:]
    return q, a

