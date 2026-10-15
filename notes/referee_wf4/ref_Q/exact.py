"""Small exact helpers written for this referee pass (independent of the theory scripts).
Integer Bareiss determinants, Lagrange interpolation over Fractions, Routh array (three-valued)."""
from fractions import Fraction as F

def bareiss_det(M):
    M = [list(r) for r in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]; sign = -sign; break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def interp(xs, ys):
    """coefficients low->high of the interpolating polynomial (Fractions)."""
    n = len(xs); coef = [F(0)] * n
    for i in range(n):
        num = [F(1)]; den = F(1)
        for j in range(n):
            if j == i: continue
            num = [F(0)] + num
            for a in range(len(num) - 1): num[a] -= xs[j] * num[a + 1]
            den *= (xs[i] - xs[j])
        for a in range(n): coef[a] += num[a] * F(ys[i]) / den
    return coef

def charpoly_int(M):
    """det(z I - M) for an integer matrix, low->high coefficients, via n+1 Bareiss determinants."""
    n = len(M); xs = list(range(n + 1)); ys = []
    for z in xs:
        ys.append(bareiss_det([[(z if i == j else 0) - M[i][j] for j in range(n)] for i in range(n)]))
    c = interp([F(v) for v in xs], ys)
    assert all(v.denominator == 1 for v in c) and c[-1] == 1
    return [int(v) for v in c]

def psi_int(D):
    """1^T adj(zI - D) 1 = det(zI - D + J) - det(zI - D) (matrix determinant lemma), low->high."""
    n = len(D); xs = list(range(n + 1)); ys = []
    for z in xs:
        Mz = [[(z if i == j else 0) - D[i][j] for j in range(n)] for i in range(n)]
        MJ = [[v + 1 for v in r] for r in Mz]
        ys.append(bareiss_det(MJ) - bareiss_det(Mz))
    c = interp([F(v) for v in xs], ys)
    assert all(v.denominator == 1 for v in c) and c[-1] == 0
    return [int(v) for v in c[:-1]]

def poly_shift_neg(c, r):
    """given p(z) low->high, return coefficients (high->low) of p(-r - x) in x."""
    n = len(c) - 1
    # p(-r-x) = sum c_i (-r - x)^i
    out = [F(0)] * (n + 1)   # low->high in x
    from math import comb
    for i, ci in enumerate(c):
        if ci == 0: continue
        for j in range(i + 1):   # (-r-x)^i = sum C(i,j) (-x)^j (-r)^(i-j)
            out[j] += F(ci) * comb(i, j) * (-1) ** j * (F(-r)) ** (i - j)
    return out[::-1]

def routh(a):
    """a high->low, real exact. Returns True (all roots Re<0), False (some Re>0 or detected non-Hurwitz), None (zero pivot)."""
    a = [F(v) for v in a]
    while a and a[0] == 0: a = a[1:]
    if a[0] < 0: a = [-v for v in a]
    n = len(a) - 1
    if n == 0: return True
    r0 = a[0::2]; r1 = a[1::2]
    firstcol = [r0[0]]
    for _ in range(n):
        if len(r1) == 0 or r1[0] == 0:
            return None
        firstcol.append(r1[0])
        L = max(len(r0), len(r1) + 1)
        r0e = r0 + [F(0)] * (L - len(r0)); r1e = r1 + [F(0)] * (L - len(r1))
        nr = [(r1e[0] * r0e[i + 1] - r0e[0] * r1e[i + 1]) / r1e[0] for i in range(L - 1)]
        while nr and nr[-1] == 0 and len(nr) > 0: nr = nr[:-1]
        r0, r1 = r1, nr
        if len(firstcol) == n + 1: break
    if len(firstcol) < n + 1: return None
    if any(v == 0 for v in firstcol): return None
    return all(v > 0 for v in firstcol)
