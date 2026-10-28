"""Referee's own exact filter (pure Fractions): F1 domination-free; F2 alpha(A|perp) < -1/2 strictly via Hurwitz minors of
h(mu) = chi_perp(mu - 1/2), chi_perp(lam) = det(lam I - P A P)/lam with P = I - J/n (an independent route to A|perp:
filter8 uses chi_A - chi_{A+J}).  Flags: S strict, B (zero coefficient or zero minor: alpha >= -1/2, not strictly below),
N non-survivor (dominated, or a negative coefficient / negative minor: alpha > -1/2 ... precisely: h not Hurwitz)."""
import sys
from fractions import Fraction as Fr
from math import comb
def decode_d6(s):
    assert s[0] == '&'; n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - b)) & 1 for b in range(6)]
    return n, [[bits[i * n + j] for j in range(n)] for i in range(n)]
def dom_free(n, A):
    for j in range(n):
        for k in range(n):
            if j == k or not (A[j][k] == 1 and A[k][j] == 0): continue
            if all(A[i][k] == 1 for i in range(n) if i not in (j, k) and A[i][j] == 1): return False
    return True
def matmul(X, Y):
    n = len(X); return [[sum(X[i][l] * Y[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
def charpoly(M):
    """Faddeev-LeVerrier: coefficients of det(lam I - M), highest first, exact Fractions."""
    n = len(M); c = [Fr(1)]; Mk = [row[:] for row in M]
    for k in range(1, n + 1):
        if k > 1:
            Mk = matmul(M, Mk)
        tr = sum(Mk[i][i] for i in range(n))
        ck = -tr / k; c.append(ck)
        for i in range(n): Mk[i][i] += ck
    return c
def det(M):
    M = [row[:] for row in M]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f: M[r] = [M[r][j] - f * M[c][j] for j in range(n)]
    return d
def hurwitz_minors(c):
    d = len(c) - 1
    H = [[c[2 * i - j + 1] if 0 <= 2 * i - j + 1 <= d else Fr(0) for j in range(d)] for i in range(d)]
    return [det([row[:k] for row in H[:k]]) for k in range(1, d + 1)]
def classify(n, A):
    if not dom_free(n, A): return 'N', 'dominated'
    P = [[(Fr(1) if i == j else Fr(0)) - Fr(1, n) for j in range(n)] for i in range(n)]
    PAP = matmul(matmul(P, [[Fr(x) for x in r] for r in A]), P)
    c = charpoly(PAP)                       # degree n, highest first
    assert c[-1] == 0, 'PAP must have eigenvalue 0 (direction 1)'
    p = c[:-1]                              # chi_perp, degree n-1, monic
    d = n - 1
    # h(mu) = p(mu - 1/2): coefficient of mu^j
    h = [Fr(0)] * (d + 1)
    for k, pk in enumerate(p):              # pk * x^(d-k), x = mu - 1/2
        e = d - k
        for j in range(e + 1):
            h[j] += pk * comb(e, j) * Fr(-1, 2) ** (e - j)
    hc = list(reversed(h))                  # highest first
    assert hc[0] == 1
    if any(x < 0 for x in hc): return 'N', 'negative coefficient'
    if any(x == 0 for x in hc): return 'B', 'zero coefficient'
    mins = hurwitz_minors(hc)
    if any(m < 0 for m in mins): return 'N', 'negative minor'
    if any(m == 0 for m in mins): return 'B', 'zero minor'
    return 'S', 'hurwitz'
if __name__ == '__main__':
    for line in sys.stdin:
        s = line.split()
        if not s: continue
        n, A = decode_d6(s[0]); arcs = sum(map(sum, A))
        f, why = classify(n, A)
        print(s[0], arcs, f, why)
