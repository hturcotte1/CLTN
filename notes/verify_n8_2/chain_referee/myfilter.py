"""Referee's own exact implementation of filter8's two necessary conditions, written from the definitions:
 F1 domination-free (Curto: k dominates j iff j->k, not k->j, and every in-neighbour i of j (i not in {j,k}) is an in-neighbour of k)
 F2 alpha(A|perp) < -1/2 strictly, decided by exact Hurwitz minors (sympy Rational) of chi_perp(mu - 1/2), where chi_perp is
    the characteristic polynomial of A compressed to 1^perp, computed INDEPENDENTLY as det(lam I - P A P)/lam with P = I - J/n.
Flags: S (strict), B (a minor or coefficient vanishes; g not Hurwitz => alpha >= -1/2), N (non-survivor: F1 fails or alpha > -1/2 decidedly)."""
import sys, sympy as sp
from sympy import Matrix, Rational, symbols, Poly, eye, ones
lam, mu = symbols('lam mu')
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
def hurwitz_minors(coeffs):
    """coeffs highest first, a0 > 0.  Returns list of leading principal minors of the Hurwitz matrix (exact)."""
    d = len(coeffs) - 1
    H = [[coeffs[2 * i - j + 1] if 0 <= 2 * i - j + 1 <= d else 0 for j in range(d)] for i in range(d)]
    return [Matrix(H)[:k, :k].det() for k in range(1, d + 1)]
def classify(n, A):
    if not dom_free(n, A): return 'N', 'dominated'
    M = Matrix(A); P = eye(n) - ones(n, n) / n
    cp = (lam * eye(n) - P * M * P).det()
    chi_perp = sp.cancel(cp / lam)
    assert Poly(chi_perp, lam).degree() == n - 1
    g = Poly(sp.expand(chi_perp.subs(lam, mu - Rational(1, 2))), mu)
    c = g.all_coeffs()
    assert c[0] == 1
    if any(x < 0 for x in c): return 'N', 'negative coefficient'
    if any(x == 0 for x in c): return 'B', 'zero coefficient'
    mins = hurwitz_minors(c)
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
