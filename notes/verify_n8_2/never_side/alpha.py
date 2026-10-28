"""Own exact classification of alpha := spectral abscissa of A restricted to 1^perp, relative to -1/2.

Construction (independent of filter8's identity): B = [e_k - e_(k+1)]_(k=0..n-2) is a rational basis of 1^perp;
the compression of A to 1^perp in that basis is Cm = (B^T B)^{-1} B^T A B, similar to U^T A U for any
orthonormal U, so chi_perp = charpoly(Cm) exactly (sympy over QQ).  Cross-check: n*chi_perp == chi_A - chi_{A+J}
(the identity the C comment states), both sides by own integer Faddeev-LeVerrier.
g(mu) := n 2^(n-1) chi_perp(mu - 1/2) (integer; equals the C code g).  alpha < -1/2 iff g is Hurwitz (all coefficients and all leading
principal minors > 0).  Otherwise let u = gcd(g(mu), g(-mu)) (contains every imaginary-axis root of g with full
multiplicity, plus the +-r pairs), g1 = g/u (no axis root).  alpha <= -1/2 iff g1 is Hurwitz and u = mu^a U(mu^2)
with every root of U real and <= 0 (then every root of u lies on the axis).  Combined with "not Hurwitz" this gives
alpha = -1/2 exactly; otherwise alpha > -1/2.
"""
import sympy as sp
from sympy import Poly, QQ, Rational, symbols, gcd, Matrix
from engine import charpoly_int, hurwitz
mu, z = symbols('mu z')

def chi_perp_two_ways(A):
    n = len(A)
    B = Matrix(n, n - 1, lambda i, k: 1 if i == k else (-1 if i == k + 1 else 0))
    Am = Matrix(A)
    Cm = (B.T * B).inv() * (B.T * Am * B)
    cp = Cm.charpoly(mu)                      # monic, degree n-1, rational coefficients
    chi_perp = Poly(cp.as_expr(), mu, domain=QQ)
    cA = charpoly_int(A)
    cAJ = charpoly_int([[A[i][j] + 1 for j in range(n)] for i in range(n)])
    diff = [cA[k] - cAJ[k] for k in range(n + 1)]   # coefficients of chi_A - chi_{A+J}, highest first
    assert diff[0] == 0
    ident = Poly(diff[1:], mu, domain=QQ)             # degree n-1, leading coefficient n
    assert ident == chi_perp * n, "identity n*chi_perp == chi_A - chi_{A+J} FAILED"
    return chi_perp

def classify_alpha(A):
    """returns ('<' | '=' | '>', detail)  for alpha vs -1/2, plus the Hurwitz detail on g."""
    n = len(A)
    chi = chi_perp_two_ways(A)
    g = Poly((chi.as_expr().subs(mu, mu - Rational(1, 2)) * (n * 2 ** (n - 1))).expand(), mu, domain=QQ)
    gc = [sp.Rational(c) for c in g.all_coeffs()]
    assert all(c.q == 1 for c in gc), "g not integral"
    gi = [int(c.p) for c in gc]
    hv, hdet = hurwitz(gi)
    if hv is True:
        return '<', hdet, gi
    gneg = Poly(g.as_expr().subs(mu, -mu).expand(), mu, domain=QQ)
    u = gcd(g, gneg)
    g1, rem = sp.div(g, u)
    assert rem.is_zero
    g1i = [int(c * sp.lcm([sp.Rational(x).q for x in g1.all_coeffs()])) for c in g1.all_coeffs()]
    if g1i[0] < 0: g1i = [-c for c in g1i]
    h1 = hurwitz(g1i)[0] if g1.degree() >= 1 else True
    if h1 is not True:
        return '>', f"g1 (axis-free part, degree {g1.degree()}) not Hurwitz: {hurwitz(g1i)[1]}", gi
    # u = mu^a U(mu^2): check structure and that every root of U is real and <= 0
    ucoeffs = u.all_coeffs()[::-1]         # low -> high
    a = 0
    while ucoeffs[a] == 0: a += 1
    rest = ucoeffs[a:]
    if any(rest[k] != 0 for k in range(1, len(rest), 2)):
        return '>', "u is not mu^a U(mu^2) -- unexpected", gi
    U = Poly(rest[::2][::-1], z, domain=QQ)
    if U.degree() >= 1:
        for f, mult in U.sqf_list()[1]:
            d = f.degree()
            if f.count_roots() != d or f.count_roots(sup=0) != d:
                return '>', "u has a root off the imaginary axis (U has a non-real or positive root)", gi
    return '=', f"g not Hurwitz ({hdet}); axis-free part Hurwitz; all {u.degree()} roots of u on the axis", gi

if __name__ == "__main__":
    import sys, time
    from engine import parse_d6
    t0 = time.time(); cnt = {}
    for line in open(sys.argv[1]):
        p = line.split()
        if not p: continue
        c, det, _ = classify_alpha(parse_d6(p[0]))
        flag = p[2] if len(p) > 2 else 'none'
        cnt[(flag, c)] = cnt.get((flag, c), 0) + 1
    print(cnt, f"{time.time()-t0:.1f}s")
