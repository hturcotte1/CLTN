"""
fixlib.py -- exact algebra for clique-class families at a FIXED legal parameter point (eps, delta).

Conventions: CLAUDE.md section 2.  A[i][j] = 1 iff i -> j.  D := J - I - A (non-edge digraph, D[i][j] = 1 iff i != j
and i -/-> j).  t = (1-eps)/eps, beta = 1 + delta/eps, r = 1/beta, a := t/beta.  Legal iff t(beta - 2) > 1.

Key identity (proof in notes_fixed.md, Lemma S):  the full support of a digraph G' is positive-stable at (eps, delta)
iff every eigenvalue of  a J + D'  has real part > -r.

Sink-class families (Lemma X1' of notes_theory.md, re-derived in notes_fixed.md, Lemma P): core G (m nodes) plus
classes V_j of sizes k_j (bidirected cliques, mutually bidirected, every class node sends to every other node,
node i of G sends to V_j iff S_j[i] = 1; b_j := 1 - S_j).  With n = m + p, p = sum k_j:

    det(z I_n - a J_n - D') = z^(p-1) * Psi(z),
    Psi(z) := z * det(zI - D - aJ) - a * sum_j k_j * det(zI - D + b_j 1^T)            (degree m+1, monic)

No floating point in any decision.  Functions whose names start with heur_ use floats to steer searches only.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
import ctln

z, k, x, a_s, r_s = sp.symbols('z k x a r')


def from_d6(s):
    """digraph6 decoder (nauty format); handles n <= 62 (one header byte) and 63 <= n <= 258047 ('~' + 3 bytes)."""
    if s[1] == '~':
        n = ((ord(s[2]) - 63) << 12) | ((ord(s[3]) - 63) << 6) | (ord(s[4]) - 63); p = s[5:]
    else:
        n = ord(s[1]) - 63; p = s[2:]
    A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def to_d6(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    if n <= 62:
        s = '&' + chr(n + 63)
    else:
        s = '&~' + chr(63 + ((n >> 12) & 63)) + chr(63 + ((n >> 6) & 63)) + chr(63 + (n & 63))
    for kk in range(0, len(bits), 6):
        v = 0
        for bb in bits[kk:kk + 6]: v = 2 * v + bb
        s += chr(v + 63)
    return s


def nonedge(A):
    m = len(A)
    return [[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)]


def build(A, classes, sizes):
    """explicit graph: core A plus sink classes; classes = list of in-sets S_j (0/1 tuples), sizes = list of k_j."""
    m = len(A); n = m + sum(sizes)
    G = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): G[i][j] = A[i][j]
    idx = m
    for S, kj in zip(classes, sizes):
        for _ in range(kj):
            v = idx; idx += 1
            for i in range(m):
                G[i][v] = S[i]      # i -> v iff S[i] = 1
                G[v][i] = 1         # class node sends to every core node
    for u in range(m, n):
        for v in range(m, n):
            if u != v: G[u][v] = 1  # classes are bidirected cliques, mutually bidirected
    return G


def charpoly_sym(M):
    """det(z I - M) for a sympy Matrix M (entries may contain symbols)."""
    m = M.shape[0]
    return sp.expand((z * sp.eye(m) - M).det(method='berkowitz'))


def psi_parts(A, classes):
    """returns (P0, [P_j]) with Psi = P0 - a * sum_j k_j * P_j,
    P0 = z * det(zI - D - aJ) (polynomial in z and a), P_j = det(zI - D + b_j 1^T) = chi_{L_j} (integer poly in z)."""
    m = len(A)
    D = sp.Matrix(nonedge(A))
    P0 = sp.expand(z * charpoly_sym(D + a_s * sp.ones(m, m)))
    Ps = []
    for S in classes:
        b = sp.Matrix([1 - S[i] for i in range(m)])
        Ps.append(charpoly_sym(D - b * sp.ones(1, m)))
    return P0, Ps


def Psi_expr(A, classes, rhos):
    """Psi as an expression in z, a, k with class sizes k_j = rho_j * k."""
    P0, Ps = psi_parts(A, classes)
    return sp.expand(P0 - a_s * k * sum(rho * P for rho, P in zip(rhos, Ps)))


def check_identity(A, classes, sizes, avals=(sp.Rational(3, 7), sp.Integer(5))):
    """exact check of det(zI - aJ_n - D') = z^(p-1) Psi(z) on the explicit graph, for given sizes and a few rational a."""
    G = build(A, classes, sizes); n = len(G); p = sum(sizes)
    Dp = sp.Matrix(nonedge(G))
    P0, Ps = psi_parts(A, classes)
    for av in avals:
        lhs = sp.Poly(charpoly_sym(Dp + av * sp.ones(n, n)), z)
        Psi = sp.expand((P0 - a_s * sum(kj * P for kj, P in zip(sizes, Ps))).subs(a_s, av))
        rhs = sp.Poly(sp.expand(z ** (p - 1) * Psi), z) if p >= 1 else None
        if lhs != rhs:
            return False
    return True


def hurwitz_minors(coeffs):
    d = len(coeffs) - 1
    aa = lambda i: coeffs[i] if 0 <= i <= d else 0
    H = sp.Matrix(d, d, lambda i, j: aa(2 * (i + 1) - (j + 1)))
    return [sp.expand(H[:j, :j].det(method='berkowitz')) for j in range(1, d + 1)]


def shifted_coeffs(Psi, r, extra_subs=None):
    """h(x) := (-1)^deg Psi(-r - x): roots x = -r - z; Re z > -r  <=>  Re x < 0.  Returns hi->lo coefficient list
    (expressions in the remaining symbols).  Psi is monic in z, so h is monic in x."""
    e = Psi.subs(extra_subs or {})
    deg = sp.Poly(e, z).degree()
    h = sp.Poly(sp.expand((-1) ** deg * e.subs(z, -r - x)), x)
    co = [sp.expand(c) for c in h.all_coeffs()]
    assert co[0] == 1
    return co


def positive_for_all_k_ge(p, k0):
    """EXACT: polynomial p(k) (sympy expr) is > 0 for every real k >= k0.  Returns (bool, info)."""
    P = sp.Poly(p, k)
    if P.is_zero: return False, 'identically zero'
    if P.degree() == 0: return (P.LC() > 0), 'constant %s' % P.LC()
    if P.LC() <= 0: return False, 'leading coefficient not positive'
    if P.eval(k0) <= 0: return False, 'value at k0 not positive'
    nr = P.count_roots(k0, None)          # exact (Sturm) count of real roots in [k0, inf)
    return nr == 0, 'roots in [k0,inf): %d' % nr


def smallest_k0(p, kmin=0):
    """EXACT: smallest integer k0 >= kmin such that p(k) > 0 for all real k >= k0 (None if LC <= 0)."""
    P = sp.Poly(p, k)
    if P.is_zero: return None
    if P.degree() == 0: return kmin if P.LC() > 0 else None
    if P.LC() <= 0: return None
    rr = P.real_roots()
    k0 = kmin
    if rr:
        top = max(rr)
        k0 = max(kmin, int(sp.floor(top)) + 1)
    ok, _ = positive_for_all_k_ge(p, k0)
    assert ok
    return k0


def fixed_point_minors(A, classes, rhos, beta, t):
    """Hurwitz data (as polynomials in k) of the family at the fixed point (beta, t); exact rationals."""
    beta = sp.Rational(beta); t = sp.Rational(t)
    r = 1 / beta; av = t / beta
    Psi = Psi_expr(A, classes, rhos).subs(a_s, av)
    co = shifted_coeffs(Psi, r)
    mins = hurwitz_minors(co)
    return co, mins


def k_threshold(A, classes, rhos, beta, t, kmin=0):
    """smallest integer k0 >= kmin such that every Hurwitz minor is > 0 for ALL real k >= k0, i.e. the family is
    stable at (beta, t) for every integer k >= k0.  None if some minor has non-positive leading coefficient in k
    (then the family is NOT stable for all large k: see notes, Lemma E)."""
    co, mins = fixed_point_minors(A, classes, rhos, beta, t)
    k0 = kmin
    for mn in mins:
        kk = smallest_k0(mn, kmin)
        if kk is None: return None, mins
        k0 = max(k0, kk)
    return k0, mins


def stable_at_k(A, classes, rhos, beta, t, kval):
    """exact three-valued Hurwitz decision of the family member at integer k (via Psi)."""
    beta = F(beta); t = F(t)
    co, mins = fixed_point_minors(A, classes, rhos, beta, t)
    coeffs = [F(str(sp.Rational(c.subs(k, kval)))) for c in co]
    return ctln.hurwitz_stable(coeffs)


def positivity_core(A, classes, beta):
    """EXACT (Lemma X2'): every member of the family is positive at beta iff w_G(beta) > 0 and sigma_{S_j}(beta) > 1 for all j.
    Returns (bool, w, sigmas)."""
    beta = F(beta); g = 1 - 1 / beta; m = len(A)
    M = [[F(A[j][i]) + (g if i == j else 0) for j in range(m)] for i in range(m)]
    w = ctln.solve_exact(M, [F(1)] * m)
    if w is None: return False, None, None
    sig = [sum((w[i] for i in range(m) if S[i]), F(0)) for S in classes]
    ok = all(v > 0 for v in w) and all(s > 1 for s in sig)
    return ok, w, sig


def eps_delta(beta, t):
    beta = F(beta); t = F(t)
    return 1 / (1 + t), (beta - 1) / (1 + t)


def certify_member(A, classes, sizes, beta, t):
    """ctln.analyze_support on the explicit graph at the exact point.  Returns the record (dict)."""
    G = build(A, classes, sizes); n = len(G)
    eps, delta = eps_delta(beta, t)
    rec = ctln.analyze_support(G, list(range(n)), eps, delta)
    rec['legal'] = ctln.legal(eps, delta)
    rec['n'] = n; rec['d6'] = to_d6(G); rec['sizes'] = list(sizes)
    return rec


# ------------------------------------------------------------- float heuristics (search steering only)
def heur_minre(poly_expr):
    import numpy as np
    P = sp.Poly(poly_expr, z)
    rts = np.roots([float(c) for c in P.all_coeffs()])
    return min(rr.real for rr in rts)


# ------------------------------------------------------------- region statements (uniform in t and k)
u_s, v_s = sp.symbols('u v')


def reduced_psi(A, classes, rhos):
    """Psi with the common factor z^j0 removed (roots z = 0 are eigenvalues 0 > -r, harmless).  Returns (g, j0)."""
    Psi = Psi_expr(A, classes, rhos)
    P = sp.Poly(Psi, z)
    co = P.all_coeffs()[::-1]           # lo -> hi
    j0 = 0
    while sp.expand(co[j0]) == 0: j0 += 1
    g = sp.expand(sum(c * z ** (i - j0) for i, c in enumerate(co) if i >= j0))
    return g, j0


def legal_a_bound(r):
    """for fixed r = 1/beta in (0, 1/2): (beta, t) legal  <=>  a = t/beta > r^2/(1 - 2r)."""
    r = sp.Rational(r)
    return r ** 2 / (1 - 2 * r)


def region_minors(g, r1):
    """Hurwitz minors (polynomials in a, k) of h(x) = (-1)^deg g(-r1 - x)."""
    deg = sp.Poly(g, z).degree()
    h = sp.Poly(sp.expand((-1) ** deg * g.subs(z, -sp.Rational(r1) - x)), x)
    co = [sp.expand(c) for c in h.all_coeffs()]
    assert co[0] == 1
    return co, hurwitz_minors(co)


def sign_pattern_certificate(p, a1, k0):
    """EXACT sufficient test: p(a1 + u, k0 + v) has all coefficients >= 0 (as a polynomial in u, v) and a positive constant
    term  =>  p > 0 on {a >= a1, k >= k0}.  Returns (ok, n_negative_coeffs, constant)."""
    q = sp.Poly(sp.expand(p.subs({a_s: sp.Rational(a1) + u_s, k: sp.Rational(k0) + v_s})), u_s, v_s)
    cs = q.coeffs()
    const = q.coeff_monomial(1)
    neg = sum(1 for c in cs if c < 0)
    return (neg == 0 and const > 0), neg, const


def check_identity_fast(A, classes, sizes, avals=(F(3, 7), F(5), F(200, 21))):
    """EXACT check of Lemma P on the explicit graph by evaluation: both sides are monic polynomials of degree n in z, so
    agreement at n + 1 distinct rational points z proves the identity (for each rational a tested).  Fractions only."""
    G = build(A, classes, sizes); n = len(G); p = sum(sizes)
    Dp = nonedge(G)
    P0, Ps = psi_parts(A, classes)
    for av in avals:
        av_s = sp.Rational(av.numerator, av.denominator)
        PsiP = sp.Poly(sp.expand((P0 - a_s * sum(kj * P for kj, P in zip(sizes, Ps))).subs(a_s, av_s)), z)
        for zi in range(-(n // 2) - 1, n - n // 2 + 1):          # n + 2 distinct points zi + 1/3 (never 0)
            zz = F(zi) + F(1, 3)
            M = [[(zz if i == j else F(0)) - av - Dp[i][j] for j in range(n)] for i in range(n)]
            lhs = ctln.det_exact(M)
            rhs = (zz ** (p - 1)) * F(str(PsiP.eval(sp.Rational(zz.numerator, zz.denominator))))
            if lhs != rhs:
                return False
    return True
