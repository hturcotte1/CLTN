"""
ref_lib.py -- referee's INDEPENDENT exact tools (no prover code imported; no floating point in any decision).
Conventions CLAUDE.md s2: A[i][j]=1 iff i->j; D = J - I - A (non-edge matrix); beta = 1 + delta/eps, r = 1/beta.
"""
import sys
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')

z = sp.Symbol('z')
x = sp.Symbol('x')
kk = sp.Symbol('k')


# ---------------------------------------------------------------- construction (own implementation)
def arms_core(c):
    """Own builder.  c = [1,1,c2,...,cd].  Node order: c2 FIRST (index 0), then c1, c0, then arms in a different order
    than the prover's, so that agreement is not an artefact of the labelling.
    Returns (D as list of lists, parent dict (D-out-neighbour) for non-root nodes, root index, N set)."""
    d = len(c) - 1
    assert c[0] == 1 and c[1] == 1 and c[2] >= 1
    names = ['c2', 'c1', 'c0']
    succ = {'c1': 'c2', 'c0': 'c1'}          # D-out-neighbour of non-root nodes
    rootout = ['c0']                         # D-out-neighbours of the root c2
    # arms grouped by (exit, depth): exit c2 depth 1 for length 2; odd j: exit c1 depth (j-1)/2; even j>=4: exit c0 depth j/2-1
    groups = {}
    groups[('c2', 1)] = c[2] - 1
    for j in range(3, d + 1):
        if j % 2: groups[('c1', (j - 1) // 2)] = groups.get(('c1', (j - 1) // 2), 0) + c[j]
        else: groups[('c0', j // 2 - 1)] = groups.get(('c0', j // 2 - 1), 0) + c[j]
    spine = {}
    # create deepest-first is not possible (targets must exist) -> go by depth, but name arms by a counter in reverse
    cnt = 0
    for (ex, dep) in sorted(groups, key=lambda g: (g[1], g[0])):
        for i in range(groups[(ex, dep)]):
            a = 'A%d' % cnt; zz = 'Z%d' % cnt; cnt += 1
            if dep == 1: tgt = ex
            else: tgt = spine[(ex, dep - 1)]      # KeyError = forest constraint violated
            names += [zz, a]                       # z before a (different order from prover)
            succ[a] = zz; succ[zz] = tgt; rootout.append(a)
            spine.setdefault((ex, dep), a)
    idx = {v: i for i, v in enumerate(names)}
    m = len(names)
    D = [[0] * m for _ in range(m)]
    for u, v in succ.items(): D[idx[u]][idx[v]] = 1
    for v in rootout: D[0][idx[v]] = 1
    return D, idx


def random_intree(m, rng):
    """random in-tree rooted at 0: every node u != 0 gets one D-out-neighbour among earlier nodes (so paths reach 0);
    root gets a random nonempty out-set (excluding itself). Then randomly permute labels."""
    par = {u: rng.randrange(0, u) for u in range(1, m)}
    N = [v for v in range(1, m) if rng.random() < 0.4] or [rng.randrange(1, m)]
    perm = list(range(m)); rng.shuffle(perm)
    D = [[0] * m for _ in range(m)]
    for u, p in par.items(): D[perm[u]][perm[p]] = 1
    for v in N: D[perm[0]][perm[v]] = 1
    return D, perm[0]


def A_of_D(D):
    m = len(D); return [[0 if i == j else 1 - D[i][j] for j in range(m)] for i in range(m)]


def D_of_A(A):
    m = len(A); return [[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)]


def charpoly_sym(M):
    """sympy berkowitz charpoly as a sympy Poly in z (independent of the repo's Faddeev-LeVerrier)."""
    return sp.Matrix(M).charpoly(z)


def psi_sym(D):
    """psi_D(z) = 1^T adj(zI - D) 1 = chi_D - chi_{D+J}  (matrix determinant lemma, T2 with tau = -1)."""
    m = len(D)
    DJ = [[D[i][j] + 1 for j in range(m)] for i in range(m)]
    return sp.Poly(charpoly_sym(D).as_expr() - charpoly_sym(DJ).as_expr(), z)


def ell_and_N(D, root):
    m = len(D)
    out = {u: [v for v in range(m) if D[u][v]] for u in range(m)}
    for u in range(m):
        if u != root: assert len(out[u]) == 1
    ell = {root: 0}
    def L(u):
        if u in ell: return ell[u]
        seen = set(); path = []
        while u not in ell:
            assert u not in seen, 'cycle avoiding root'
            seen.add(u); path.append(u); u = out[u][0]
        b = ell[u]
        for v in reversed(path):
            b += 1; ell[v] = b
        return ell[path[0]]
    for u in range(m): L(u)
    return ell, set(out[root])


def Lam_Gam(D, root):
    ell, N = ell_and_N(D, root)
    m = len(D)
    Lam = 1 + sum(x ** (1 + ell[v]) for v in range(m) if v not in N)
    Gam = sum(x ** (1 + ell[w]) for w in N)
    return sp.expand(Lam), sp.expand(Gam), N


def to_z_poly(expr_x, deg):
    """z^deg * expr(1/z) as Poly in z"""
    return sp.Poly(sp.expand(z ** deg * expr_x.subs(x, 1 / z)), z)


# ---------------------------------------------------------------- reduced polynomials from c (Lemmas T, Q)
def reduced_from_c(c):
    """p (hi->lo), G (hi->lo), m from Lam/Gam computed from the arm counts (formula derived independently by referee:
    Lam = 1 + x + x^2 + sum_a x^len(a);  Gam = x^3 + sum_a x^(len(a)+1);  |N| = 1 + #arms)."""
    d = len(c) - 1
    narms = (c[2] - 1) + sum(c[3:])
    m = 3 + 2 * narms
    Lam = 1 + x + x ** 2 + (c[2] - 1) * x ** 2 + sum(c[j] * x ** j for j in range(3, d + 1))
    Gam = x ** 3 + (c[2] - 1) * x ** 3 + sum(c[j] * x ** (j + 1) for j in range(3, d + 1))
    Nn = 1 + narms
    num = sp.expand((m + 1 - Nn) * (1 - Gam) + (Nn - 1) * Lam)
    Q, rem = sp.div(sp.Poly(num, x), sp.Poly(1 - x, x))
    assert rem.is_zero
    Qx = Q.as_expr()
    G = to_z_poly(Qx, d)       # degree d, leading coeff m
    p = sp.Poly(sum(c[j] * z ** (d - j) for j in range(d + 1)), z)
    return p, G, m


# ---------------------------------------------------------------- exact strip test (Hurwitz)
def hurwitz_minors(coeffs):
    """coeffs hi->lo (sympy exprs), returns list of leading principal Hurwitz minors Delta_1..Delta_n."""
    n = len(coeffs) - 1
    a = lambda i: coeffs[i] if 0 <= i <= n else 0
    H = sp.Matrix(n, n, lambda i, j: a(2 * (i + 1) - (j + 1)))
    return [sp.expand(H[:i, :i].det(method='berkowitz')) for i in range(1, n + 1)]


def shifted(poly_z, s):
    """Q(y) = (-1)^d poly(-y - s): roots y = -z - s. All roots of poly have Re z > -s iff Q Hurwitz-stable."""
    y = sp.Symbol('y')
    d = poly_z.degree()
    Q = sp.Poly(sp.expand((-1) ** d * poly_z.as_expr().subs(z, -y - s)), y)
    return Q


def strip_right_of(poly_z, s):
    """three-valued exact: all roots Re z > -s ? (True/False/None)"""
    Q = shifted(poly_z, s)
    co = Q.all_coeffs()
    if co[0] < 0: co = [-a for a in co]
    if any(a == 0 for a in co): return None
    if any(a < 0 for a in co): return False
    for mnr in hurwitz_minors(co):
        if mnr == 0: return None
        if mnr < 0: return False
    return True


def k_threshold(p, G, s):
    """For H_k = G + k p, find exact least integer k0 >= 1 such that for EVERY real k >= k0 all coefficients and all
    Hurwitz minors of the shifted polynomial are > 0 (=> all roots of H_k have Re > -s for every integer k >= k0).
    Uses exact real-root isolation of each minor as a polynomial in k."""
    H = sp.Poly(G.as_expr() + kk * p.as_expr(), z)
    Q = shifted(H, s)
    y = Q.gens[0]
    co = [sp.expand(a) for a in sp.Poly(Q.as_expr(), y).all_coeffs()]
    polys = co + hurwitz_minors(co)
    bound = F(1)
    detail = []
    for e in polys:
        P = sp.Poly(e, kk)
        lc = P.LC()
        if lc <= 0:
            return None, 'leading coefficient in k not positive: %s' % lc
        rr = P.real_roots() if P.degree() > 0 else []
        top = None
        if rr:
            # exact: isolate largest real root with rational upper bound
            ivs = P.intervals()
            top = max(F(str(b)) for (a, b), mult in ivs)
            # refine to make sure the bound is tight enough
            ivs = P.intervals(eps=sp.Rational(1, 10 ** 6))
            top = max(F(str(b)) for (a, b), mult in ivs)
        detail.append((P.degree(), str(top)))
        if top is not None and top >= bound: bound = top
    # every poly positive for real k > bound (leading coeff > 0 and no real root above bound)
    k0 = int(bound) + 1 if bound == int(bound) else -(-bound.numerator // bound.denominator)
    # bound may be a root itself; verify integer k0 strictly passes, and find the exact least integer by testing downward
    def passes(k):
        vals = [sp.Poly(e, kk).eval(k) for e in polys]
        return all(v > 0 for v in vals)
    assert passes(k0)
    # integers above `bound` all pass (no real root above it, positive leading coefficients); below it test integers
    while k0 > 1 and passes(k0 - 1):
        k0 -= 1
    return k0, passes, detail
