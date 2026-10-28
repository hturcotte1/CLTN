"""
engine.py -- agent interval_all_t: independent exact decision engine.

Everything here is written from the definitions in the task / CLAUDE.md s2:
  A[i][j] = 1 iff arc i->j
  W[i][j] = 0 (i==j); -1+eps if A[j][i]==1; -1-delta otherwise
  full support {0..n-1} stable fixed point iff
    (1) (I-W) x = 1 has a solution with x_i > 0 for all i
    (3) all eigenvalues of M = -I + W have Re < 0  (Hurwitz minors of charpoly of M)
  three-valued: a vanishing Hurwitz minor or x_i == 0 or singular (I-W) => UNDECIDED.

Limit ray: eps = 1/(t+1), delta = (beta-1) eps, beta = b/d rational.
Then (t+1)(I-W) = K(t) := I + beta N + t J with N = J - I - A^T  (checked below
by build_K_vs_W_check against the definition, symbolically, for every graph).
M = -I + W = -(I-W) = -K/(t+1), so d(t+1) M = -K' with K' = d K integer-polynomial.
Char poly of M, p_M(l) = sum a_k(M) l^{n-k} relates to that of -K' by
  a_k(M) = a_k(-K') / (d(t+1))^k   and   H_k(M) = H_k(-K') / (d(t+1))^{k(k+1)/2}
(each term of the k-th Hurwitz determinant has index sum k(k+1)/2), so for t > 0
sign H_k(M) = sign H_k(-K').  x = (t+1) K^{-1} 1 = (t+1) d adj(K') 1 / det K'.
All decisions are therefore made on integer/rational polynomials in t.
No floats anywhere.
"""
from fractions import Fraction as F
import sympy
from sympy import Symbol, Poly, QQ, ZZ

T = Symbol('t')

# ----------------------------------------------------------------- digraph6
def parse_digraph6(s):
    assert s[0] == '&', s
    s = s[1:]
    n = ord(s[0]) - 63
    assert 0 <= n <= 62, "only small n handled"
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        assert 0 <= v < 64
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    assert len(bits) >= n * n
    A = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        assert A[i][i] == 0
    return n, A

# ----------------------------------------------------------------- polynomials in t (coeff lists, low->high, Fractions)
def pnorm(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p

def padd(p, q):
    r = [F(0)] * max(len(p), len(q))
    for i, c in enumerate(p): r[i] += c
    for i, c in enumerate(q): r[i] += c
    return pnorm(r)

def psub(p, q):
    r = [F(0)] * max(len(p), len(q))
    for i, c in enumerate(p): r[i] += c
    for i, c in enumerate(q): r[i] -= c
    return pnorm(r)

def pmul(p, q):
    if not p or not q: return []
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0: continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return pnorm(r)

def pscale(p, c):
    return pnorm([a * c for a in p])

def pdivexact(p, q):
    """exact division p / q in Q[t]; asserts zero remainder"""
    p = pnorm(p); q = pnorm(q)
    assert q, "division by zero polynomial"
    if not p: return []
    assert len(p) >= len(q)
    rem = list(p)
    quot = [F(0)] * (len(p) - len(q) + 1)
    lq = q[-1]
    for k in range(len(p) - len(q), -1, -1):
        c = rem[k + len(q) - 1] / lq
        quot[k] = c
        if c != 0:
            for j, b in enumerate(q):
                rem[k + j] -= c * b
    assert all(x == 0 for x in rem), "non-exact division"
    return pnorm(quot)

def peval(p, x):
    r = F(0)
    for c in reversed(p):
        r = r * x + c
    return r

def pconst(c):
    return pnorm([F(c)])

def to_sympy(p):
    return Poly(sum(sympy.Rational(c.numerator, c.denominator) * T**i for i, c in enumerate(p)), T, domain=QQ)

# ----------------------------------------------------------------- polynomial matrices
def mat_det_bareiss(Mx):
    """Fraction-free Bareiss determinant over Q[t] (entries = poly lists)."""
    n = len(Mx)
    M = [[list(e) for e in row] for row in Mx]
    sign = 1
    prev = [F(1)]
    for k in range(n - 1):
        # pivot
        if not M[k][k]:
            sw = None
            for i in range(k + 1, n):
                if M[i][k]:
                    sw = i; break
            if sw is None:
                return []
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = psub(pmul(M[i][j], M[k][k]), pmul(M[i][k], M[k][j]))
                M[i][j] = pdivexact(num, prev) if num else []
            M[i][k] = []
        prev = M[k][k]
    d = M[n - 1][n - 1]
    return pscale(d, sign) if sign == -1 else d

def charpoly_FL(Mx):
    """Faddeev-LeVerrier: returns [c_0=1, c_1, ..., c_n] with det(l I - Mx) = sum c_k l^{n-k}; entries poly lists."""
    n = len(Mx)
    I = [[pconst(1) if i == j else [] for j in range(n)] for i in range(n)]
    Mk = [[[] for _ in range(n)] for _ in range(n)]
    c = [pconst(1)]
    for k in range(1, n + 1):
        # Mk = A*M_{k-1} + c_{k-1} I
        AM = [[[] for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                s = []
                for l in range(n):
                    if Mx[i][l] and Mk[l][j]:
                        s = padd(s, pmul(Mx[i][l], Mk[l][j]))
                if i == j:
                    s = padd(s, c[k - 1])
                AM[i][j] = s
        Mk = AM
        tr = []
        for i in range(n):
            for l in range(n):
                if Mx[i][l] and Mk[l][i]:
                    tr = padd(tr, pmul(Mx[i][l], Mk[l][i]))
        c.append(pscale(tr, F(-1, k)))
    return c

def hurwitz_matrix(a, n):
    """a = [a_0..a_n] coefficient polys of p(l) = sum a_k l^{n-k}. Hurwitz matrix H[i][j] = a_{2(i+1)-(j+1)} (0 if out of range)."""
    H = [[[] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            idx = 2 * (i + 1) - (j + 1)
            if 0 <= idx <= n:
                H[i][j] = list(a[idx])
    return H

def hurwitz_minors(a, n):
    H = hurwitz_matrix(a, n)
    return [mat_det_bareiss([row[:k] for row in H[:k]]) for k in range(1, n + 1)]

# ----------------------------------------------------------------- symbolic-in-t build
def build_Kprime(A, b, d):
    """K' = d*K, K = I + beta N + t J, N = J - I - A^T, beta = b/d. Entries poly lists (ints as Fractions)."""
    n = len(A)
    Kp = []
    for i in range(n):
        row = []
        for j in range(n):
            const = d if i == j else 0
            if i != j and A[j][i] == 0:
                const += b
            row.append(pnorm([F(const), F(d)]))
        Kp.append(row)
    return Kp

def build_K_vs_W_check(A, b, d):
    """Symbolic check against the definition: (t+1)*d*(I - W) == K' with eps=1/(t+1), delta=(beta-1)eps."""
    n = len(A)
    beta = sympy.Rational(b, d)
    eps = 1 / (T + 1); delta = (beta - 1) * eps
    ok = True
    Kp = build_Kprime(A, b, d)
    for i in range(n):
        for j in range(n):
            if i == j: w = 0
            elif A[j][i] == 1: w = -1 + eps
            else: w = -1 - delta
            lhs = sympy.together((T + 1) * d * ((1 if i == j else 0) - w))
            rhs = to_sympy(Kp[i][j]).as_expr()
            if sympy.simplify(lhs - rhs) != 0:
                ok = False
    return ok

def symbolic_conditions(A, beta):
    """Returns dict of named condition polynomials (poly lists) that must all be > 0 for stability+positivity+legality,
    plus bookkeeping."""
    n = len(A)
    b, d = beta.numerator, beta.denominator
    Kp = build_Kprime(A, b, d)
    negKp = [[pscale(e, -1) for e in row] for row in Kp]
    a = charpoly_FL(negKp)             # char poly of -K'
    H = hurwitz_minors(a, n)           # H_1..H_n of -K'
    D = mat_det_bareiss(Kp)            # det K'
    # check: a_n = det(-( -K')) = det K'  (p(0) = det(0 I - (-K')) = det(K'))
    assert psub(a[n], D) == [], "a_n != det K'"
    assert psub(H[n - 1], pmul(D, H[n - 2])) == [], "H_n != a_n H_{n-1}"
    # Cramer numerators N_i = det(K' with column i replaced by ones)
    N = []
    for i in range(n):
        Mi = [[list(Kp[r][c]) if c != i else pconst(1) for c in range(n)] for r in range(n)]
        N.append(mat_det_bareiss(Mi))
    # legality: q = (beta-1) t/(t+1) > 1  <=>  (b-2d) t - d > 0
    L = pnorm([F(-d), F(b - 2 * d)])
    conds = {}
    for k in range(1, n):
        conds['H%d' % k] = H[k - 1]
    conds['D'] = D                       # det K' > 0 (together with H_{n-1} > 0 gives H_n > 0)
    for i in range(n):
        conds['N%d' % i] = N[i]          # x_i = (t+1) d N_i / D > 0  <=> N_i > 0 given D > 0
    conds['L'] = L
    return conds, {'a': a, 'H': H, 'D': D, 'N': N}

# ----------------------------------------------------------------- real root isolation / sign set
def positive_root_intervals(conds):
    """Factor every condition polynomial over Q; isolate the positive real roots of every distinct irreducible factor
    with sympy's exact isolation; refine until intervals of different factors are disjoint.
    Returns (roots, facts) where roots = sorted list of (lo, hi, factor_index) with lo<=root<=hi rationals and
    facts = list of sympy Poly (primitive, integer)."""
    facts = []
    fact_of = {}
    for name, p in conds.items():
        assert p, "condition polynomial %s is identically zero" % name
        sp = to_sympy(p)
        _, fl = sp.factor_list()
        names = []
        for f, m in fl:
            f = Poly(f, T, domain=ZZ)
            if f.degree() < 1: continue
            # normalise sign: positive leading coefficient
            if f.LC() < 0: f = -f
            key = tuple(f.all_coeffs())
            if key not in fact_of:
                fact_of[key] = len(facts); facts.append(f)
            names.append((fact_of[key], m))
        # keep for reporting
        conds[name] = p
    roots = []   # (lo, hi, fi)
    eps = sympy.Rational(1, 2**20)
    ivs = {}
    for fi, f in enumerate(facts):
        ivs[fi] = [(sympy.Rational(lo), sympy.Rational(hi)) for (lo, hi), m in f.intervals(inf=0, eps=eps)]
    # refine until pairwise disjoint across different factors
    for _ in range(200):
        allr = sorted((lo, hi, fi) for fi, L in ivs.items() for (lo, hi) in L)
        bad = set()
        for u in range(len(allr) - 1):
            lo1, hi1, f1 = allr[u]; lo2, hi2, f2 = allr[u + 1]
            if hi1 >= lo2:   # overlap or touching (same factor: sympy guarantees disjointness, touching is refined away)
                bad.add(f1); bad.add(f2)
        if not bad:
            break
        eps = eps / 2**20
        for fi in bad:
            ivs[fi] = [(sympy.Rational(lo), sympy.Rational(hi)) for (lo, hi), m in facts[fi].intervals(inf=0, eps=eps)]
    else:
        raise RuntimeError("could not separate roots (common root between supposedly distinct irreducibles?)")
    allr = sorted((lo, hi, fi) for fi, L in ivs.items() for (lo, hi) in L)
    # drop roots at t <= 0 (intervals returned with inf=0 may include 0 itself)
    allr = [(lo, hi, fi) for (lo, hi, fi) in allr if hi > 0]
    return [(F(int(lo.p), int(lo.q)), F(int(hi.p), int(hi.q)), fi) for lo, hi, fi in allr], facts

def sign_set(conds):
    """Exact set S = {t > 0 : every condition polynomial > 0}, as a list of open intervals between consecutive roots.
    Returns dict."""
    roots, facts = positive_root_intervals(conds)
    # gap sample points
    pts = []
    bounds = []
    prev_hi = F(0)
    for (lo, hi, fi) in roots:
        if lo > prev_hi:
            pts.append((prev_hi + lo) / 2); bounds.append((prev_hi, lo))
        else:
            # root interval touching 0 (lo == 0 == prev_hi): no gap
            pts.append(None); bounds.append((prev_hi, lo))
        prev_hi = hi
    pts.append(prev_hi + 1); bounds.append((prev_hi, None))
    gaps = []
    for g, x in enumerate(pts):
        if x is None:
            gaps.append(None); continue
        bad = [name for name, p in conds.items() if peval(p, x) <= 0]
        gaps.append(bad)
    # gap g is between root g-1 (or 0) and root g (or infinity)
    S = []
    for g, bad in enumerate(gaps):
        if bad is None or bad: continue
        left = (roots[g - 1] if g > 0 else None)
        right = (roots[g] if g < len(roots) else None)
        S.append((left, right))
    # which condition polys vanish at each root
    root_conds = []
    for (lo, hi, fi) in roots:
        f = facts[fi]
        vanish = []
        for name, p in conds.items():
            q, r = to_sympy(p).div(Poly(f, T, domain=QQ))
            if r.is_zero: vanish.append(name)
        root_conds.append(vanish)
    return {'roots': roots, 'facts': facts, 'gaps': gaps, 'S': S, 'root_conds': root_conds}

# ----------------------------------------------------------------- pointwise exact engine (Fractions), from the definition
def build_W(A, eps, delta):
    n = len(A)
    W = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            W[i][j] = (-1 + eps) if A[j][i] == 1 else (-1 - delta)
    return W

def solve_frac(Mx, rhs):
    n = len(Mx)
    M = [list(map(F, row)) + [F(r)] for row, r in zip(Mx, rhs)]
    for k in range(n):
        piv = None
        for i in range(k, n):
            if M[i][k] != 0: piv = i; break
        if piv is None: return None
        M[k], M[piv] = M[piv], M[k]
        for i in range(n):
            if i != k and M[i][k] != 0:
                f = M[i][k] / M[k][k]
                M[i] = [a - f * b for a, b in zip(M[i], M[k])]
    return [M[i][n] / M[i][i] for i in range(n)]

def det_frac(Mx):
    n = len(Mx)
    M = [list(map(F, row)) for row in Mx]
    det = F(1)
    for k in range(n):
        piv = None
        for i in range(k, n):
            if M[i][k] != 0: piv = i; break
        if piv is None: return F(0)
        if piv != k:
            M[k], M[piv] = M[piv], M[k]; det = -det
        det *= M[k][k]
        for i in range(k + 1, n):
            if M[i][k] != 0:
                f = M[i][k] / M[k][k]
                M[i] = [a - f * b for a, b in zip(M[i], M[k])]
    return det

def charpoly_frac(Mx):
    n = len(Mx)
    Mk = [[F(0)] * n for _ in range(n)]
    c = [F(1)]
    for k in range(1, n + 1):
        AM = [[sum(Mx[i][l] * Mk[l][j] for l in range(n)) + (c[k - 1] if i == j else 0) for j in range(n)] for i in range(n)]
        Mk = AM
        tr = sum(Mx[i][l] * Mk[l][i] for i in range(n) for l in range(n))
        c.append(-tr / k)
    return c

def pointwise(A, eps, delta):
    """Three-valued exact decision at a point. Returns dict(stable in {True,False,None}, positive in {True,False,None}, ...)."""
    n = len(A)
    W = build_W(A, eps, delta)
    IW = [[(1 if i == j else 0) - W[i][j] for j in range(n)] for i in range(n)]
    x = solve_frac(IW, [F(1)] * n)
    if x is None:
        positive = None
    elif all(v > 0 for v in x): positive = True
    elif any(v == 0 for v in x): positive = None
    else: positive = False
    M = [[-1 + W[i][j] if i == j else W[i][j] for j in range(n)] for i in range(n)]
    a = charpoly_frac(M)
    H = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            idx = 2 * (i + 1) - (j + 1)
            if 0 <= idx <= n: H[i][j] = a[idx]
    minors = [det_frac([row[:k] for row in H[:k]]) for k in range(1, n + 1)]
    stable = True; detail = 'all Hurwitz minors positive'
    for k, m in enumerate(minors, 1):
        if m == 0: stable = None; detail = 'Hurwitz minor %d vanishes (UNDECIDED)' % k; break
        if m < 0: stable = False; detail = 'Hurwitz minor %d negative' % k; break
    q = delta * (1 - eps) / eps
    return {'stable': stable, 'positive': positive, 'detail': detail, 'x': x, 'minors': minors, 'q': q, 'legal': q > 1}

# ----------------------------------------------------------------- beta vs stable interval (bisection of stored minpoly)
def minpoly_eval(coeffs_hi_to_lo, x):
    r = F(0)
    for c in coeffs_hi_to_lo:
        r = r * x + c
    return r

def compare_beta_to_root(beta, interval, minpoly_hi_to_lo, max_steps=400):
    """Return (-1 if beta < root, +1 if beta > root, 0 if beta == root), steps used.
    root = unique root of minpoly in the isolating interval [lo, hi]."""
    lo, hi = F(interval[0]), F(interval[1])
    if lo > hi: lo, hi = hi, lo
    if lo == hi:
        assert minpoly_eval(minpoly_hi_to_lo, lo) == 0, "degenerate interval endpoint is not a root"
        return (-1 if beta < lo else (1 if beta > lo else 0)), 0
    flo = minpoly_eval(minpoly_hi_to_lo, lo); fhi = minpoly_eval(minpoly_hi_to_lo, hi)
    if flo == 0 and fhi == 0:
        raise ValueError("both endpoints are roots: interval not isolating")
    if flo == 0:
        root = lo; return (-1 if beta < root else (1 if beta > root else 0)), 0
    if fhi == 0:
        root = hi; return (-1 if beta < root else (1 if beta > root else 0)), 0
    assert (flo < 0) != (fhi < 0), "no sign change on isolating interval"
    steps = 0
    while lo <= beta <= hi:
        if minpoly_eval(minpoly_hi_to_lo, beta) == 0:
            return 0, steps
        mid = (lo + hi) / 2
        fm = minpoly_eval(minpoly_hi_to_lo, mid)
        if fm == 0:
            root = mid
            return (-1 if beta < root else (1 if beta > root else 0)), steps
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
        steps += 1
        if steps > max_steps: raise RuntimeError("bisection did not separate")
    return (-1 if beta < lo else 1), steps

def beta_in_interval(beta, iv):
    """beta strictly inside the stored stable interval iv = {'lo':..., 'hi':...}? returns (bool, info)"""
    info = {}
    lo = iv['lo']
    if lo['kind'] == 'beta_min':
        lo_ok = beta > F(lo['beta']); info['lo'] = ('rational', str(F(lo['beta'])), lo_ok)
    else:
        s, st = compare_beta_to_root(beta, lo['beta_interval'], lo['beta_minpoly_hi_to_lo'])
        lo_ok = s > 0; info['lo'] = ('root', lo['beta_interval'], st, lo_ok)
    hi = iv['hi']
    if hi['kind'] == 'beta_min':
        hi_ok = beta < F(hi['beta']); info['hi'] = ('rational', str(F(hi['beta'])), hi_ok)
    else:
        s, st = compare_beta_to_root(beta, hi['beta_interval'], hi['beta_minpoly_hi_to_lo'])
        hi_ok = s < 0; info['hi'] = ('root', hi['beta_interval'], st, hi_ok)
    return lo_ok and hi_ok, info
