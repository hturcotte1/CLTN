"""Independent exact verifier (label v_roots).
No project code imported. All decisions in Fraction / integer arithmetic.
Stability: count roots of p(mu)=det(mu I - M), M = -I + W, in Re(mu) >= 0 by
  (R) exact Routh array (zero pivot -> UNDECIDED),
  (C) Cauchy index of p(iy) via an own Sturm chain (gcd(Re,Im) nonconstant -> UNDECIDED),
  (S) sympy exact complex root isolation (Poly.intervals(all=True)), refined until
      every isolating rectangle lies strictly on one side of Re = 0.
"""
import json, sys, itertools
from fractions import Fraction as F
import sympy as sp

# ---------------- digraph6 -----------------
def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    assert 0 < n < 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    need = n * n
    assert len(bits) >= need and len(bits) - need < 6 and not any(bits[need:]), s
    A = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        assert A[i][i] == 0
    return A

def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def canon(A):
    n = len(A)
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best

# ---------------- matrices -----------------
def weight_matrix(A, eps, delta):
    n = len(A)
    W = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                W[i][j] = F(0)
            elif A[j][i] == 1:          # j sends to i: weak
                W[i][j] = -1 + eps
            else:
                W[i][j] = -1 - delta
    return W

def solve(Mx, b):
    """Gaussian elimination in Fractions. Returns (det, x) ; x None if singular."""
    n = len(Mx)
    M = [row[:] + [b[i]] for i, row in enumerate(Mx)]
    det = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return F(0), None
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    x = [F(0)] * n
    for r in range(n - 1, -1, -1):
        s = M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))
        x[r] = s / M[r][r]
    return det, x

def det(Mx):
    return solve(Mx, [F(0)] * len(Mx))[0]

def charpoly(M):
    """Faddeev-LeVerrier: returns [c0=1, c1, ..., cn] with p(mu)=sum c_k mu^(n-k) = det(mu I - M)."""
    n = len(M)
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    def mul(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    c = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = mul(M, Mk)
        Mk = [[Mk[i][j] + c[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AM = mul(M, Mk)
        ck = -sum(AM[i][i] for i in range(n)) / k
        c.append(ck)
    return c

def peval(c, x):
    v = F(0)
    for a in c:
        v = v * x + a
    return v

# ---------------- (R) Routh -----------------
def routh_rhp(c):
    """Number of roots with Re>0 from first-column sign changes; None if a zero pivot occurs."""
    n = len(c) - 1
    r0 = [c[i] for i in range(0, n + 1, 2)]
    r1 = [c[i] for i in range(1, n + 1, 2)]
    rows = [r0, r1]
    width = len(r0)
    rows = [r + [F(0)] * (width - len(r)) for r in rows]
    for _ in range(n - 1):
        a, b = rows[-2], rows[-1]
        if b[0] == 0:
            return None, [str(r[0]) for r in rows]
        new = [(b[0] * a[k + 1] - a[0] * b[k + 1]) / b[0] for k in range(width - 1)] + [F(0)]
        rows.append(new)
    first = [r[0] for r in rows]
    if any(x == 0 for x in first):
        return None, [str(x) for x in first]
    changes = sum(1 for i in range(len(first) - 1) if (first[i] > 0) != (first[i + 1] > 0))
    return changes, None

# ---------------- (C) Cauchy index -----------------
def pdeg(p):
    while p and p[-1] == 0:
        p = p[:-1]
    return p
def prem(a, b):  # polys low->high, remainder of a by b
    a = pdeg(a[:]); b = pdeg(b[:])
    assert b
    while a and len(a) >= len(b):
        f = a[-1] / b[-1]; s = len(a) - len(b)
        for i in range(len(b)):
            a[i + s] -= f * b[i]
        assert a[-1] == 0
        a = pdeg(a[:-1])
    return a
def sign_inf(p, plus):
    d = len(p) - 1; lc = p[-1]
    s = 1 if lc > 0 else -1
    if not plus and d % 2 == 1:
        s = -s
    return s
def var(seq):
    s = [x for x in seq if x != 0]
    return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
def cauchy_index_line(Q, P):
    """Ind_{-inf}^{+inf} Q/P via Sturm chain P, Q mod P, ...; returns (index, gcd_degree)."""
    Q = prem(Q, P)
    chain = [pdeg(P), Q]
    while chain[-1]:
        r = prem(chain[-2], chain[-1])
        chain.append([-x for x in r])
    chain = [p for p in chain if p]
    g = chain[-1]
    V = lambda plus: var([sign_inf(p, plus) for p in chain])
    return V(False) - V(True), len(g) - 1

def cauchy_rhp(c):
    """c high->low. p(iy) = R(y) + i I(y). Delta arg / pi = n_L - n_R.
    Delta arg/pi = -Ind(I/R) + [lim_{+inf} - lim_{-inf}] of arctan(I/R)/pi."""
    n = len(c) - 1
    low = list(reversed(c))  # low[k] coefficient of mu^k
    R = [F(0)] * (n + 1); I = [F(0)] * (n + 1)
    for k, a in enumerate(low):
        # (iy)^k = i^k y^k
        m = k % 4
        if m == 0: R[k] += a
        elif m == 1: I[k] += a
        elif m == 2: R[k] -= a
        else: I[k] -= a
    R = pdeg(R); I = pdeg(I)
    if not R or not I:
        return None, 'R or I identically zero'
    ind, gdeg = cauchy_index_line(I, R)
    if gdeg > 0:
        return None, 'gcd(Re p(iy), Im p(iy)) nonconstant: root on imaginary axis'
    dR, dI = len(R) - 1, len(I) - 1
    if dI < dR:
        bnd = F(0)
    else:  # I/R -> +-inf at the ends
        sp_ = sign_inf(I, True) * sign_inf(R, True)
        sm_ = sign_inf(I, False) * sign_inf(R, False)
        bnd = F(sp_ - sm_, 2)
    diff = -ind + bnd
    assert diff.denominator == 1
    diff = int(diff)
    assert (n - diff) % 2 == 0
    return (n - diff) // 2, None

# ---------------- (S) sympy root isolation -----------------
mu = sp.Symbol('mu')
def sympy_rhp(c):
    """Square-free decomposition, then exact isolation of each square-free factor;
    roots counted with multiplicity."""
    P = sp.Poly([sp.Rational(x.numerator, x.denominator) for x in c], mu, domain='QQ')
    _, facs = P.sqf_list()
    assert sum(f.degree() * m for f, m in facs) == P.degree()
    total = 0
    for f, mult in facs:
        if f.degree() == 0:
            continue
        eps_ = sp.Rational(1, 10**6); done = False
        for it in range(12):
            real, cplx = f.intervals(all=True, eps=eps_)
            cnt = sum(m for _, m in real) + sum(m for _, m in cplx)
            assert cnt == f.degree(), (cnt, f.degree())
            rhp = 0; ambiguous = False
            for (a, b), m in real:
                a = sp.Rational(a); b = sp.Rational(b)
                if a > 0: rhp += m
                elif b < 0: pass
                else: ambiguous = True
            for (z1, z2), m in cplx:
                x1 = sp.re(z1); x2 = sp.re(z2)
                assert x1.is_Rational and x2.is_Rational and x1 <= x2
                if x1 > 0: rhp += m
                elif x2 < 0: pass
                else: ambiguous = True
            if not ambiguous:
                total += rhp * mult; done = True; break
            eps_ = eps_ / 10**4
        if not done:
            return None, 'isolation did not separate from imaginary axis'
    note = None if all(m == 1 for _, m in facs) else 'repeated roots: multiplicities %s' % sorted(set(m for _, m in facs))
    return total, note

# ---------------- a point check -----------------
def check_point(A, eps, delta, use_sympy=True):
    n = len(A)
    out = {'eps': str(eps), 'delta': str(delta)}
    q = delta * (1 - eps) / eps
    out['q'] = str(q)
    out['legal'] = bool(0 < eps < 1 and delta > 0 and q > 1)
    W = weight_matrix(A, eps, delta)
    IW = [[F(int(i == j)) - W[i][j] for j in range(n)] for i in range(n)]
    # closed-form cross-check with A read transposed (edge j -> i)
    for i in range(n):
        for j in range(n):
            cf = -delta * (i == j) + (1 + delta) - (eps + delta) * (A[j][i] if i != j else 0)
            assert cf == IW[i][j]
    D, x = solve(IW, [F(1)] * n)
    out['det_IminusW'] = str(D)
    out['nondegenerate'] = D != 0
    if x is None:
        out['positive'] = None
    else:
        # verify the solution exactly
        assert all(sum(IW[i][j] * x[j] for j in range(n)) == 1 for i in range(n))
        out['x'] = [str(v) for v in x]
        out['positive'] = all(v > 0 for v in x)
    M = [[W[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)]
    c = charpoly(M)
    # cross-check charpoly against direct determinants at 3 rational points
    for r in (F(0), F(1, 3), F(-5, 2), F(7)):
        assert peval(c, r) == det([[r * (i == j) - M[i][j] for j in range(n)] for i in range(n)])
    out['charpoly'] = [str(v) for v in c]
    rR, detR = routh_rhp(c)
    rC, detC = cauchy_rhp(c)
    rS, detS = sympy_rhp(c) if use_sympy else (None, 'skipped')
    out['rhp_routh'] = rR; out['rhp_cauchy'] = rC; out['rhp_sympy'] = rS
    if detR: out['routh_undecided_first_col'] = detR
    if detC: out['cauchy_note'] = detC
    if detS: out['sympy_note'] = detS
    decided = [v for v in (rR, rC, rS) if v is not None]
    if not decided:
        out['stable'] = None
    elif len(set(decided)) > 1:
        out['stable'] = 'DISAGREE'
    else:
        # stable iff no roots with Re >= 0: Routh decides only when no zero pivot (so no
        # imaginary-axis roots); Cauchy returns None if imaginary-axis roots exist; sympy
        # flags ambiguity.  Require all three decided for a 'True'.
        out['stable'] = (decided[0] == 0) if len(decided) == 3 else None
    out['methods_agree'] = (len(set(decided)) <= 1)
    return out
