"""
tlib.py -- exact helpers for the general-n theory of full-support CTLN counterexamples.

Conventions: CLAUDE.md section 2.  A[i][j] = 1 iff i -> j.  N = J - I - A^T.
K = I + beta N + t J,  t = (1-eps)/eps,  beta = 1 + delta/eps,  gamma = 1 - 1/beta.
Legal iff t(beta-2) > 1.

Limit criterion (notes/lemma_I_spectral.md, Lemma I.3/I.4, PROVED there):
  limit-stable at beta  <=>  alpha(A|perp) < -gamma(beta)
  positivity at beta (every t)  <=>  w = (A^T + gamma I)^{-1} 1 > 0
  m chi_perp = chi_A - chi_{A+J}   (Corollary I.1.1)
No floating point in any decision made here.  Floats appear only in functions
whose name starts with 'heur_' and are used to steer searches.
"""
import sys, os
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
sys.path.insert(0, '/home/user/CLTN/n7')
import ctln
from fastregion import icharpoly

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A

def to_d6(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k+6]: v = 2 * v + b
        s += chr(v + 63)
    return s

def edges(A):
    n = len(A); return sum(A[i][j] for i in range(n) for j in range(n) if i != j)

def is_complete(A):
    n = len(A); return all(A[i][j] for i in range(n) for j in range(n) if i != j)

# ---------------------------------------------------------------- polynomials, hi -> lo integer lists

def chi(A):
    """det(z I - A), hi->lo integer coefficients."""
    return icharpoly([row[:] for row in A])

def chi_perp_times_m(A):
    """m * chi_{A|perp}(z), hi->lo, degree m-1, leading coeff m.  = chi_A - chi_{A+J}."""
    m = len(A)
    c1 = chi(A)
    AJ = [[A[i][j] + 1 for j in range(m)] for i in range(m)]
    c2 = chi(AJ)
    d = [x - y for x, y in zip(c1, c2)]
    assert d[0] == 0 and d[1] == m
    return d[1:]

def taylor_shift(p, c):
    """p hi->lo (Fractions ok); return coefficients of p(z + c), hi->lo."""
    p = [F(x) for x in p]
    n = len(p) - 1
    q = p[:]
    for i in range(n):
        for j in range(1, n - i + 1):
            q[j] += c * q[j - 1]
    return q

def primitive(p):
    from math import gcd
    from functools import reduce
    den = 1
    for x in p:
        den = den * F(x).denominator // gcd(den, F(x).denominator)
    ip = [int(F(x) * den) for x in p]
    g = reduce(gcd, [abs(x) for x in ip if x], 0) or 1
    ip = [x // g for x in ip]
    if ip[0] < 0: ip = [-x for x in ip]
    return ip

def alpha_less_than(A, c, P=None):
    """Exact three-valued test of alpha(A|perp) < c  (c rational).
    Roots of chi_perp all have Re < c  <=>  chi_perp(z + c) Hurwitz.
    Returns (True|False|None, detail).  False means alpha >= c (some root with Re >= c):
    by the Hurwitz theorem a polynomial with positive leading coeff is Hurwitz iff all
    Hurwitz minors > 0; if not Hurwitz then some root has Re >= 0.  ctln.hurwitz_stable returns
    None only when a minor vanishes -- then 'not Hurwitz' is still certain (a vanishing minor
    excludes Hurwitz), so None is mapped to False with detail 'boundary' (alpha >= c)."""
    P = P or chi_perp_times_m(A)
    q = primitive(taylor_shift(P, F(c)))
    v, det = ctln.hurwitz_stable([F(x) for x in q])
    if v is None:
        return False, 'not Hurwitz (vanishing minor/zero coeff): alpha >= c, boundary possible: ' + det
    return v, det

def w_vec(A, beta):
    """w = (A^T + gamma I)^{-1} 1 at rational beta (Lemma I.4).  None if singular."""
    beta = F(beta); g = 1 - 1 / beta; m = len(A)
    M = [[F(A[j][i]) + (g if i == j else 0) for j in range(m)] for i in range(m)]
    return ctln.solve_exact(M, [F(1)] * m)

def positive_at(A, beta):
    w = w_vec(A, beta)
    return w is not None and all(x > 0 for x in w)

def limit_ok_at(A, beta, P=None):
    """exact: beta > 2, limit-stable (alpha < -gamma) and positive at beta."""
    beta = F(beta)
    if beta <= 2: return False
    g = 1 - 1 / beta
    v, _ = alpha_less_than(A, -g, P)
    return bool(v) and positive_at(A, beta)

def sfp_witness(A, beta, ts=(10, 30, 100, 300, 1000, 3000, 10**4, 10**5, 10**6)):
    """Search t in ts for an exact legal stable positive full-support fixed point at this beta,
    certified by ctln.analyze_support.  Returns (t, eps, delta, record) or None."""
    beta = F(beta)
    n = len(A)
    for t in ts:
        t = F(t)
        if t * (beta - 2) <= 1: continue
        eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        assert ctln.legal(eps, delta)
        rec = ctln.analyze_support(A, list(range(n)), eps, delta)
        if rec['stable_fixed_point'] and rec['stable'] is True and not rec['target_free']:
            return t, eps, delta, rec
    return None

# ---------------------------------------------------------------- float heuristics (search only)
def heur_alpha(A):
    import numpy as np
    P = chi_perp_times_m(A)
    r = np.roots([float(x) for x in P])
    return max(r.real) if len(r) else -1.0

def heur_w(A, beta):
    import numpy as np
    g = 1 - 1 / beta
    M = np.array(A, dtype=float).T + g * np.eye(len(A))
    return np.linalg.solve(M, np.ones(len(A)))
