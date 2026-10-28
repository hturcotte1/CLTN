"""v_sympy for n = 8: independent exact verifier, sympy Rational/Poly arithmetic only, no floats.
Written from scratch (decoder, CTLN matrix, Hurwitz minors, Routh cross-check, Sturm counting);
the n = 7 verifier notes/verify_124/v_sympy/verify.py was used as the starting point.
"""
import json
from sympy import Rational, Matrix, eye, zeros, symbols, Poly, S, sturm, sign

lam = symbols('lam')
bt = symbols('bt')

def decode_d6(s):
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    need = n * n
    assert len(bits) >= need and len(bits) - need < 6, (s, len(bits), need)
    assert all(b == 0 for b in bits[need:]), ('nonzero padding', s)
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def encode_d6(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    out = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out += chr(63 + v)
    return out

def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def build_W(A, eps, delta):
    n = len(A)
    W = zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i == j:
                W[i, j] = S(0)
            elif A[j][i] == 1:          # j sends to i: weak inhibition
                W[i, j] = -1 + eps
            else:
                W[i, j] = -1 - delta
    return W

def hurwitz_minors(coeffs):
    """coeffs = [a0, a1, ..., an] of monic p(lam) = a0 lam^n + ... + an (a0 = 1)."""
    n = len(coeffs) - 1
    a = lambda k: coeffs[k] if 0 <= k <= n else S(0)
    H = Matrix(n, n, lambda i, j: a(2 * (j + 1) - (i + 1)))
    return [H[:k, :k].det(method='bareiss') for k in range(1, n + 1)]

def routh_first_column(coeffs):
    """Independent cross-check: exact Routh array first column, or None on a zero pivot."""
    n = len(coeffs) - 1
    r0 = [coeffs[k] for k in range(0, n + 1, 2)]
    r1 = [coeffs[k] for k in range(1, n + 1, 2)]
    width = len(r0)
    r1 = r1 + [S(0)] * (width - len(r1))
    rows = [r0, r1]
    for _ in range(n - 1):
        p, q = rows[-2], rows[-1]
        if q[0] == 0:
            return None
        new = [(q[0] * p[k + 1] - p[0] * q[k + 1]) / q[0] if k + 1 < width else S(0) for k in range(width)]
        rows.append(new)
    return [r[0] for r in rows]

def check(A, eps, delta, crosscheck=True):
    """Exact decision at one (eps, delta). stable in {True, False, 'UNDECIDED'}; positive in {True, False, None}."""
    n = len(A)
    eps, delta = Rational(eps), Rational(delta)
    res = {'eps': str(eps), 'delta': str(delta)}
    q = delta * (1 - eps) / eps
    res['q'] = str(q)
    res['legal'] = bool(eps > 0 and eps < 1 and delta > 0 and q > 1)
    W = build_W(A, eps, delta)
    M = eye(n) - W
    if crosscheck:
        At = Matrix(A).T
        J = Matrix(n, n, lambda i, j: 1)
        res['closed_form_ok'] = bool((-delta * eye(n) + (1 + delta) * J - (eps + delta) * At) - M == zeros(n, n))
    D = M.det(method='bareiss')
    res['det_I_minus_W'] = str(D)
    if D == 0:
        res['positive'] = None
        res['x'] = None
    else:
        x = M.LUsolve(Matrix([1] * n))
        assert M * x == Matrix([1] * n)
        res['x'] = [str(v) for v in x]
        res['positive'] = bool(all(v > 0 for v in x))
    Jac = -eye(n) + W
    cp = Jac.charpoly(lam)
    coeffs = [S(c) for c in cp.all_coeffs()]
    assert coeffs[0] == 1 and len(coeffs) == n + 1
    if crosscheck:
        P = Poly(cp.as_expr(), lam)
        res['charpoly_crosscheck'] = bool(all(P.eval(k) == (k * eye(n) - Jac).det(method='bareiss') for k in (-2, -1, 0, 1, 2)))
    res['charpoly_hi_to_lo'] = [str(c) for c in coeffs]
    mins = hurwitz_minors(coeffs)
    res['hurwitz_minors_sign'] = [(1 if m > 0 else (-1 if m < 0 else 0)) for m in mins]
    if any(m == 0 for m in mins):
        res['stable'] = 'UNDECIDED'
    elif all(m > 0 for m in mins):
        res['stable'] = True
    else:
        res['stable'] = False
    rf = routh_first_column(coeffs)
    if rf is None:
        res['routh'] = 'zero_pivot'
    else:
        res['routh'] = bool(all(v > 0 for v in rf))
    res['routh_agrees'] = (res['routh'] == res['stable']) if res['stable'] != 'UNDECIDED' else None
    return res

def point_from_beta(beta, t):
    """Limit-ray parametrisation: eps = 1/(t+1), delta = (beta-1)*eps."""
    beta = Rational(beta); t = Rational(t)
    eps = 1 / (t + 1)
    return eps, (beta - 1) * eps

# ---------- Sturm machinery (own implementation on top of sympy.sturm) ----------
def sturm_seq(coeffs_hi_to_lo):
    P = Poly([int(c) for c in coeffs_hi_to_lo], bt)
    return P, sturm(P)

def sign_variations(seq, x):
    vals = [p.eval(Rational(x)) for p in seq]
    vals = [v for v in vals if v != 0]
    return sum(1 for u, v in zip(vals, vals[1:]) if (u > 0) != (v > 0))

def count_roots_open(P, seq, a, b):
    """Number of distinct real roots in the open interval (a, b), a < b, assuming P(a) != 0 != P(b)."""
    a, b = Rational(a), Rational(b)
    assert P.eval(a) != 0 and P.eval(b) != 0
    return sign_variations(seq, a) - sign_variations(seq, b)

def refine_root(P, seq, a, b, iters=40):
    """Bisect [a,b] (containing exactly one root) until width <= 2^-iters * (b-a). Exact rationals."""
    a, b = Rational(a), Rational(b)
    for _ in range(iters):
        m = (a + b) / 2
        if P.eval(m) == 0:
            return m, m
        if count_roots_open(P, seq, a, m) == 1:
            b = m
        else:
            assert count_roots_open(P, seq, m, b) == 1
            a = m
    return a, b
