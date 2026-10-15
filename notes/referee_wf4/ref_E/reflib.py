"""Independent exact helpers for refereeing Theorem E / Lemma B. No floating point in any decision."""
from fractions import Fraction as F
import sympy as sp

z, l, a = sp.symbols('z l a')

def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63; p = s[2:]; A = [[0]*n for _ in range(n)]
    bits = []
    for ch in p:
        v = ord(ch) - 63
        bits += [(v >> (5-k)) & 1 for k in range(6)]
    for i in range(n):
        for j in range(n):
            A[i][j] = bits[i*n+j]
    return A

def d6_encode(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(n+63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k+6]: v = 2*v + b
        s += chr(v+63)
    return s

def blowup(H, mu):
    """closed-twin blow-up: node c replaced by a bidirected clique of size mu[c]; classes listed in order."""
    cls = [c for c in range(len(H)) for _ in range(mu[c])]
    n = len(cls)
    return [[0 if i == j else (1 if cls[i] == cls[j] else H[cls[i]][cls[j]]) for j in range(n)] for i in range(n)], cls

def charpoly_int(M, var):
    return sp.Matrix(M).charpoly(var).as_expr()

def nchi_perp(A, var=l):
    """n * chi_{A|perp}(var) = 1^T adj(var I - A) 1 = det(var I - A + J) - det(var I - A)  (matrix det lemma)."""
    n = len(A)
    AmJ = [[A[i][j] - 1 for j in range(n)] for i in range(n)]
    return sp.expand(charpoly_int(AmJ, var) - charpoly_int(A, var))

def chi_psi_A(A, var=z):
    """chi(z)=det(zI+A), psi(z)=1^T adj(zI+A) 1 = det(zI+A+J)-det(zI+A)."""
    n = len(A)
    negA = [[-A[i][j] for j in range(n)] for i in range(n)]
    negAJ = [[-A[i][j]-1 for j in range(n)] for i in range(n)]
    c = sp.expand(charpoly_int(negA, var)); p = sp.expand(charpoly_int(negAJ, var) - c)
    return c, p

def det_frac(M):
    M = [[F(x) for x in r] for r in M]; n = len(M); d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c+1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n): M[r][k] -= f*M[c][k]
    return d

def hurwitz_minors(coeffs):
    """coeffs hi->lo (Fractions). returns list of leading Hurwitz minors Delta_1..Delta_m."""
    co = [F(x) for x in coeffs]; m = len(co) - 1
    if co[0] < 0: co = [-x for x in co]
    ai = lambda k: co[k] if 0 <= k <= m else F(0)
    Hm = [[ai(2*(i+1)-(j+1)) for j in range(m)] for i in range(m)]
    return [det_frac([r[:k] for r in Hm[:k]]) for k in range(1, m+1)], co

def all_roots_left_of(poly_expr, var, x0):
    """True iff every root of poly has Re < x0 (exact Hurwitz on p(x0 + y)); returns (verdict, minors). verdict None if a minor vanishes."""
    y = sp.Symbol('y_')
    P = sp.Poly(sp.expand(poly_expr.subs(var, x0 + y)), y)
    co = [F(int(sp.numer(c)), int(sp.denom(c))) for c in P.all_coeffs()]
    mins, co2 = hurwitz_minors(co)
    if any(c <= 0 for c in co2):
        return (None if any(c == 0 for c in co2) else False), mins
    if any(mm == 0 for mm in mins): return None, mins
    return all(mm > 0 for mm in mins), mins

# arithmetic in Q(z0), z0 = 1/2 + i sqrt(3)/2, z0^2 = z0 - 1 ; element (x, y) = x + y z0
def red(poly_expr, var):
    r = sp.rem(sp.Poly(poly_expr, var), sp.Poly(var**2 - var + 1, var))
    co = r.all_coeffs()
    if len(co) == 1: return (sp.simplify(co[0]), sp.Integer(0))
    if len(co) == 0: return (sp.Integer(0), sp.Integer(0))
    return (sp.simplify(co[1]), sp.simplify(co[0]))
def qmul(u, v):
    x1, y1 = u; x2, y2 = v
    # (x1 + y1 z)(x2 + y2 z) = x1x2 + (x1y2+x2y1) z + y1y2 (z - 1)
    return (sp.simplify(x1*x2 - y1*y2), sp.simplify(x1*y2 + x2*y1 + y1*y2))
def qinv(u):
    x, y = u; N = sp.simplify(x*x + x*y + y*y)
    return (sp.simplify((x + y)/N), sp.simplify(-y/N))
def qre(u):
    return sp.simplify(u[0] + u[1]/2)
