"""
classext.py -- exact algebra of the CLIQUE-CLASS EXTENSION  G[S,T;k].

G: m nodes, adjacency A (A[i][j]=1 iff i->j).  Add k new nodes V, bidirected among themselves;
i -> v for every v in V iff S[i] = 1;  v -> i for every v in V iff T[i] = 1.   n = m + k.
Non-edge digraph D := J - I - A  (D[i][j] = 1 iff i != j and i -/-> j);  b := 1 - S,  c := 1 - T.

Lemma X1 (proved in notes, checked here exactly):
    n * chi_{D'|perp}(z) = z^(k-1) * Phi_k(z),   Phi_k(z) = z*psi_D(z) + k*chi_L(z),
    psi_D(z) = 1^T adj(zI - D) 1 = m chi_{D|perp}(z),   L := D - b 1^T - 1 c^T.
and A'|perp = -I - D'|perp, so  alpha(A'|perp) = -1 - min Re{ roots of Phi_k, and 0 if k >= 2 }.
All polynomials are integer polynomials, hi -> lo lists.  No floating point.
"""
import sys
from fractions import Fraction as F
import sympy as sp
import tlib
from fastregion import icharpoly

z, kk = sp.symbols('z k')

def blow(A, S, T, k):
    m = len(A); n = m + k
    B = [[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    for i in range(m):
        for v in range(m, n):
            B[i][v] = S[i]; B[v][i] = T[i]
    for a in range(m, n):
        for b in range(m, n):
            if a != b: B[a][b] = 1
    return B

def nonedge_digraph(A):
    m = len(A); return [[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)]

def chi_poly(M):
    """sympy Poly det(zI - M) for an integer matrix M."""
    return sp.Poly(icharpoly([r[:] for r in M]), z)

def psi_poly(M):
    """1^T adj(zI - M) 1 = chi_M - chi_{M - J}  (T2 with tau = -1)."""
    m = len(M)
    MJ = [[M[i][j] - 1 for j in range(m)] for i in range(m)]   # M - J: det(zI - M + J) = chi_M + psi
    # det(zI - (M - J)) = det(zI - M + J) = chi_M(z) + 1^T adj(zI - M) 1
    return chi_poly(MJ) - chi_poly(M)

def L_matrix(A, S, T):
    D = nonedge_digraph(A); m = len(A)
    return [[D[i][j] - (1 - S[i]) - (1 - T[j]) for j in range(m)] for i in range(m)]

def Phi(A, S, T):
    """returns (Q, P) with Phi_k = Q + k P, Q = z psi_D, P = chi_L (sympy Polys in z)."""
    D = nonedge_digraph(A)
    Q = sp.Poly(z, z) * psi_poly(D)
    P = chi_poly(L_matrix(A, S, T))
    return Q, P

def check_X1(A, S, T, ks=range(1, 9)):
    """exact check of Lemma X1 against the explicit graphs."""
    Q, P = Phi(A, S, T)
    for k in ks:
        B = blow(A, S, T, k); n = len(B)
        Dp = nonedge_digraph(B)
        lhs = psi_poly(Dp)                         # = n chi_{D'|perp}
        rhs = sp.Poly(z**(k - 1), z) * (Q + k * P)
        if lhs != rhs:
            return False, k
    return True, None

def hurwitz_minors_poly(coeffs):
    """coeffs: list of sympy expressions (hi->lo) of a polynomial in x; returns Hurwitz minors Delta_1..Delta_d."""
    d = len(coeffs) - 1
    a = lambda i: coeffs[i] if 0 <= i <= d else 0
    H = sp.Matrix(d, d, lambda i, j: a(2 * (i + 1) - (j + 1)))
    return [sp.expand(H[:j, :j].det(method='berkowitz')) for j in range(1, d + 1)]

def shifted_family(A, S, T, c):
    """h_k(x) := (-1)^deg * Phi_k(-x + c)... we want roots zeta of Phi_k with Re zeta > c  <=>  x := c - zeta has Re x < 0.
    Returns coefficient list (hi->lo, in x) with entries polynomials in k, leading coefficient positive."""
    Q, P = Phi(A, S, T)
    x = sp.symbols('x')
    expr = sp.expand((Q.as_expr() + kk * P.as_expr()).subs(z, c - x))
    poly = sp.Poly(expr, x)
    co = [sp.expand(cf) for cf in poly.all_coeffs()]
    lead = co[0]
    if sp.Poly(lead, kk).LC() < 0:
        co = [-cf for cf in co]
    return co
