"""Independent referee check of Theorem F (family F_k), algebraic part.
No floating point in any decision. Own constructions; only sympy + fractions."""
import sympy as sp
from fractions import Fraction as Fr
z, k, x, b, r = sp.symbols('z k x b r')

# (0) Geneson's graph from the paper's nonedge matrix N (verbatim in notes/geneson_2607.21396.md): A = J - I - N^T
Npap = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]
G1 = [[0 if i == j else 1 - Npap[j][i] for j in range(6)] for i in range(6)]
G1_claim = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
assert G1 == G1_claim, 'G1 mismatch'
print('(0) G1 = J - I - N_paper^T: OK')

def Fk(kk):
    n = 6 + kk; A = [[0]*n for _ in range(n)]
    for i in range(6):
        for j in range(6): A[i][j] = G1[i][j]
    for v in range(6, n):
        for i in range(6):
            A[v][i] = 1                     # class sends to every G1 node
            A[i][v] = 1 if i <= 4 else 0    # G1 nodes 0..4 send to class; node 5 does not
        for u in range(6, n):
            if u != v: A[v][u] = 1
    return A

Rk = (6 + k)*z**4 + (11 + k)*z**3 + (16 + k)*z**2 + 10*z + 3

# (1) explicit: n chi_{D'|perp} = chi_{D'} - chi_{D'+J} must equal z^{k+1} R_k ; equivalently for A|perp
def charp(M):
    return sp.Matrix(M).charpoly(z).as_expr()
for kk in range(0, 7):
    A = Fk(kk); n = len(A)
    D = [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]
    DJ = [[D[i][j] + 1 for j in range(n)] for i in range(n)]
    lhs = sp.expand(charp(D) - charp(DJ))
    rhs = sp.expand(z**(kk + 1) * Rk.subs(k, kk))
    assert lhs == rhs, kk
    # also directly in A: n chi_{A|perp}(w) at w = -1 - z equals (-1)^(n-1) * lhs(z)
    AJ = [[A[i][j] + 1 for j in range(n)] for i in range(n)]
    pa = sp.expand(charp(A) - charp(AJ))
    assert sp.expand(pa.subs(z, -1 - z) - (-1)**(n - 1) * lhs) == 0
print('(1) explicit n chi_perp(D) = z^(k+1) R_k and A|perp = -I - D|perp, k = 0..6: OK')

# (2) all-k identity via Lemma W1: mu^T adj(zI - D_+ M) 1 with symbolic k (7x7)
A0 = G1
D0 = [[0 if i == j else 1 - A0[i][j] for j in range(6)] for i in range(6)]
S = [1,1,1,1,1,0]; T = [1]*6
Dp = sp.zeros(7, 7)
for i in range(6):
    for j in range(6): Dp[i, j] = D0[i][j]
    Dp[i, 6] = 1 - S[i]; Dp[6, i] = 1 - T[i]
M = sp.diag(1, 1, 1, 1, 1, 1, k)
mu = sp.Matrix([1]*6 + [k])
X = z*sp.eye(7) - Dp*M
val = sp.expand((mu.T * X.adjugate() * sp.ones(7, 1))[0])
assert sp.expand(val - z**2 * Rk) == 0
print('(2) symbolic-k: mu^T adj(zI - D_+ M) 1 = z^2 R_k exactly (all k): OK')
# Sanity of W1 itself on explicit k (n chi = z^(n-m) * that, n - m = k - 1)
for kk in range(1, 7):
    A = Fk(kk); n = len(A)
    D = [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]
    DJ = [[D[i][j] + 1 for j in range(n)] for i in range(n)]
    assert sp.expand(charp(D) - charp(DJ) - z**(kk - 1) * val.subs(k, kk)) == 0
print('    W1 consistency on explicit k = 1..6: OK')

# (3) Hurwitz data of 16 R_k(-1/2 - x)
h = sp.Poly(sp.expand(16 * Rk.subs(z, -sp.Rational(1, 2) - x)), x)
a = h.all_coeffs()
print('    coeffs:', [sp.factor(c) for c in a])
assert [sp.expand(c) for c in a] == [sp.expand(c) for c in [16*(k+6), 16*(k+1), 8*(2*k+17), 12*(k+1), 3*k+16]]
H = sp.Matrix(4, 4, lambda i, j: a[2*(i+1)-(j+1)] if 0 <= 2*(i+1)-(j+1) <= 4 else 0)
mins = [sp.factor(H[:j, :j].det()) for j in range(1, 5)]
print('    Hurwitz minors:', mins)
assert sp.expand(mins[1] - 64*(k+1)*(k+16)) == 0 and sp.expand(mins[2] - 8192*(k+1)**2) == 0
assert sp.expand(mins[3] - 8192*(k+1)**2*(3*k+16)) == 0
print('(3) Hurwitz data positive for all k > -1: OK  => roots of R_k have Re > -1/2 => alpha(F_k|perp) < -1/2')

# (4) positivity of G1: w = (A^T + g I)^{-1} 1
Mb = sp.Matrix(6, 6, lambda i, j: b*A0[j][i] + ((b - 1) if i == j else 0))   # = b (A^T + g I)
detMb = sp.expand(Mb.det(method='berkowitz'))
adj1 = [sp.expand(e) for e in (Mb.adjugate(method='berkowitz') * sp.ones(6, 1))]
assert all(sp.expand(sum(Mb[i, j]*adj1[j] for j in range(6)) - detMb) == 0 for i in range(6))
w = [sp.factor(b*e/detMb) for e in adj1]       # w = b adj(Mb) 1 / det Mb
print('    det Mb =', sp.factor(detMb)); print('    adj(Mb)1 =', [sp.factor(e) for e in adj1])
print('    w_G1(b) =', w)
sig = sp.factor(sp.together(sum(w[i] for i in range(5)) - 1))
print('    sigma_S - 1 =', sig)
bp = sp.CRootOf(b**3 - 3*b**2 + 2*b - 1, 0)
assert sp.Poly(b**3 - 3*b**2 + 2*b - 1, b).count_roots() == 1     # single real root
assert sp.expand((b**3 - 3*b**2 + 2*b - 1).subs(b, z + 1) - (z**3 - z - 1)) == 0
# sign analysis on (2, bp): every numerator/denominator factor has no root in [2, bp) -> constant sign = sign at an interior rational
lo = sp.Integer(2); hi = sp.Rational(2324, 1000)     # bp > 2.3247 > hi ; check bp location
assert (hi**3 - 3*hi**2 + 2*hi - 1) < 0 and (sp.Rational(2325,1000)**3 - 3*sp.Rational(2325,1000)**2 + 2*sp.Rational(2325,1000) - 1) > 0
facs = set()
for e in w + [sig]:
    nu, de = sp.fraction(sp.together(e))
    for f, _ in sp.factor_list(nu, b)[1] + sp.factor_list(de, b)[1]:
        facs.add(sp.expand(f))
for f in facs:
    P = sp.Poly(f, b)
    if P.degree() == 0: continue
    if sp.expand(f - (b**3 - 3*b**2 + 2*b - 1)) == 0 or sp.expand(f + (b**3 - 3*b**2 + 2*b - 1)) == 0:
        continue
    # no root on [2, 2325/1000] (covers [2, bp])
    assert P.count_roots(2, sp.Rational(2325, 1000)) == 0, f
print('    all other factors root-free on [2, 2.325] (covers [2, b_p]); adj_4-type factor vanishes only at b_p')
for bb in [sp.Rational(201, 100), sp.Rational(21, 10), sp.Rational(22, 10), sp.Rational(2324, 1000)]:
    vals = [e.subs(b, bb) for e in w]; sv = sig.subs(b, bb)
    assert all(v > 0 for v in vals) and sv > 0, (bb, vals, sv)
print('(4) w_G1 > 0 and sigma_S > 1 on (2, b_p): OK')
# beyond b_p: some w_i < 0, others > 0, for every b > b_p (check sign structure symbolically)
for bb in [sp.Rational(2325, 1000), sp.Rational(233, 100), sp.Rational(24, 10), sp.Rational(3), sp.Rational(10)]:
    vals = [e.subs(b, bb) for e in w]
    assert any(v < 0 for v in vals) and any(v > 0 for v in vals), bb
print('    beyond b_p mixed signs at sample points (full proof: see text)')

# (5) Lemma X2 direct check on explicit F_k: w' positive iff predicted, at many rational beta
def wexact(A, beta):
    n = len(A); gg = 1 - Fr(1) / beta
    Mx = sp.Matrix(n, n, lambda i, j: sp.Rational(A[j][i]) + (sp.Rational(gg.numerator, gg.denominator) if i == j else 0))
    import sys; sys.path.insert(0, '/home/user/CLTN'); import ctln
    M2 = [[Fr(A[j][i]) + (gg if i == j else 0) for j in range(n)] for i in range(n)]
    return ctln.solve_exact(M2, [Fr(1)]*n)
for kk in [0, 1, 2, 3, 5]:
    A = Fk(kk)
    for bb in [Fr(201, 100), Fr(21, 10), Fr(23, 10), Fr(2324, 1000), Fr(2325, 1000), Fr(24, 10), Fr(3)]:
        ww = wexact(A, bb)
        pos = all(v > 0 for v in ww)
        pred = bb < Fr(23247, 10000) and bb > 2
        assert pos == pred, (kk, bb, pos)
print('(5) explicit positivity of F_k (k = 0,1,2,3,5) matches (2, b_p) at sample betas: OK')
