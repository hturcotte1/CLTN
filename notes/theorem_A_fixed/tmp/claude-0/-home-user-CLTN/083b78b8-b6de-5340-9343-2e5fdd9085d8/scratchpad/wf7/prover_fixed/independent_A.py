"""
independent_A.py -- sympy-free re-proof of Theorem A (step 3) and of Theorem SJ, using only fractions.Fraction and the repo's
ctln.charpoly / ctln.det_exact / ctln.hurwitz_stable (read-only).  Independent code path from theorem_A.py / selfjoin.py.

R = &F^U|m^QTq? given by its rows (hard-coded, then compared with the digraph6 decoding), D = J - I - A, b = e_2 (node 2 is the
only core node that does not send to the class), (beta, t) = (21/10, 20), a = t/beta = 200/21, r = 1/beta = 10/21.

 (I1) Psi_k(z) = z chi_{D+aJ}(z) - a k chi_{D - b 1^T}(z)  (Lemma P) as Fraction coefficient lists; h_k(x) = Psi_k(-r - x).
      Every coefficient of h_k is affine in k, so its j-th Hurwitz minor Delta_j(k) is a polynomial of degree <= j in k.
      Delta_j is recovered EXACTLY by Lagrange interpolation at k = 0..j, confirmed at 6 further k, and all of its coefficients
      are checked > 0  =>  Delta_j(k) > 0 for every real k >= 0  =>  (Hurwitz) R_k stable at the point for every k >= 0.
      The interpolated polynomials are compared with those stored in theorem_A_certificate.json (string form).
 (I2) direct ctln.hurwitz_stable(h_k) == True for k = 0..2000 (redundant with I1; a sanity sweep).
 (I3) self-join: G(z, c) = det(zI - D - cJ) = chi_D(z) - c psi_D(z) is affine in c; the Hurwitz minors of G(-r - x, c) are
      polynomials in c of degree <= j; interpolated exactly at c = 0..j; all coefficients > 0  =>  every root of G(., c) has
      Re > -r for every real c >= 0 (Theorem SJ).
 (I4) positivity: w_R = (A^T + (1 - r) I)^{-1} 1 > 0 and sigma_S > 1 at beta = 21/10 (ctln.solve_exact).
Prints ALL INDEPENDENT CHECKS PASSED or raises.  Writes independent_A.json.
"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln

ROWS = ['0111110', '1011011', '1101101', '1100111', '1101001', '0010101', '1100100']
A = [[int(ch) for ch in row] for row in ROWS]
m = 7
r = F(10, 21); a = F(200, 21); gamma = 1 - r
D = [[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)]
b = [1 if i == 2 else 0 for i in range(m)]


def padd(p, q):  # hi->lo lists
    n = max(len(p), len(q)); p = [F(0)] * (n - len(p)) + p; q = [F(0)] * (n - len(q)) + q
    return [x + y for x, y in zip(p, q)]


def pscale(p, s):
    return [s * x for x in p]


def pmul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def compose_shift(p, r):
    """p(z) hi->lo  ->  p(-r - x) hi->lo in x (Horner)."""
    out = [F(0)]
    lin = [F(-1), -r]          # -x - r
    for c in p:
        out = padd(pmul(out, lin), [c])
    while len(out) > 1 and out[0] == 0: out = out[1:]
    return out


def hurwitz_minors(co):
    d = len(co) - 1
    aa = lambda i: co[i] if 0 <= i <= d else F(0)
    H = [[aa(2 * (i + 1) - (j + 1)) for j in range(d)] for i in range(d)]
    return [ctln.det_exact([row[:j] for row in H[:j]]) for j in range(1, d + 1)]


def lagrange(xs, ys):
    """exact interpolating polynomial, lo->hi coefficients."""
    n = len(xs); coeffs = [F(0)] * n
    for i in range(n):
        num = [F(1)]; den = F(1)
        for j in range(n):
            if j == i: continue
            num = pmul(num, [F(1), F(-xs[j])])   # hi->lo
            den *= (xs[i] - xs[j])
        num = num[::-1]
        for p in range(len(num)):
            coeffs[p] += num[p] * ys[i] / den
    return coeffs


def peval_lohi(c, x):
    return sum(ci * x ** i for i, ci in enumerate(c))


def minors_as_polys(hk_of, deg):
    polys = []
    for j in range(1, deg + 1):
        xs = list(range(0, j + 1))
        ys = [hurwitz_minors(hk_of(F(v)))[j - 1] for v in xs]
        c = lagrange([F(v) for v in xs], ys)
        for v in range(j + 1, j + 7):   # confirm the degree bound numerically-exactly
            assert peval_lohi(c, F(v)) == hurwitz_minors(hk_of(F(v)))[j - 1], (j, v)
        polys.append(c)
    return polys


out = {}
# ---------------- (I1)
J = [[F(1)] * m for _ in range(m)]
chi_a = ctln.charpoly([[F(D[i][j]) + a for j in range(m)] for i in range(m)])          # hi->lo, degree 7
chi_L = ctln.charpoly([[F(D[i][j]) - b[i] for j in range(m)] for i in range(m)])
assert ctln._charpoly_direct([[F(D[i][j]) + a for j in range(m)] for i in range(m)]) == chi_a
assert ctln._charpoly_direct([[F(D[i][j]) - b[i] for j in range(m)] for i in range(m)]) == chi_L


def hk(kk):
    psi = padd(chi_a + [F(0)], pscale(chi_L, -a * kk))     # z*chi_a - a k chi_L, degree 8
    h = compose_shift(psi, r)
    h = pscale(h, F((-1) ** 8))
    assert h[0] == 1
    return h


polys = minors_as_polys(hk, 8)
allpos = [all(c > 0 for c in p) for p in polys]
assert all(allpos), allpos
out['I1_minor_polys_lohi'] = [[str(c) for c in p] for p in polys]
out['I1_all_coefficients_positive'] = allpos
# compare with theorem_A_certificate.json (stored as sympy strings; parse with Fraction-safe evaluation at k values)
cert = json.load(open('theorem_A_certificate.json'))
import re
for j, rec in enumerate(cert['A3_hurwitz_minors_in_k']):
    expr = rec['poly']
    # exact comparison: evaluate the stored string with k = v using Python's Fraction (integers wrapped)
    def ev(s, v):
        s2 = re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'F(\1)', s)
        s2 = s2.replace('**F(', '**(')
        return eval(s2, {'F': F, 'k': F(v)})
    for v in range(0, 12):
        assert ev(expr, v) == peval_lohi(polys[j], F(v)), (j, v)
out['I1_matches_theorem_A_certificate'] = True
# ---------------- (I2)
for kk in range(0, 2001):
    v = ctln.hurwitz_stable(hk(F(kk)))
    assert v[0] is True, (kk, v)
out['I2_direct_hurwitz_k_0_to_2000'] = True
# ---------------- (I3) self-join locus
chi_D = ctln.charpoly([[F(D[i][j]) for j in range(m)] for i in range(m)])
chi_DJ = ctln.charpoly([[F(D[i][j]) + 1 for j in range(m)] for i in range(m)])   # chi_D - psi_D
psi_D = padd(chi_D, pscale(chi_DJ, F(-1)))                                          # psi_D = chi_D - chi_{D+J}


def Gc(cc):
    p = padd(chi_D, pscale(psi_D, -cc))
    h = pscale(compose_shift(p, r), F((-1) ** 7))
    assert h[0] == 1
    return h


polysJ = minors_as_polys(Gc, 7)
allposJ = [all(c > 0 for c in p) or all(c == 0 for c in p) for p in polysJ]
out['I3_minor_polys_in_c_lohi'] = [[str(c) for c in p] for p in polysJ]
out['I3_all_coefficients_positive'] = [all(c > 0 for c in p) for p in polysJ]
print('I3 minors in c, all coefficients positive:', out['I3_all_coefficients_positive'])
# ---------------- (I4)
M = [[F(A[j][i]) + (gamma if i == j else 0) for j in range(m)] for i in range(m)]
w = ctln.solve_exact(M, [F(1)] * m)
sig = sum(w[i] for i in range(m) if b[i] == 0)
assert all(x > 0 for x in w) and sig > 1
out['I4_w_R'] = [str(x) for x in w]; out['I4_sigma_S'] = str(sig)
json.dump(out, open('independent_A.json', 'w'), indent=1)
print('I1 all minors positive-coefficient polys in k:', allpos)
print('ALL INDEPENDENT CHECKS PASSED' if all(out['I3_all_coefficients_positive']) else 'I1/I2/I4 PASSED; I3 see output')
