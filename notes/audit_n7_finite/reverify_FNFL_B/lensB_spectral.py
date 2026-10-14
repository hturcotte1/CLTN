"""Lens B, part 1: adjacency consistency and exact spectral abscissa of A|perp."""
import sys, json
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')

cert = json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
g = cert['graph']
rows = g['adjacency_rows']
A = [[int(c) for c in r] for r in rows]
n = len(A)

def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    M = [[bits[i*n+j] for j in range(n)] for i in range(n)]
    pad = bits[n*n:]
    return n, M, pad

def d6_encode(M):
    n = len(M)
    bits = [M[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(63+n) + ''.join(chr(63 + int(''.join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))

n6, M6, pad = d6_decode(g['digraph6_as_enumerated'])
print("d6 n =", n6, "pad bits", pad)
print("adjacency_rows == digraph6 decode:", M6 == A)
print("re-encode:", d6_encode(A), "== cert:", d6_encode(A) == g['digraph6_as_enumerated'])
print("diag zero:", all(A[i][i] == 0 for i in range(n)), "edges:", sum(map(sum, A)), "cert edges:", g['edges'])
clique = all(A[i][j] and A[j][i] for i in range(n) for j in range(n) if i != j)
print("clique:", clique)

l = sp.symbols('l')
Am = sp.Matrix(A)
Jm = sp.ones(n, n)
chiA = (l*sp.eye(n) - Am).det(method='berkowitz')
chiAJ = (l*sp.eye(n) - (Am + Jm)).det(method='berkowitz')
diff = sp.expand(chiA - chiAJ)
chiperp = sp.Poly(sp.expand(diff / n), l)
print("chi_A - chi_{A+J} =", sp.factor(diff))
print("chi_perp monic? lc =", chiperp.LC(), " degree", chiperp.degree())

# independent route: explicit orthonormal-free compression. Use basis V of 1^perp (e_i - e_n);
# A|perp is similar to (V^+ P A V) with P=I-J/n: chi of the map x -> P A x on 1^perp.
V = sp.Matrix(n, n-1, lambda i, j: (1 if i == j else 0) - (1 if i == n-1 else 0))
P = sp.eye(n) - Jm / n
Mred = (V.T*V).inv() * V.T * P * Am * V   # matrix of P A restricted to 1^perp in basis V
chi2 = sp.Poly((l*sp.eye(n-1) - Mred).det(method='berkowitz'), l)
print("chi_perp via explicit compression equals (chi_A - chi_{A+J})/n:", sp.expand(chi2.as_expr() - chiperp.as_expr()) == 0)

fac = sp.factor_list(chiperp.as_expr())
print("factor over Q:", fac)
claimed = (l+1)*(l**2+l+1)*(7*l**3+14*l**2+14*l+5)
print("7*chi_perp == claimed product:", sp.expand(7*chiperp.as_expr() - claimed) == 0)

def hurwitz_minors(coeffs):
    # coeffs a0..ad with a0 leading (>0); Hurwitz matrix H[i][j] = a_{2(i+1)-(j+1)}
    d = len(coeffs) - 1
    a = lambda k: coeffs[k] if 0 <= k <= d else 0
    H = sp.Matrix(d, d, lambda i, j: a(2*(i+1) - (j+1)))
    return [H[:k, :k].det() for k in range(1, d+1)]

# Exact real-part analysis per irreducible factor.
alpha_candidates = []
x = sp.symbols('x')
for f, mult in fac[1]:
    fp = sp.Poly(f, l)
    print("\nfactor", f, "mult", mult, "deg", fp.degree())
    # (1) sympy exact isolation of all complex roots (rational boxes)
    ivs = fp.intervals(all=True, eps=F(1, 10**6))
    reals, cplx = ivs
    for (a, b), k in reals:
        print("  real root in [%s, %s]" % (a, b))
    for (lo, hi), k in cplx:
        print("  complex root in Re[%s,%s] x Im[%s,%s]" % (sp.re(lo), sp.re(hi), sp.im(lo), sp.im(hi)))
    # (2) exact real parts: Re z for roots z of f are roots of R(s) = res_y(f(y), f(s - y)) / 2 for the pair (z, zbar)
    y, s = sp.symbols('y s')
    R = sp.Poly(sp.resultant(sp.Poly(f.subs(l, y), y), sp.Poly(f.subs(l, s - y), y)), s)
    # (3) shift test: does f(x - 1/2) have all roots in open LHP (Hurwitz, strict)? or any on the line?
    g_sh = sp.Poly(sp.expand(f.subs(l, x - sp.Rational(1, 2))), x)
    cs = g_sh.all_coeffs()
    if cs[0] < 0: cs = [-c for c in cs]
    hm = hurwitz_minors(cs)
    print("  f(x-1/2) coeffs", cs, "Hurwitz minors", hm, "all coeff>0:", all(c > 0 for c in cs))
    # roots on Re = -1/2 line: f(-1/2 + i w) = 0 for real w -> gcd of real and imag parts
    w = sp.symbols('w', real=True)
    val = sp.expand(f.subs(l, -sp.Rational(1, 2) + sp.I*w))
    re_, im_ = sp.Poly(sp.re(val), w), sp.Poly(sp.im(val), w)
    gg = sp.gcd(re_, im_)
    print("  gcd(Re,Im) of f(-1/2+iw) in w:", gg.as_expr())
    print("  sum-of-pairs resultant R(s) factors:", sp.factor_list(R.as_expr()))

# Assemble exact verdict: per factor
print()
