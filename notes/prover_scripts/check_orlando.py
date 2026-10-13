"""Sanity checks (NOT proofs): Orlando's formula, the kernel-vector argument for
Delta_{m-1}=0, and Geneson's eps->0 critical value via the compression G."""
import itertools, random
import sympy as sp

def hurwitz_minor(coeffs, k):
    m = len(coeffs) - 1
    a = lambda l: coeffs[l] if 0 <= l <= m else 0
    return sp.Matrix(k, k, lambda i, j: a(2 * (i + 1) - (j + 1))).det()

z = sp.symbols('z')
# 1. Orlando, symbolic roots, m = 2..5
for m in range(2, 6):
    r = sp.symbols(f'r1:{m+1}')
    p = sp.Poly(sp.prod([z - ri for ri in r]), z)
    co = p.all_coeffs()
    lhs = sp.expand(hurwitz_minor(co, m - 1))
    rhs = sp.expand((-1) ** (m * (m - 1) // 2) * sp.prod([r[i] + r[j] for i, j in itertools.combinations(range(m), 2)]))
    assert sp.expand(lhs - rhs) == 0, m
print("Orlando formula: symbolic identity holds for m = 2..5")
# m = 6,7: random Gaussian-integer roots in conjugate pairs / reals
random.seed(3)
for m in (6, 7):
    for _ in range(30):
        roots = []
        while len(roots) < m:
            if m - len(roots) >= 2 and random.random() < 0.5:
                a_, b_ = random.randint(-5, 5), random.randint(1, 5)
                roots += [a_ + b_ * sp.I, a_ - b_ * sp.I]
            else:
                roots.append(sp.Integer(random.randint(-6, 6)))
        p = sp.Poly(sp.expand(sp.prod([z - x for x in roots])), z)
        lhs = hurwitz_minor(p.all_coeffs(), m - 1)
        rhs = sp.expand((-1) ** (m * (m - 1) // 2) * sp.prod([roots[i] + roots[j] for i, j in itertools.combinations(range(m), 2)]))
        assert sp.expand(lhs - rhs) == 0
print("Orlando formula: 60 exact random instances m = 6,7 hold")

# 2. Kernel vector: p = (z^2+w^2) r(z); P(x)=x^m p(1/x)=(1+w^2x^2)R(x); V(x)=x R(-x)
x = sp.symbols('x')
for m in range(2, 8):
    for _ in range(20):
        w2 = sp.Rational(random.randint(1, 20), random.randint(1, 5))
        rr = sp.Poly([1] + [random.randint(-9, 9) for _ in range(m - 2)], z)
        p = sp.Poly(sp.expand((z**2 + w2) * rr.as_expr()), z)
        a = p.all_coeffs()                 # a[0] = 1 ... a[m]
        Px = sp.Poly(sum(a[k] * x**k for k in range(m + 1)), x)
        Rx, rem = sp.div(Px, sp.Poly(1 + w2 * x**2, x))
        assert rem.is_zero
        Vx = sp.Poly(x * Rx.as_expr().subs(x, -x), x)
        v = [Vx.coeff_monomial(x**j) for j in range(1, m)]
        assert any(c != 0 for c in v) and Vx.coeff_monomial(1) == 0 and Vx.degree() <= m - 1
        H = sp.Matrix(m - 1, m - 1, lambda i, j: a[2*(i+1)-(j+1)] if 0 <= 2*(i+1)-(j+1) <= m else 0)
        assert H * sp.Matrix(v) == sp.zeros(m - 1, 1)
        assert H.det() == 0
print("kernel-vector construction verified, m = 2..7")

# 3. Geneson graph: N from the paper (receiver-first nonedge matrix), compression G
N = sp.Matrix([[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]])
m = 6
V = sp.Matrix(m, m - 1, lambda i, j: 1 if i == j else (-1 if i == m - 1 else 0))
G = (V.T * V).inv() * V.T * N * V
chiG = sp.factor(sp.Matrix(G).charpoly(z).as_expr())
print("charpoly of G (compression of N to 1-perp):", chiG)
ev = sp.Poly(sp.Matrix(G).charpoly(z).as_expr(), z).nroots(n=30)
print("eigenvalues of G:", ev)
mre = min(sp.re(e) for e in ev)
print("min Re =", mre, "  -1/minRe =", -1 / mre)
sext = lambda b: 797*b**6 - 6468*b**5 + 22544*b**4 - 44440*b**3 + 53280*b**2 - 38016*b + 13824
print("sextic at -1/minRe:", sp.N(sext(-1 / mre), 20))
# exact: b_c is root of the sextic; check Re(nu) = -1/b_c via resultant: nu, nubar roots of chi; 2x = nu + nubar
