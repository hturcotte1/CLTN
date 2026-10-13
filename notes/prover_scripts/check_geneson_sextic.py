"""Exact: beta = -2/(nu_i+nu_j), i<j, nu_i roots of the quartic factor g of charpoly(G) (Geneson's graph),
are the roots of (a constant times) the sextic recorded in ANOMALIES.md / notes."""
import sympy as sp
x, y, b = sp.symbols('x y b')
g = 6*y**4 + 11*y**3 + 16*y**2 + 10*y + 3
R = sp.Poly(sp.resultant(g, g.subs(y, x - y), y), x)       # roots: nu_i + nu_j over all ordered pairs (i,j)
diag = sp.Poly(g.subs(y, x / 2) * 16, x)                     # roots: 2 nu_i
Q, rem = sp.div(R, diag); assert rem.is_zero                 # roots: nu_i + nu_j, i != j  (each twice)
c, fl = sp.factor_list(Q.as_expr())
print("Q = const * prod f^k with (deg, k):", [(sp.degree(f, x), k) for f, k in fl])
assert len(fl) == 1 and fl[0][1] == 2 and sp.degree(fl[0][0], x) == 6
S = fl[0][0]                                                  # roots nu_i + nu_j, i < j
Sb = sp.Poly(sp.expand(S.subs(x, -2 / b) * b**6), b)
sext = sp.Poly(797*b**6 - 6468*b**5 + 22544*b**4 - 44440*b**3 + 53280*b**2 - 38016*b + 13824, b)
ratio = sp.cancel(Sb.as_expr() / sext.as_expr())
print("S(-2/b) b^6 / sextic =", ratio)
assert ratio.is_number
print("sextic irreducible over Q:", sext.is_irreducible)
