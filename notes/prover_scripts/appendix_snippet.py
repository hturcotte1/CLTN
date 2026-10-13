import sympy as sp
z, x, y, b = sp.symbols('z x y b')
N = sp.Matrix([[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]])  # Geneson (19)
m = 6
V = sp.Matrix(m, m-1, lambda i, j: 1 if i == j else (-1 if i == m-1 else 0))      # basis of 1-perp
G = (V.T*V).inv() * V.T * N * V                                                  # similar to U^T N U  (C'(g))
print(sp.factor(G.charpoly(z).as_expr()))                                        # z*(6z^4+11z^3+16z^2+10z+3)/6
g = 6*y**4 + 11*y**3 + 16*y**2 + 10*y + 3
R = sp.Poly(sp.resultant(g, g.subs(y, x - y), y), x)                             # roots nu_i+nu_j, all ordered pairs
Q, r = sp.div(R, sp.Poly(16*g.subs(y, x/2), x)); assert r.is_zero                 # drop i = j
c, fl = sp.factor_list(Q.as_expr()); assert len(fl) == 1 and fl[0][1] == 2        # Q = c * S^2
S = fl[0][0]
sext = 797*b**6 - 6468*b**5 + 22544*b**4 - 44440*b**3 + 53280*b**2 - 38016*b + 13824
print(sp.cancel(sp.expand(S.subs(x, -2/b)*b**6) / sext))                         # 1
