import sympy as sp
Nl = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]   # test_ctln.GENESON_N
m = 6; l = sp.symbols('l'); I, J = sp.eye(m), sp.ones(m, m)
A = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - Nl[j][i])         # A[i][j] = 1 iff i -> j
chiA, chiAJ = (l*I - A).det(), (l*I - A - J).det()
chiperp = sp.factor(sp.expand((chiA - chiAJ)/m))                          # Cor. I.1.1: m chi_{A|perp} = chi_A - chi_{A+J}
print(chiperp)                                                            # (l+1)(6l^4+13l^3+19l^2+13l+4)/6
alpha = max(sp.re(r) for r in sp.Poly(sp.expand(chiperp), l).nroots(n=30)); print(alpha, 1/(1 + alpha))
b = sp.Rational(54, 25); g = 1 - 1/b
w = (A.T + g*I).LUsolve(sp.ones(m, 1)); print(list(w), sum(w))            # I.4: w > 0, W1 > 1
h = sp.gcd(sp.Poly(chiA, l), sp.Poly(chiAJ, l)); print(h)                 # I.7: h = l + 1 (bidirected twins), harmless
