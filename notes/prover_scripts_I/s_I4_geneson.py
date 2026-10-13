"""Sanity (not proof): Lemmas I.3-I.6 on Geneson's graph (test_ctln.py GENESON_N; A[b][a] = 1 - N[a][b])."""
import sympy as sp, mpmath as mp
mp.mp.dps = 50
Nl = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]
m = 6
A = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - Nl[j][i])
E = sum(A); print("|E| =", E, " e =", m*(m-1) - E, " density =", sp.Rational(E, m*(m-1)))
print("A =", A.tolist())
pairs0 = [(i, j) for i in range(m) for j in range(i+1, m) if A[i, j] == 0 and A[j, i] == 0]
pairs1 = [(i, j) for i in range(m) for j in range(i+1, m) if A[i, j] + A[j, i] == 1]
print("non-adjacent pairs (0-based):", pairs0, "  one-way pairs:", pairs1)
lam = sp.symbols('lam')
chiA = (lam*sp.eye(m) - A).det(); chiAJ = (lam*sp.eye(m) - A - sp.ones(m, m)).det()
chiperp = sp.factor(sp.expand((chiA - chiAJ)/m))
print("chi_{A|perp}(lam) =", chiperp)
roots = sp.Poly(sp.expand(chiperp), lam).nroots(n=30)
alpha = max(sp.re(r) for r in roots)
print("spec(A|perp) =", [sp.N(r, 12) for r in roots]); print("alpha(A|perp) =", sp.N(alpha, 15))
bc = 1/(1 + alpha); print("beta_c = 1/(1+alpha) =", sp.N(bc, 15), "  q_c =", sp.N(bc - 1, 15))
sext = lambda b: 797*b**6 - 6468*b**5 + 22544*b**4 - 44440*b**3 + 53280*b**2 - 38016*b + 13824
print("lemmas.md sextic at beta_c:", sp.N(sext(bc), 10))
# trace bound
print("trace(A|perp) =", -sp.Rational(E, m), " mean Re =", sp.N(-sp.Rational(E, m*(m-1)), 10))
# positivity at beta = 54/25
b = sp.Rational(54, 25); g = 1 - 1/b
Mg = A.T + g*sp.eye(m); w = Mg.LUsolve(sp.ones(m, 1)); W1 = sum(w)
print("gamma =", g, " w =", list(w), " W1 =", W1, " (-1)^(m-1) det M =", sp.sign((-1)**(m-1)*Mg.det()))
v = [21597025, 34234375, 34234375, 46065775, 5850625, 67159525]
print("w proportional to paper's v:", all(w[i]*v[0] == w[0]*v[i] for i in range(m)))
print("x_inf/theta = w/W1 =", [sp.nsimplify(x) for x in (w/W1)])
# symmetric part of A|perp: lambda_max, lambda_min (Helmert basis, floats for sanity)
U = sp.zeros(m, m-1)
for k in range(1, m):
    for i in range(k): U[i, k-1] = 1
    U[k, k-1] = -k
    U[:, k-1] = U[:, k-1]/sp.sqrt(k*(k+1))
S = (U.T*(A + A.T)*U)/2
ev = sorted(mp.eig(mp.matrix(S.evalf(40).tolist()), left=False, right=False), key=lambda z: mp.re(z))
print("spec sym(A|perp):", [mp.nstr(mp.re(z), 8) for z in ev])
