"""Sanity (not proof): Lemma I.1 and Corollary I.1.1, symbolic for m=2,3 and exact random for m=4,5."""
import sympy as sp, random
random.seed(11)
def orth_complement(u):
    m = len(u)
    # basis of u^perp, then Gram-Schmidt (exact, with sqrt)
    M = sp.Matrix([list(u)])
    ns = M.nullspace()
    on = sp.GramSchmidt(ns, orthonormal=True)
    return sp.Matrix.hstack(*on)
lam, tau = sp.symbols('lam tau')
cnt = 0
for m in [2, 3]:
    X = sp.Matrix(m, m, lambda i, j: sp.Symbol(f'x{i}{j}'))
    for u in [sp.Matrix([1]*m)/sp.sqrt(m), sp.Matrix([3, 4] + [0]*(m-2))/5, sp.Matrix([1, 2] + [2]*(m-2))/sp.sqrt(1 + 4*(m-1))]:
        U = orth_complement(u)
        assert sp.simplify(U.T*U - sp.eye(m-1)) == sp.zeros(m-1, m-1)
        lhs = sp.expand(sp.simplify((U.T*X*U).det()))
        rhs = sp.expand((u.T*X.adjugate()*u)[0])
        assert sp.simplify(lhs - rhs) == 0, (m, u); cnt += 1
for m in [4, 5]:
    for _ in range(5):
        X = sp.Matrix(m, m, lambda i, j: sp.Rational(random.randint(-9, 9), random.randint(1, 4)))
        u = sp.Matrix([1]*m)/sp.sqrt(m); U = orth_complement(u)
        lhs = sp.expand(sp.simplify((lam*sp.eye(m-1) - U.T*X*U).det()))
        rhs = sp.expand((u.T*(lam*sp.eye(m)-X).adjugate()*u)[0])
        assert sp.simplify(lhs - rhs) == 0; cnt += 1
        # Corollary I.1.1
        J = sp.ones(m, m)
        c1 = sp.expand((lam*sp.eye(m) - X - tau*J).det())
        c2 = sp.expand((lam*sp.eye(m) - X).det() - tau*m*lhs)
        assert sp.simplify(c1 - c2) == 0; cnt += 1
print("I.1 / I.1.1 sanity checks passed:", cnt)
