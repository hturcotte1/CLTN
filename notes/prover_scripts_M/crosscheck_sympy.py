"""EXACT, independent of ctln: sympy rational arithmetic for K = I + beta N + t J (full support).
Prints positivity a = adj(C)1, D = det K, the Hurwitz minors of det(mu I + K), and the limit test on omega."""
import sys, sympy as sp
from certify_M_violation import from_d6
s, beta = sys.argv[1], sp.Rational(sys.argv[2])
ts = [sp.Rational(x) for x in sys.argv[3:]]
A = from_d6(s); m = len(A)
N = sp.Matrix(m, m, lambda i, j: 1 if (i != j and A[j][i] == 0) else 0)
I = sp.eye(m); J = sp.ones(m, m); one = sp.ones(m, 1)
C = I + beta * N
a = C.adjugate() * one
print('a = adj(C)1 all > 0:', all(x > 0 for x in a), ' det C =', C.det())
mu = sp.symbols('mu')
def hurwitz(coeffs):  # high->low
    n = len(coeffs) - 1
    c = lambda k: coeffs[k] if 0 <= k <= n else 0
    return [sp.Matrix(k, k, lambda i, j: c(2 * (i + 1) - (j + 1))).det() for k in range(1, n + 1)]
for t in ts:
    K = C + t * J
    p = sp.Poly((mu * I + K).det(), mu).all_coeffs()
    H = hurwitz(p)
    print(f't={t}: D=det K={K.det()}  Hurwitz minor signs {[sp.sign(h) for h in H]}  -> stable={all(h > 0 for h in H)}')
# limit: omega = det(mu I + C + J) - det(mu I + C)
om = sp.Poly(sp.expand((mu * I + C + J).det() - (mu * I + C).det()), mu).all_coeffs()
H = hurwitz(om)
print('limit omega Hurwitz minor signs', [sp.sign(h) for h in H])
