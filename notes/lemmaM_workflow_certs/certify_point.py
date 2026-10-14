"""Exact certification of a (graph, beta, t1, t2) Lemma-M flag with ctln.analyze_support.
eps = 1/(t+1), delta = (beta-1)*eps (Fractions).  Also exact alpha-side check via sympy:
limit-stability at beta (Lemma I) decided exactly by Hurwitz on omega(mu; beta) via fastregion."""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln

def point(A, beta, t):
    beta = F(beta); t = F(t)
    eps = 1 / (t + 1); delta = (beta - 1) * eps
    r = ctln.analyze_support(A, list(range(len(A))), eps, delta)
    q = ctln.q_of(eps, delta)
    return dict(beta=str(beta), t=str(t), eps=str(eps), delta=str(delta), q=str(q), legal=bool(q > 1),
                stable=r['stable'], detail=r['stability_detail'], positive=r['positive'],
                x=r['x'], stable_fixed_point=r['stable_fixed_point'], charpoly=r['charpoly'])

def limit_stable(A, beta):
    """Exact: Hurwitz test on omega(mu) = 1^T adj(mu I + C) 1 = (det(mu I + C + J) - det(mu I + C)), C = I + beta N."""
    import sympy as sp
    n = len(A); beta = sp.Rational(str(F(beta)))
    mu = sp.symbols('mu')
    N = sp.Matrix(n, n, lambda i, j: 0 if i == j else (1 - A[j][i]))
    C = sp.eye(n) + beta * N
    J = sp.ones(n, n)
    om = sp.expand((mu * sp.eye(n) + C + J).det() - (mu * sp.eye(n) + C).det())
    P = sp.Poly(om, mu)
    coeffs = [F(int(sp.fraction(c)[0]), int(sp.fraction(c)[1])) for c in P.all_coeffs()]
    return ctln.hurwitz_stable(coeffs), [str(c) for c in coeffs]

if __name__ == '__main__':
    A = json.loads(sys.argv[1]); beta = sys.argv[2]
    for t in sys.argv[3:]:
        print(json.dumps(point(A, beta, t)))
    print('limit (t=inf):', limit_stable(A, beta))
