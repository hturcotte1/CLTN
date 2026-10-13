"""Sanity (not proof): Lemma I.7.
 (i) chi_K(lam) == (-beta)^m [chi_A(nu) + (1 + t/beta) m chi_{A|perp}(nu)],  nu = (1-beta-lam)/beta   (symbolic in beta,t)
 (ii) h := gcd(chi_A, chi_{A+J}); every root mu of h gives eigenvalue 1-beta-beta*mu of K at all (beta,t)
 (iii) if h has a root with Re mu >= -1/2: ctln.analyze_support never reports stable at random legal points."""
import sys, random
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import sympy as sp
import ctln
random.seed(99)
lam, nu, be, tt = sp.symbols('lam nu beta t')
def chi(X, v): return sp.expand((v*sp.eye(X.shape[0]) - X).det(method='berkowitz'))
st = dict(ident=0, h_nontriv=0, killed=0, killed_checks=0, twin_cases=0)
for trial in range(300):
    m = random.randint(2, 6)
    p = random.choice([0.5, 0.7, 0.85])
    A = sp.Matrix(m, m, lambda i, j: 0 if i == j else int(random.random() < p))
    if random.random() < 0.4 and m >= 3:      # plant non-adjacent or bidirected in-twins 0,1
        c = random.randint(0, 1); A[0, 1] = A[1, 0] = c
        for k in range(2, m): A[k, 1] = A[k, 0]
        st['twin_cases'] += 1
    I, J = sp.eye(m), sp.ones(m, m)
    K = I + be*(J - I - A.T) + tt*J
    chiA = chi(A, nu); chiAJ = chi(A + J, nu); chip = sp.expand((chiA - chiAJ)/m)
    rhs = sp.expand((-be)**m * (chiA + (1 + tt/be)*m*chip).subs(nu, (1 - be - lam)/be))
    if m <= 4:
        assert sp.expand(chi(K, lam) - rhs) == 0; st['ident'] += 1
    h = sp.gcd(sp.Poly(chiA, nu), sp.Poly(chiAJ, nu))
    if h.degree() == 0: continue
    st['h_nontriv'] += 1
    roots = [r for f, _ in sp.factor_list(h.as_expr())[1] for r in sp.Poly(f, nu).all_roots()]
    # (ii) at random rational (beta, t): 1-beta-beta*mu annihilates chi_K
    b0 = sp.Rational(random.randint(101, 500), 100); t0 = sp.Rational(random.randint(1, 1000), 7)
    cK = chi(K.subs({be: b0, tt: t0}), lam)
    for r in roots:
        val = sp.N(cK.subs(lam, 1 - b0 - b0*r), 40)
        assert abs(val) < 1e-20, val
    if max(sp.re(sp.N(r, 30)) for r in roots) >= -sp.Rational(1, 2) - sp.Rational(1, 10**20):
        st['killed'] += 1
        Al = [[int(A[i, j]) for j in range(m)] for i in range(m)]
        for _ in range(4):
            eps = F(1, random.randint(2, 2000)); q = F(random.randint(101, 600), 100)
            delta = q*eps/(1 - eps)
            assert ctln.legal(eps, delta)
            rec = ctln.analyze_support(Al, list(range(m)), eps, delta)
            assert rec['stable'] is not True, (Al, eps, delta)
            st['killed_checks'] += 1
print("I.7 sanity checks passed:", st)
