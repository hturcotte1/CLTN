"""Sanity (not proof): Lemma I.3 (a),(b),(c),(g) with exact rationals.
 (a) chi_{C|perp}(lam) == (-beta)^(m-1) chi_{A|perp}((1-beta-lam)/beta), using m chi_{X|perp} = chi_X - chi_{X+J} (Cor I.1.1)
 (g) K == (1-beta) I - beta (A^T - tau J), tau = 1 + t/beta
 (b,c) limit verdict [chi_{A|perp}(z - gamma) Hurwitz] vs exact Hurwitz verdict of det(mu I + K) at t = 10^9."""
import sys, random
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import sympy as sp
import ctln
random.seed(3)
lam, z = sp.symbols('lam z')
def chi(X):
    return sp.expand((lam*sp.eye(X.shape[0]) - X).det(method='berkowitz'))
def chi_perp(X):
    m = X.shape[0]
    return sp.expand((chi(X) - chi(X + sp.ones(m, m)))/m)
def Q(x): return sp.Rational(x.numerator, x.denominator)
def toF(c): c = sp.Rational(c); return F(int(c.p), int(c.q))
stats = dict(a=0, g=0, agree=0, disagree=0, boundary=0, limit_stable=0)
for trial in range(250):
    m = random.randint(2, 6)
    p = random.choice([0.5, 0.75, 0.9])
    A = sp.Matrix(m, m, lambda i, j: 0 if i == j else int(random.random() < p))
    beta = F(random.randint(101, 700), 100); gam = 1 - 1/beta; b = Q(beta)
    I, J = sp.eye(m), sp.ones(m, m)
    N = J - I - A.T; C = I + b*N
    # (a)
    lhs = chi_perp(C)
    rhs = sp.expand((-b)**(m-1) * chi_perp(A).subs(lam, (1 - b - lam)/b))
    assert sp.expand(lhs - rhs) == 0; stats['a'] += 1
    # (g)
    t = F(10)**9; tau = 1 + t/beta
    K = C + Q(t)*J
    assert K == (1 - b)*I - b*(A.T - Q(tau)*J); stats['g'] += 1
    # limit verdict
    h = sp.Poly(chi_perp(A).subs(lam, z - Q(gam)), z)
    lim, _ = ctln.hurwitz_stable([toF(c) for c in h.all_coeffs()])
    mu = sp.symbols('mu')
    pt = sp.Poly(sp.expand((mu*I + K).det(method='berkowitz')), mu)
    fin, _ = ctln.hurwitz_stable([toF(c) for c in pt.all_coeffs()])
    if lim is None: stats['boundary'] += 1; continue
    stats['limit_stable'] += bool(lim)
    if lim == fin: stats['agree'] += 1
    else: stats['disagree'] += 1; print("disagree (could be T(beta) > 1e9):", A.tolist(), beta, lim, fin)
print(stats)
