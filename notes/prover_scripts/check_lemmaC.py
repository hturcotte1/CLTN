"""Sanity checks (NOT proofs) for Lemma C (t -> infinity) and Lemma F (ray monotonicity).
Decisions use exact Hurwitz tests (ctln.hurwitz_stable); nroots only for the printed asymptotics."""
import random, sys
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import ctln

mu = sp.symbols('mu')
random.seed(11)

def charpoly_K(N, beta, t):
    m = len(N)
    K = sp.Matrix(m, m, lambda i, j: (1 if i == j else 0) + beta * N[i][j] + t)
    return [sp.Rational(c) for c in sp.Poly((mu * sp.eye(m) + K).det(), mu).all_coeffs()]

def omega_coeffs(N, beta):
    m = len(N)
    C = sp.Matrix(m, m, lambda i, j: (1 if i == j else 0) + beta * N[i][j])
    om = (sp.ones(1, m) * (mu * sp.eye(m) + C).adjugate() * sp.ones(m, 1))[0]
    return [sp.Rational(c) for c in sp.Poly(sp.expand(om), mu).all_coeffs()]

agree = disagree_small_t = undecided = 0
for trial in range(150):
    m = random.randint(2, 6)
    N = [[0 if i == j else int(random.random() < 0.3) for j in range(m)] for i in range(m)]
    beta = sp.Rational(random.randint(201, 400), 100)
    om_st, _ = ctln.hurwitz_stable([F(int(c.p), int(c.q)) for c in omega_coeffs(N, beta)])
    if om_st is None:
        undecided += 1; continue
    T = sp.Integer(10) ** 7
    p_st, _ = ctln.hurwitz_stable([F(int(c.p), int(c.q)) for c in charpoly_K(N, beta, T)])
    if p_st == om_st: agree += 1
    else: disagree_small_t += 1; print("needs larger t:", N, beta, om_st, p_st)
print(f"Lemma C: omega-Hurwitz verdict equals p_t-Hurwitz verdict at t=1e7 in {agree} cases; {disagree_small_t} disagreements; {undecided} omega boundary cases skipped")

# eigenvalue asymptotics lambda_max(K) - m t -> 1 + beta e/m
N = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]
beta = sp.Rational(54, 25); m = 6; e = sum(map(sum, N))
for t in (10**3, 10**5):
    r = sp.Poly(charpoly_K(N, beta, t), mu).nroots(n=30)
    lam = [-x for x in r]
    big = max(lam, key=lambda z: sp.re(z))
    print(f"t={t}: largest eigenvalue - m t = {sp.N(sp.re(big) - m*t, 12)}, predicted limit 1+beta e/m = {sp.N(1 + beta*e/m, 12)}")

# Lemma F: stability of K(beta,t) is preserved under (beta,t) -> (lam beta, lam t), 0<lam<1
cnt = 0
for trial in range(200):
    m = random.randint(2, 6)
    N = [[0 if i == j else int(random.random() < 0.4) for j in range(m)] for i in range(m)]
    beta = sp.Rational(random.randint(100, 600), 100); t = sp.Rational(random.randint(1, 300), 10)
    st, _ = ctln.hurwitz_stable([F(int(c.p), int(c.q)) for c in charpoly_K(N, beta, t)])
    if st is True:
        lam = sp.Rational(random.randint(1, 99), 100)
        st2, _ = ctln.hurwitz_stable([F(int(c.p), int(c.q)) for c in charpoly_K(N, lam * beta, lam * t)])
        assert st2 is True, (N, beta, t, lam)
        cnt += 1
print(f"Lemma F: {cnt} stable (beta,t) pairs, all remain stable after scaling toward the origin")
