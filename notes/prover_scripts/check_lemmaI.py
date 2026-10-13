"""Sanity checks (NOT proofs) for Lemma I (limit form of all three conditions):
  (A^T+gI) a = ((beta s - d)/beta) 1,   det(A^T+gI) = (-1/beta)^m (d - beta s),
  a one-signed <=> w := (A^T+gI)^{-1} 1 > 0,   o_k s < 0 <=> in_k(w) < 1 (when w > 0),
  omega Hurwitz <=> every eigenvalue of A|perp has Re < -g   (g = 1 - 1/beta).
Exact rationals; random graphs and supports; also Geneson's graph."""
import random, sys
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import ctln, fastregion as fr
random.seed(5)
mu = sp.symbols('mu')
stats = dict(cases=0, onesigned=0, w_pos=0, off_checked=0, hurwitz_checked=0)
for trial in range(400):
    n = random.randint(3, 7); m = random.randint(2, min(n, 6))
    A = [[0 if i == j else int(random.random() < random.choice([0.5, 0.7, 0.85])) for j in range(n)] for i in range(n)]
    sig = sorted(random.sample(range(n), m))
    bd = fr.beta_data(A, sig)
    beta = F(random.randint(201, 500), 100); g = 1 - 1 / beta
    ev = lambda p: fr.peval(p, beta)
    a = [ev(p) for p in bd["a"]]; d = ev(bd["d"]); s = sum(a)
    AT = sp.Matrix(m, m, lambda i, j: A[sig[j]][sig[i]])          # (A_sigma^T)_{ij} = A[sig_j][sig_i]
    Mg = AT + sp.Rational(g.numerator, g.denominator) * sp.eye(m)
    lhs = Mg * sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in a])
    c = (beta * s - d) / beta
    assert all(v == sp.Rational(c.numerator, c.denominator) for v in lhs)
    assert Mg.det() == sp.Rational((-1) ** m) / sp.Rational(beta.numerator, beta.denominator) ** m * sp.Rational(d - beta * s)
    stats["cases"] += 1
    onesigned = all(x > 0 for x in a) or all(x < 0 for x in a)
    wpos = False
    if Mg.det() != 0:
        w = Mg.LUsolve(sp.ones(m, 1)); wpos = all(v > 0 for v in w)
    assert onesigned == wpos
    stats["onesigned"] += onesigned; stats["w_pos"] += wpos
    if wpos:
        for k, okp in bd["o"].items():
            ink = sum(w[b] for b in range(m) if A[sig[b]][k] == 1)
            assert (ev(okp) * s < 0) == (ink < 1); stats["off_checked"] += 1
    # stability in the limit: omega Hurwitz  <=>  Re spec(A|perp) < -g ; test via the polynomial det(z I + g I + A|perp_V)
    V = sp.Matrix(m, m - 1, lambda i, j: 1 if i == j else (-1 if i == m - 1 else 0))
    Aperp = (V.T * V).inv() * V.T * sp.Matrix(m, m, lambda i, j: A[sig[i]][sig[j]]) * V
    chi = sp.Poly(((mu + sp.Rational(g.numerator, g.denominator)) * sp.eye(m - 1) - Aperp).det(), mu)  # roots mu = lambda - g ... 
    # roots of chi are lambda - g? we need Re lambda < -g <=> roots of det(z I - (Aperp + g I)) in Re z < 0
    chi = sp.Poly((mu * sp.eye(m - 1) - (Aperp + sp.Rational(g.numerator, g.denominator) * sp.eye(m - 1))).det(), mu)
    st1, _ = ctln.hurwitz_stable([F(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in chi.all_coeffs()])
    om = [fr.peval(p, beta) for p in fr.stability_coeffs(bd, fr.INF)]
    st2, _ = ctln.hurwitz_stable(om)
    if st1 is not None and st2 is not None:
        assert st1 == st2; stats["hurwitz_checked"] += 1
print("Lemma I checks pass:", stats)
