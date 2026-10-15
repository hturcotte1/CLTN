"""Independent exact test of Lemma D1 / D2 (notes_theory.md section 1). No floats in any decision."""
import sys, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
import ctln, fastregion as fr

z = sp.symbols('z')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def rand_graph(n, p):
    return [[0 if i == j else (1 if random.random() < p else 0) for j in range(n)] for i in range(n)]

def psi(X):  # 1^T adj(zI - X) 1, via chi_X - chi_{X+J}
    m = X.shape[0]
    J = sp.ones(m, m)
    return sp.expand((z*sp.eye(m) - X).det(method='berkowitz') - (z*sp.eye(m) - X - J).det(method='berkowitz'))

def compress_charpoly(X):  # det(zI - P X | 1^perp) using basis V = e_i - e_m
    m = X.shape[0]
    V = sp.zeros(m, m-1)
    for i in range(m-1):
        V[i, i] = 1; V[m-1, i] = -1
    Xc = (V.T*V).inv() * V.T * X * V
    return sp.expand((z*sp.eye(m-1) - Xc).det(method='berkowitz'))

def shifted_hurwitz(p, r):
    """all roots of p have Re > -r  <=> q(y) := p(-y - r) has all roots Re<0. Exact three-valued."""
    q = sp.Poly(sp.expand(p.subs(z, -sp.Symbol('y') - r)), sp.Symbol('y'))
    c = [F(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in q.all_coeffs()]
    if c[0] < 0: c = [-x for x in c]
    return ctln.hurwitz_stable(c)[0]

stats = {"D1_poly": 0, "D1_compress": 0, "D1_stab_agree": 0, "D1_stab_none": 0, "D2": 0,
         "D2_pos": 0, "D2_singular": 0}
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 300):
    n = random.randint(2, 7)
    A = rand_graph(n, random.choice([0.5, 0.7, 0.85, 0.95]))
    As = sp.Matrix(A); m = n
    J = sp.ones(m, m); I = sp.eye(m)
    D = J - I - As
    # D1 identity: A|perp = -I - D|perp as charpolys: chi_{A|perp}(z) = (-1)^(m-1) chi_{D|perp}(-1-z)
    pA = psi(As); pD = psi(D)
    assert sp.expand(pA - (-1)**(m-1) * pD.subs(z, -1 - z)) == 0, (A,)
    stats["D1_poly"] += 1
    if m >= 2:
        assert sp.expand(pD - m*compress_charpoly(D)) == 0
        assert sp.expand(pA - m*compress_charpoly(As)) == 0
        stats["D1_compress"] += 1
        # limit stability: repo omega (Hurwitz in mu) vs psi_D roots Re > -1/beta
        bd = fr.beta_data(A, list(range(n)))
        for _ in range(3):
            beta = F(random.randint(101, 600), random.randint(50, 200))
            if beta <= 1: continue
            coeffs = [fr.peval(p, beta) for p in fr.stability_coeffs(bd, fr.INF)]
            v_repo = ctln.hurwitz_stable(coeffs)[0]
            v_D = shifted_hurwitz(pD, sp.Rational(beta.denominator, beta.numerator))
            # also direct: alpha(A|perp) < -(1-1/beta) via psi_A roots Re < -gamma
            gam = 1 - 1/beta
            qA = sp.Poly(sp.expand(pA.subs(z, sp.Symbol('y') - sp.Rational(gam.numerator, gam.denominator))), sp.Symbol('y'))
            cA = [F(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in qA.all_coeffs()]
            v_A = ctln.hurwitz_stable(cA)[0]
            assert v_repo == v_D == v_A, (A, beta, v_repo, v_D, v_A)
            stats["D1_stab_agree"] += 1
            if v_repo is None: stats["D1_stab_none"] += 1
    # D2 at rational beta
    for _ in range(3):
        beta = F(random.randint(101, 900), random.randint(50, 300))
        if beta <= 1: continue
        r = 1/beta; gam = 1 - r
        Mg = [[F(A[j][i]) + (gam if i == j else 0) for j in range(m)] for i in range(m)]   # A^T + gamma I
        Mu = [[F(int(D[j, i])) + (r if i == j else 0) for j in range(m)] for i in range(m)] # D^T + r I
        MuS = sp.Matrix(m, m, lambda i, j: sp.Rational(Mu[i][j].numerator, Mu[i][j].denominator))
        assert MuS.det() != 0
        u = MuS.LUsolve(sp.ones(m, 1)); U = sum(u)
        MgS = sp.Matrix(m, m, lambda i, j: sp.Rational(Mg[i][j].numerator, Mg[i][j].denominator))
        if MgS.det() == 0:
            assert U == 1; stats["D2_singular"] += 1
            continue
        assert U != 1
        w = MgS.LUsolve(sp.ones(m, 1))
        assert all(sp.simplify(w[i] - u[i]/(U-1)) == 0 for i in range(m))
        lhs = all(x > 0 for x in w)
        rhs = all(x > 0 for x in u) and U > 1
        assert lhs == rhs, (A, beta, list(w), list(u), U)
        # cross-check against repo positivity (adj(C)1 one-signed) at t = INF
        bd = fr.beta_data(A, list(range(n)))
        s0 = fr.peval(bd["s"], beta)
        if s0 != 0:
            rp = all(fr.peval(ai, beta)/s0 > 0 for ai in bd["a"])
            assert rp == lhs, (A, beta)
        stats["D2"] += 1
        if lhs: stats["D2_pos"] += 1
print(stats)
