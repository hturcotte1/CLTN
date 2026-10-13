"""Sanity checks (NOT proofs) for Lemma A and the compression identity (Lemma C').
Exact sympy arithmetic only. Random graphs, all support sizes 1..6, n up to 7."""
import random, sys
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import ctln

beta, t, mu, eps, delta = sp.symbols('beta t mu eps delta')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def rand_graph(n, p=0.5):
    return [[0 if i == j else int(random.random() < p) for j in range(n)] for i in range(n)]

def check(A, sigma):
    m, n = len(sigma), len(A)
    N = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - A[sigma[j]][sigma[i]])
    I, J, one = sp.eye(m), sp.ones(m, m), sp.ones(m, 1)
    C = I + beta * N
    K = C + t * J
    # (i) (I - W)_sigma == eps K  with t=(1-eps)/eps, beta = 1 + delta/eps
    W = sp.Matrix(n, n, lambda i, j: 0 if i == j else (-1 + eps if A[j][i] == 1 else -1 - delta))
    IW = sp.eye(m) - W.extract(sigma, sigma)
    lhs = (eps * K).subs({t: (1 - eps) / eps, beta: 1 + delta / eps})
    assert sp.simplify(IW - lhs) == sp.zeros(m, m), "(i)"
    # (ii) adj(K) 1 == adj(C) 1
    a = C.adjugate() * one
    assert sp.expand(K.adjugate() * one - a) == sp.zeros(m, 1), "(ii)"
    for ai in a:
        assert sp.Poly(ai, beta).degree() <= m - 1 if ai != 0 else True
    # (iii)
    d, s = C.det(), sum(a)
    assert sp.expand(K.det() - (d + t * s)) == 0, "(iii)"
    # (iv)
    phi = (mu * I + C).det()
    omega = sp.expand((one.T * (mu * I + C).adjugate() * one)[0])
    assert sp.expand((mu * I + K).det() - (phi + t * omega)) == 0, "(iv)"
    if m >= 1:
        P = sp.Poly(omega, mu)
        assert P.degree() == m - 1 and P.LC() == m, "(iv) degree/leading coeff"
    # (v) off-support drive, at random rational parameters with det K != 0
    th = sp.Rational(random.randint(1, 9), random.randint(1, 9))
    e0 = sp.Rational(random.randint(1, 30), 31); dl0 = sp.Rational(random.randint(1, 40), 7)
    t0, b0 = (1 - e0) / e0, 1 + dl0 / e0
    Kn = K.subs({t: t0, beta: b0}); Dn = Kn.det()
    if Dn != 0:
        x = Kn.LUsolve((th / e0) * one)
        Wn = W.subs({eps: e0, delta: dl0})
        for k in range(n):
            if k in sigma:
                continue
            drive = th + sum(Wn[k, sigma[j]] * x[j] for j in range(m))
            Nk = [1 - A[sigma[j]][k] for j in range(m)]
            ok = (d - beta * sum(Nk[j] * a[j] for j in range(m))).subs(beta, b0)
            assert sp.simplify(drive - th * ok / Dn) == 0, "(v)"
    # Lemma C': omega(mu) = m det((mu+1) I_{m-1} + beta G), G = matrix of P N on 1-perp
    if m >= 2:
        V = sp.Matrix(m, m - 1, lambda i, j: 1 if i == j else (-1 if i == m - 1 else 0))
        G = (V.T * V).inv() * V.T * N * V
        rhs = m * ((mu + 1) * sp.eye(m - 1) + beta * G).det()
        assert sp.expand(omega - rhs) == 0, "C' compression identity"
        # s = omega(0) = m det(I + beta G)
        assert sp.expand(s - m * (sp.eye(m - 1) + beta * G).det()) == 0, "s identity"
        # trace of G = -(#ones in N)/m
        assert G.trace() == -sp.Rational(sum(N), m), "trace"
    return True

cnt = 0
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 40):
    n = random.randint(2, 7)
    A = rand_graph(n, random.choice([0.3, 0.5, 0.7]))
    m = random.randint(1, min(n, 6))
    sigma = sorted(random.sample(range(n), m))
    check(A, sigma); cnt += 1
print("Lemma A (i)-(v) and C' identities: all", cnt, "random cases pass")
