"""Sanity (not proof): Lemma I.4 with exact rationals, against ctln.analyze_support (independent code path).
Checks (I.4.0)(ii), (b1), (c) P1<=>P2<=>P3 + all formulas, (e) sign rule, (f) off-support formula."""
import sys, random
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import sympy as sp
import ctln
random.seed(41)
def Q(x): return sp.Rational(x.numerator, x.denominator)
st = dict(cases=0, P_true=0, w_pos_a_neg=0, offsupport=0, d_zero=0)
for trial in range(600):
    n = random.randint(2, 7); m = random.randint(1, n)
    p = random.choice([0.5, 0.7, 0.85, 0.95])
    G = [[0 if i == j else int(random.random() < p) for j in range(n)] for i in range(n)]
    sig = sorted(random.sample(range(n), m))
    A = sp.Matrix(m, m, lambda a, b: G[sig[a]][sig[b]])
    beta = F(random.randint(101, 600), 100); t = F(random.randint(1, 4000), random.choice([1, 7, 100]))
    eps = 1/(1 + t); delta = (beta - 1)/(1 + t); gam = 1 - 1/beta
    b, g, tt = Q(beta), Q(gam), Q(t)
    I, J = sp.eye(m), sp.ones(m, m)
    N = J - I - A.T; C = I + b*N; K = C + tt*J; Mg = A.T + g*I
    a = C.adjugate()*sp.ones(m, 1); d = C.det()
    # (I.4.0)(ii)
    assert a == (-b)**(m-1) * Mg.adjugate()*sp.ones(m, 1)
    detM = Mg.det()
    assert detM != 0            # rational beta: -gamma in (-1,0) is never an eigenvalue of the integer matrix A^T
    w = Mg.LUsolve(sp.ones(m, 1)); W1 = sum(w)
    if d == 0: st['d_zero'] += 1
    # (b1)
    assert d == (-b)**m * detM * (1 - W1)
    if W1 != 1:
        u = w / (b*(W1 - 1)); assert C*u == sp.ones(m, 1)
    P2 = all(v > 0 for v in a) or all(v < 0 for v in a)
    P3 = all(v > 0 for v in w)
    P3p = P3 or all(v < 0 for v in w)
    rec = ctln.analyze_support(G, sig, eps, delta)
    P1 = rec['x'] is not None and rec['positive']
    assert P1 == P2 == P3 == P3p, (P1, P2, P3, P3p)
    st['cases'] += 1
    if P3:
        st['P_true'] += 1
        Den = b*(W1 - 1) + tt*W1
        assert W1 > 1 and C*w == b*(W1 - 1)*sp.ones(m, 1) and K*w == Den*sp.ones(m, 1)
        x = [Q(F(v)) for v in rec['x']]
        assert all(x[i] == (1 + tt)*w[i]/Den for i in range(m))          # theta = 1
        assert a == (-b)**(m-1)*detM*w
        assert d == (-b)**(m-1)*detM*b*(W1 - 1)
        assert K.det() == (-b)**(m-1)*detM*Den
        if all(v < 0 for v in a):
            st['w_pos_a_neg'] += 1
            assert K.det() < 0 and rec['stable'] is False
        # (f) off-support
        for k in range(n):
            if k in sig: continue
            ink = sum(w[j] for j in range(m) if G[sig[j]][k] == 1)
            Wfull = ctln.weight_matrix(G, eps, delta)
            y = 1 + sum(Q(Wfull[k][sig[j]])*x[j] for j in range(m))
            assert y == b*(ink - 1)/Den; st['offsupport'] += 1
print("I.4 sanity checks passed:", st)
