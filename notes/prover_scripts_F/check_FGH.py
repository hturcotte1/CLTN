"""
Sanity checks (NOT proofs) for notes/lemma_F_reduction.md, Lemmas F, G, H, R2, P, E.
Exact arithmetic only (Fraction / sympy). Run:  python3 -I check_FGH.py
Every assertion failure aborts with a message.
"""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
from itertools import combinations, permutations
import sympy as sp
import ctln

rng = random.Random(20261013)

def rand_graph(n, p=None):
    p = rng.random() if p is None else p
    return [[0 if i == j else (1 if rng.random() < p else 0) for j in range(n)] for i in range(n)]

def rand_legal():
    # random rational eps in (0,1), then delta with q > 1
    while True:
        eps = F(rng.randint(1, 99), 100) * F(1, rng.choice([1, 1, 2, 10, 100]))
        delta = F(rng.randint(1, 400), rng.randint(1, 100))
        if ctln.legal(eps, delta):
            return eps, delta

def I_minus(M):
    m = len(M)
    return [[(F(1) if i == j else F(0)) - M[i][j] for j in range(m)] for i in range(m)]

def induced(A, sigma):
    return [[A[i][j] for j in sigma] for i in sigma]

# ---------------------------------------------------------------- H0: closed form I - W = eps K
def check_H0(trials=300):
    for _ in range(trials):
        n = rng.randint(1, 7); A = rand_graph(n); eps, delta = rand_legal()
        W = ctln.weight_matrix(A, eps, delta)
        t = (1 - eps) / eps; beta = 1 + delta / eps
        N = [[0 if i == j else 1 - A[j][i] for j in range(n)] for i in range(n)]   # receiver-first nonedge
        K = [[(1 if i == j else 0) + beta * N[i][j] + t for j in range(n)] for i in range(n)]
        assert I_minus(W) == [[eps * K[i][j] for j in range(n)] for i in range(n)], "H0 failed"
    print("H0  I - W = eps (I + beta N + t J), N = J - I - A^T : OK", trials)

# ---------------------------------------------------------------- F: restriction lemma
def check_F(trials=400):
    for _ in range(trials):
        n = rng.randint(2, 7); A = rand_graph(n); eps, delta = rand_legal()
        k = rng.randint(1, n); sigma = sorted(rng.sample(range(n), k))
        W = ctln.weight_matrix(A, eps, delta)
        H = induced(A, sigma)
        WH = ctln.weight_matrix(H, eps, delta)
        assert ctln.submatrix(W, sigma) == WH, "W_sigma != W(G|sigma)"
        rG = ctln.analyze_support(A, sigma, eps, delta)
        rH = ctln.analyze_support(H, list(range(k)), eps, delta)
        for key in ("x", "positive", "stable", "charpoly", "clique", "BOUNDARY"):
            assert rG.get(key) == rH.get(key), ("F mismatch", key)
        # stable fixed point in G  ==>  full support of H satisfies (1)&(3)
        if rG["stable_fixed_point"]:
            assert rH["positive"] and rH["stable"] is True
    print("F   W_sigma = W(G|sigma); (1),(3) verdicts identical : OK", trials)

def check_perm(trials=200):
    for _ in range(trials):
        n = rng.randint(2, 6); A = rand_graph(n); eps, delta = rand_legal()
        p = list(range(n)); rng.shuffle(p)
        Ap = [[A[p[i]][p[j]] for j in range(n)] for i in range(n)]   # new i = old p[i]
        r = ctln.analyze_support(A, list(range(n)), eps, delta)
        rp = ctln.analyze_support(Ap, list(range(n)), eps, delta)
        assert r["positive"] == rp["positive"] and r["stable"] == rp["stable"] and r["charpoly"] == rp["charpoly"]
        if r["x"] is not None:
            assert rp["x"] == [r["x"][p[i]] for i in range(n)]
    print("F'  relabelling invariance of (1),(3) : OK", trials)

# ---------------------------------------------------------------- R2: clique rule
def check_R2_symbolic(mmax=7):
    e, d, th = sp.symbols("epsilon delta theta", positive=True)
    for m in range(1, mmax + 1):
        IW = sp.Matrix(m, m, lambda i, j: 1 if i == j else 1 - e)      # all W_ij = -1+eps
        assert sp.simplify(IW - (e * sp.eye(m) + (1 - e) * sp.ones(m, m))) == sp.zeros(m, m)
        c = th / (e + m * (1 - e))
        assert sp.simplify(IW * sp.ones(m, 1) * c - th * sp.ones(m, 1)) == sp.zeros(m, 1)
        lam = sp.symbols("lam")
        cp = sp.factor((lam * sp.eye(m) - (-IW)).det())   # char poly of -I + W_sigma = -(I - W_sigma)
        target = (lam + e) ** (m - 1) * (lam + e + m * (1 - e))
        assert sp.expand(cp - target) == 0, ("R2 charpoly", m, cp)
        for r in range(0, m + 1):
            drive = th + c * ((m - r) * (-1 + e) + r * (-1 - d))
            assert sp.simplify(drive - th * (e - r * (e + d)) / (e + m * (1 - e))) == 0
    print("R2  clique: x = theta/(eps+m(1-eps)) 1, eigenvalues -eps^(m-1), -(eps+m(1-eps)), drive formula : OK m<=", mmax)

def check_R2_numeric(trials=300):
    for _ in range(trials):
        n = rng.randint(2, 7); A = rand_graph(n); eps, delta = rand_legal()
        m = rng.randint(1, n); sigma = sorted(rng.sample(range(n), m))
        for i in sigma:
            for j in sigma:
                if i != j: A[i][j] = 1
        r = ctln.analyze_support(A, sigma, eps, delta)
        assert r["clique"] and r["positive"] and r["stable"] is True
        assert r["stable_fixed_point"] == (not r["targets"]), "Rule 2 failed"
    print("R2  clique sigma: stable FP support iff no target (ctln) : OK", trials)

# ---------------------------------------------------------------- G: theta scaling
def check_G(trials=300):
    keys = ("positive", "off_support_ok", "stable", "stability_detail", "charpoly",
            "is_fixed_point", "stable_fixed_point", "agrees", "BOUNDARY")
    for _ in range(trials):
        n = rng.randint(1, 6); A = rand_graph(n); eps, delta = rand_legal()
        m = rng.randint(1, n); sigma = sorted(rng.sample(range(n), m))
        th = F(rng.randint(1, 1000), rng.randint(1, 1000))
        r1 = ctln.analyze_support(A, sigma, eps, delta, F(1))
        r2 = ctln.analyze_support(A, sigma, eps, delta, th)
        for k in keys:
            assert r1.get(k) == r2.get(k), ("theta changed", k)
        if r1["x"] is not None:
            assert [F(v) * th for v in r1["x"]] == [F(v) for v in r2["x"]]
    print("G   theta only rescales x; every verdict field identical : OK", trials)

# ---------------------------------------------------------------- E: embedding / inside-out domination
def check_E(trials=3000):
    hits = 0
    for _ in range(trials):
        k = rng.randint(1, 6); H = rand_graph(k); eps, delta = rand_legal()
        IW = I_minus(ctln.weight_matrix(H, eps, delta))
        x = ctln.solve_exact(IW, [F(1)] * k)
        if x is None or not all(v > 0 for v in x):
            continue
        hits += 1
        # every in-neighbourhood S of an extra node k*: test domination criterion and monotonicity
        drive = {}
        for r in range(k + 1):
            for S in combinations(range(k), r):
                S = set(S)
                drive[frozenset(S)] = 1 - sum(((1 - eps) if j in S else (1 + delta)) * x[j] for j in range(k))
        for S, dv in drive.items():
            dom = any((i not in S) and all(H[j][i] == 1 for j in S) for i in range(k))
            if dom:
                assert dv < 0, "domination criterion failed"
            for j in S:
                assert drive[S - {j}] <= dv, "monotonicity failed"
        assert drive[frozenset()] < 0
        # cross-check with ctln on the padded graph G = H + node with in-set S
        S = set(j for j in range(k) if rng.random() < 0.5)
        G = [row[:] + [1 if j in S else 0] for j, row in enumerate(H)] + [[rng.randint(0, 1) for _ in range(k)] + [0]]
        r = ctln.analyze_support(G, list(range(k)), eps, delta)
        assert r["off_support_ok"] == (drive[frozenset(S)] < 0)
    print("E   inside-out domination => drive<0; drive monotone in in-set; empty in-set always OK : OK", hits, "positive cases")

# ---------------------------------------------------------------- H: rank-one / Sherman-Morrison
def check_H_symbolic(mmax=4):
    t = sp.symbols("t")
    for m in range(1, mmax + 1):
        C = sp.Matrix(m, m, lambda i, j: sp.Symbol(f"c{i}{j}"))
        one = sp.ones(m, 1)
        K = C + t * one * one.T
        u = C.adjugate() * one
        assert sp.expand(K.adjugate() * one - u) == sp.zeros(m, 1), "adj(C+tJ)1 != adj(C)1"
        assert sp.expand(K.det() - (C.det() + t * (one.T * u)[0])) == 0, "det lemma"
    print("H   adj(C+tJJ)1 = adj(C)1 and det(C+tJ) = det C + t 1'adj(C)1 for GENERIC C : OK m<=", mmax)

def check_H_numeric(trials=500):
    for _ in range(trials):
        n = rng.randint(1, 6); A = rand_graph(n); eps, delta = rand_legal()
        t = (1 - eps) / eps; beta = 1 + delta / eps
        N = [[0 if i == j else 1 - A[j][i] for j in range(n)] for i in range(n)]
        C = sp.Matrix(n, n, lambda i, j: (1 if i == j else 0) + sp.Rational(beta.numerator, beta.denominator) * N[i][j])
        tt = sp.Rational(t.numerator, t.denominator)
        K = C + tt * sp.ones(n, n)
        u = C.adjugate() * sp.ones(n, 1); dK = K.det()
        rhs = all(ui != 0 for ui in u) and dK != 0 and all(sp.sign(ui) == sp.sign(dK) for ui in u)
        IW = I_minus(ctln.weight_matrix(A, eps, delta))
        x = ctln.solve_exact(IW, [F(1)] * n)
        lhs = x is not None and all(v > 0 for v in x)
        assert lhs == rhs, "Lemma H equivalence failed"
        assert dK == C.det() + tt * sum(u)
    print("H   x>0 (nonsingular)  <=>  adj(C)1 entries nonzero, common sign = sign det K : OK", trials)

def check_two_node():
    e, d = sp.symbols("epsilon delta", positive=True)
    q = d * (1 - e) / e
    # (1+d)(1-e) - 1 = e (q - 1)
    assert sp.simplify((1 + d) * (1 - e) - 1 - e * (q - 1)) == 0
    # 2-node: I - W = [[1,a],[b,1]], eigenvalues of -(I-W) are -1 -/+ sqrt(ab)
    a, b, lam = sp.symbols("a b lam", positive=True)
    cp = sp.expand((lam * sp.eye(2) + sp.Matrix([[1, a], [b, 1]])).det())
    assert sp.expand(cp - ((lam + 1) ** 2 - a * b)) == 0
    print("K2  (1+delta)(1-eps) - 1 = eps (q-1); char poly (lam+1)^2 - ab : OK")


# ---------------------------------------------------------------- H (strengthened) + (2) in beta-form
def check_H_strong(trials=600):
    cnt = 0
    for _ in range(trials):
        n = rng.randint(1, 6); k = rng.randint(1, n)
        A = rand_graph(n); eps, delta = rand_legal()
        sigma = sorted(rng.sample(range(n), k))
        t = (1 - eps) / eps; beta = 1 + delta / eps
        N = [[0 if i == j else 1 - A[sigma[j]][sigma[i]] for j in range(k)] for i in range(k)]
        C = [[F(1 if i == j else 0) + beta * N[i][j] for j in range(k)] for i in range(k)]
        dC = ctln.det_exact(C)
        u = []
        for i in range(k):
            Ci = [row[:] for row in C]
            for r in range(k): Ci[r][i] = F(1)
            u.append(ctln.det_exact(Ci))            # Cramer: (adj C 1)_i
        D = dC + t * sum(u)
        one_signed = all(v > 0 for v in u) or all(v < 0 for v in u)
        if one_signed:
            sg = 1 if u[0] > 0 else -1
            assert dC * sg > 0 and D * sg > 0, "strict one-signed u must force sign det C = sign det K"
        W = ctln.weight_matrix(A, eps, delta)
        x = ctln.solve_exact(I_minus(ctln.submatrix(W, sigma)), [F(1)] * k)
        pos = x is not None and all(v > 0 for v in x)
        assert pos == one_signed, "x>0 iff adj(C)1 strictly one-signed"
        if pos:
            cnt += 1
            for kk in range(n):
                if kk in sigma: continue
                drive = 1 + sum(W[kk][sigma[j]] * x[j] for j in range(k))
                nonin = [j for j in range(k) if A[sigma[j]][kk] == 0]
                assert drive * D == dC - beta * sum(u[j] for j in nonin), "drive formula"
                assert (drive < 0) == (beta * sum(abs(u[j]) for j in nonin) > abs(dC)), "(2) beta-form"
                # outside-in / inside-out domination relative to sigma
                for ii in range(k):
                    i = sigma[ii]
                    if A[i][kk] == 1 and all(A[l][i] == 0 or A[l][kk] == 1 for l in sigma if l != i):
                        assert drive >= eps * x[ii] > 0, "outside-in"
                    if A[i][kk] == 0 and all(A[l][kk] == 0 or A[l][i] == 1 for l in sigma if l != i):
                        assert drive <= -delta * x[ii] < 0, "inside-out"
    print("H+  x>0 <=> adj(C)1 strictly one-signed (then sign det C = sign det K); drive*detK/theta = detC - beta*sum_{j-/->k} u_j;"
          " outside-in => drive >= eps x_i; inside-out => drive <= -delta x_i : OK", cnt, "positive cases")

def check_inside_domination(trials=400):
    errs = [[[0,0,0,0],[0,0,0,0],[1,1,0,0],[1,1,0,0]], [[0,0,0,0],[1,0,0,0],[1,1,0,0],[1,1,1,0]],
            [[0,1,0,0],[1,0,0,0],[1,1,0,1],[1,1,1,0]]]
    import graphs
    for A in errs:
        assert graphs.dominated_nodes(A), "expected a dominated node"
    for _ in range(trials):
        A = [row[:] for row in rng.choice(errs)] if rng.random() < 0.3 else rand_graph(rng.randint(2, 6))
        n = len(A); eps, delta = rand_legal()
        dom = False
        for j in range(n):
            for k in range(n):
                if j != k and A[j][k] == 1 and A[k][j] == 0 and all(A[l][j] == 0 or A[l][k] == 1 for l in range(n) if l not in (j, k)):
                    dom = True
        if dom:
            x = ctln.solve_exact(I_minus(ctln.weight_matrix(A, eps, delta)), [F(1)] * n)
            assert x is None or not all(v > 0 for v in x), "inside domination should kill positivity"
    print("D   graph with a dominated node: no positive full-support solution (incl. the 3 region2d-error 4-node graphs) : OK")

def check_degeneracy_polys_nonzero():
    # at eps = 1, delta = -1 every off-diagonal weight is 0, so I - W = I: all principal minors and
    # all column-replacement determinants equal 1 -> none of these polynomials is identically zero.
    for n in range(1, 7):
        A = rand_graph(n)
        W = ctln.weight_matrix(A, F(1), F(-1))
        assert all(W[i][j] == 0 for i in range(n) for j in range(n))
    print("O   degeneracy polynomials are not identically zero (I - W = I at eps=1, delta=-1) : OK")

if __name__ == "__main__":
    check_H0(); check_F(); check_perm(); check_R2_symbolic(); check_R2_numeric()
    check_G(); check_E(); check_H_symbolic(); check_H_numeric(); check_two_node()
    check_H_strong(); check_inside_domination(); check_degeneracy_polys_nonzero()
    print("ALL SANITY CHECKS PASSED")
