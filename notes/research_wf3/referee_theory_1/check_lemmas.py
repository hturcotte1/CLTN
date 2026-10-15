"""Independent exact checks of Lemma W (W1,W2,W3), X1, X2, X1', X2'.
Built from scratch: graphs constructed here, actual CTLN data from ctln.weight_matrix / analyze_support.
No floating point anywhere."""
import sys, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
import ctln
z = sp.symbols('z')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 12345)

def rand_graph(m, p=None):
    p = p if p is not None else random.choice([F(1,3), F(1,2), F(2,3), F(3,4)])
    return [[0 if i == j else (1 if random.random() < p else 0) for j in range(m)] for i in range(m)]

def nonedge(A):
    m = len(A); return sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - A[i][j])

def chi(Mm):
    return Mm.charpoly(z).as_expr()

def uadjv(X, u, v):
    """v^T adj(zI - X) u  via the matrix determinant lemma: det(zI - X + u v^T) - det(zI - X)."""
    return sp.expand(chi(X - u * v.T) - chi(X))

def psi(Dm):
    n = Dm.shape[0]
    return uadjv(Dm, sp.ones(n, 1), sp.ones(n, 1))

def blowup(H, mu):
    idx = [a for a in range(len(H)) for _ in range(mu[a])]
    n = len(idx)
    return [[0 if i == j else (1 if idx[i] == idx[j] else H[idx[i]][idx[j]]) for j in range(n)] for i in range(n)], idx

def ctln_x(A, eps, delta):
    n = len(A)
    W = ctln.weight_matrix(A, eps, delta)
    IW = [[(F(1) if i == j else F(0)) - W[i][j] for j in range(n)] for i in range(n)]
    return ctln.solve_exact(IW, [F(1)] * n)

def rand_params():
    # random legal or illegal rational (eps, delta); positivity is t-free so legality not needed for W3/X2 checks
    eps = F(1, random.randint(2, 40)); delta = F(random.randint(1, 60), random.randint(1, 20))
    return eps, delta

fails = []
ncheck = {'W1': 0, 'W2a': 0, 'W2b': 0, 'W3': 0, 'X1': 0, 'X2': 0, 'X2pos_hits': 0, 'X1p': 0, 'X2p': 0, 'X2p_hits': 0}

# ---------- Lemma W ----------
for trial in range(40):
    m = random.randint(2, 5)
    H = rand_graph(m)
    mu = [random.randint(1, 3) for _ in range(m)]
    B, idx = blowup(H, mu); n = len(B)
    D = nonedge(H); Dp = nonedge(B); M = sp.diag(*mu)
    # W1
    lhs = psi(Dp)
    rhs = sp.expand(z ** (n - m) * uadjv(D * M, sp.ones(m, 1), sp.Matrix(mu)))
    ncheck['W1'] += 1
    if sp.expand(lhs - rhs) != 0: fails.append(('W1', H, mu))
    # also check the claim D' = E D E^T and A' = E(A+I)E^T - I by construction
    E = sp.Matrix(n, m, lambda i, a: 1 if idx[i] == a else 0)
    Am = sp.Matrix(H); Ap = sp.Matrix(B)
    assert Dp == E * D * E.T and Ap == E * (Am + sp.eye(m)) * E.T - sp.eye(n)
    # W2: det(zI + A') = (z-1)^(n-m) det(zI + (A+I)M - I)
    AQ = (Am + sp.eye(m)) * M - sp.eye(m)
    l2 = chi(-Ap)
    r2 = sp.expand((z - 1) ** (n - m) * chi(-AQ))
    ncheck['W2a'] += 1
    if sp.expand(l2 - r2) != 0: fails.append(('W2a', H, mu))
    l3 = uadjv(-Ap, sp.ones(n, 1), sp.ones(n, 1))
    r3 = sp.expand((z - 1) ** (n - m) * uadjv(-AQ, sp.ones(m, 1), sp.Matrix(mu)))
    ncheck['W2b'] += 1
    if sp.expand(l3 - r3) != 0: fails.append(('W2b', H, mu))
    # W3 against the ACTUAL CTLN fixed point (x proportional to w', class-constant)
    for _ in range(3):
        eps, delta = rand_params()
        beta = 1 + delta / eps; gamma = 1 - 1 / beta
        x = ctln_x(B, eps, delta)
        Y = ((Am.T + sp.eye(m)) * M + (sp.Rational(gamma) - 1) * sp.eye(m))
        if Y.det() == 0 or x is None:
            continue
        y = Y.LUsolve(sp.ones(m, 1))
        wp = [y[idx[i]] for i in range(n)]
        # w' must solve (A'^T + gamma I) w' = 1 exactly
        Pp = Ap.T + sp.Rational(gamma) * sp.eye(n)
        ok1 = (Pp * sp.Matrix(wp) - sp.ones(n, 1)) == sp.zeros(n, 1)
        # and the CTLN x must be a scalar multiple of w' (class-constant)
        ratios = set(sp.Rational(x[i]) / wp[i] for i in range(n) if wp[i] != 0)
        zeros_ok = all((wp[i] == 0) == (x[i] == 0) for i in range(n))
        ok2 = len(ratios) == 1 and zeros_ok
        # sign: x>0 iff w'>0
        ok3 = (all(v > 0 for v in x) == all(v > 0 for v in wp))
        ncheck['W3'] += 1
        if not (ok1 and ok2 and ok3): fails.append(('W3', H, mu, str(eps), str(delta), ok1, ok2, ok3))

# ---------- X1 / X2 ----------
def ext(A, S, T, k):
    m = len(A); n = m + k
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
        for v in range(m, n): B[i][v] = S[i]; B[v][i] = T[i]
    for a in range(m, n):
        for b in range(m, n):
            if a != b: B[a][b] = 1
    return B

def w_vec(A, gamma):
    m = len(A); P = sp.Matrix(A).T + sp.Rational(gamma) * sp.eye(m)
    if P.det() == 0: return None
    return P.LUsolve(sp.ones(m, 1))

for trial in range(40):
    m = random.randint(2, 6)
    A = rand_graph(m)
    S = [random.randint(0, 1) for _ in range(m)]
    T = [random.randint(0, 1) for _ in range(m)]
    D = nonedge(A)
    b = sp.Matrix([1 - s for s in S]); c = sp.Matrix([1 - t for t in T])
    L = D - b * sp.ones(1, m) - sp.ones(m, 1) * c.T
    for k in range(1, 5):
        B = ext(A, S, T, k); n = m + k
        lhs = psi(nonedge(B))
        rhs = sp.expand(z ** (k - 1) * (z * psi(D) + k * chi(L)))
        ncheck['X1'] += 1
        if sp.expand(lhs - rhs) != 0: fails.append(('X1', A, S, T, k))

# X2: T = all. Compare criterion against ctln positivity of full support of explicit G[S,all;k]
for trial in range(200):
    m = random.randint(2, 6)
    A = rand_graph(m, random.choice([F(1,2), F(2,3), F(3,4), F(4,5)]))
    S = [random.randint(0, 1) for _ in range(m)]
    T = [1] * m
    for _ in range(2):
        eps, delta = rand_params()
        beta = 1 + delta / eps; gamma = 1 - 1 / beta
        w = w_vec(A, gamma)
        if w is None:
            crit = None
        else:
            sig = sum((w[i] for i in range(m) if S[i]), sp.Integer(0))
            crit = all(v > 0 for v in w) and sig > 1
        for k in (1, 2, 3, 5):
            B = ext(A, S, T, k)
            rec = ctln.analyze_support(B, list(range(m + k)), eps, delta)
            if rec['x'] is None:
                continue
            ncheck['X2'] += 1
            pos = rec['positive']
            if pos: ncheck['X2pos_hits'] += 1
            if crit is None:
                if pos: fails.append(('X2-undefined-wG-but-positive', A, S, k, str(eps), str(delta)))
            elif pos != crit:
                fails.append(('X2', A, S, k, str(eps), str(delta), pos, crit))

# ---------- X1' / X2' ----------
def ext_multi(A, classes):
    m = len(A); sizes = [kk for S, kk in classes]; n = m + sum(sizes)
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    owner = []
    for cidx, (S, kk) in enumerate(classes): owner += [cidx] * kk
    for v in range(m, n):
        S = classes[owner[v - m]][0]
        for i in range(m): B[i][v] = S[i]; B[v][i] = 1
        for w in range(m, n):
            if w != v: B[v][w] = 1
    return B

for trial in range(25):
    m = random.randint(2, 5)
    A = rand_graph(m)
    D = nonedge(A)
    nc = random.randint(2, 3)
    classes = [([random.randint(0, 1) for _ in range(m)], random.randint(1, 3)) for _ in range(nc)]
    B = ext_multi(A, classes); n = len(B)
    p = sum(kk for S, kk in classes)
    Bv = sp.zeros(m, 1)
    for S, kk in classes:
        for i in range(m): Bv[i] += kk * (1 - S[i])
    rhs = sp.expand(z ** (p - 1) * (z * psi(D) + p * chi(D) + uadjv(D, Bv, sp.ones(m, 1))))
    ncheck['X1p'] += 1
    if sp.expand(psi(nonedge(B)) - rhs) != 0: fails.append(('X1p', A, classes))

for trial in range(150):
    m = random.randint(2, 5)
    A = rand_graph(m, random.choice([F(2,3), F(3,4), F(4,5)]))
    nc = random.randint(2, 3)
    classes = [([random.choice([0, 1, 1]) for _ in range(m)], random.randint(1, 3)) for _ in range(nc)]
    B = ext_multi(A, classes); n = len(B)
    for _ in range(2):
        eps, delta = rand_params()
        beta = 1 + delta / eps; gamma = 1 - 1 / beta
        w = w_vec(A, gamma)
        rec = ctln.analyze_support(B, list(range(n)), eps, delta)
        if rec['x'] is None: continue
        ncheck['X2p'] += 1
        pos = rec['positive']
        if pos: ncheck['X2p_hits'] += 1
        if w is None:
            if pos: fails.append(('X2p-undef', A, classes, str(eps), str(delta)))
            continue
        crit = all(v > 0 for v in w) and all(sum((w[i] for i in range(m) if S[i]), sp.Integer(0)) > 1 for S, kk in classes)
        if pos != crit: fails.append(('X2p', A, classes, str(eps), str(delta), pos, crit))

print('checks:', ncheck)
print('FAILURES:', len(fails))
for f in fails[:20]: print(f)
