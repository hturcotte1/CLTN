"""Independent exact primitives (Fractions / ints only). No project modules imported."""
from fractions import Fraction as F
from itertools import permutations

def d6_decode(s):
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    need = n * n
    assert len(bits) >= need and len(bits) - need < 6, (s, len(bits))
    assert not any(bits[need:]), "nonzero padding"
    A = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    return A

def d6_encode(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    out = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        out += chr(63 + v)
    return out

def n_arcs(A):
    return sum(A[i][j] for i in range(len(A)) for j in range(len(A)) if i != j)

def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def missing_arcs(A):
    n = len(A)
    return sum(1 for i in range(n) for j in range(n) if i != j and A[i][j] == 0)

def dominates(A, k, j):
    """Curto Def 1.1 as in notes/lemmas.md Lemma D: k dominates j."""
    n = len(A)
    if k == j: return False
    if not (A[j][k] == 1 and A[k][j] == 0): return False
    for i in range(n):
        if i in (j, k): continue
        if A[i][j] == 1 and A[i][k] == 0: return False
    return True

def dominated_pairs(A):
    n = len(A)
    return [(k, j) for k in range(n) for j in range(n) if dominates(A, k, j)]

def canon(A):
    n = len(A)
    best = None
    for p in permutations(range(n)):
        # relabel: new node a = old node p[a]
        s = tuple(A[p[a]][p[b]] for a in range(n) for b in range(n))
        if best is None or s > best: best = s
    return best

# ---------- exact linear algebra ----------
def det(M):
    M = [[F(x) for x in row] for row in M]
    n = len(M); d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n): M[r][k] -= f * M[c][k]
    return d

def solve(M, b):
    n = len(M)
    Ab = [[F(x) for x in M[i]] + [F(b[i])] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if Ab[r][c] != 0), None)
        if p is None: return None  # singular
        Ab[c], Ab[p] = Ab[p], Ab[c]
        for r in range(n):
            if r != c and Ab[r][c] != 0:
                f = Ab[r][c] / Ab[c][c]
                for k in range(c, n + 1): Ab[r][k] -= f * Ab[c][k]
    return [Ab[i][n] / Ab[i][i] for i in range(n)]

def matmul(X, Y):
    n = len(X); m = len(Y[0]); k = len(Y)
    return [[sum(X[i][t] * Y[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def charpoly(M):
    """det(x I - M), coefficients high->low, Faddeev-LeVerrier, exact."""
    n = len(M)
    M = [[F(x) for x in row] for row in M]
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    c = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]  # M_0 = 0
    for k in range(1, n + 1):
        # M_k = M (M_{k-1} + c_{k-1} I)
        T = [[Mk[i][j] + (c[-1] if i == j else 0) for j in range(n)] for i in range(n)]
        Mk = matmul(M, T)
        ck = -sum(Mk[i][i] for i in range(n)) / k
        c.append(ck)
    return c

def poly_sub(p, q):
    L = max(len(p), len(q))
    p = [F(0)] * (L - len(p)) + list(p); q = [F(0)] * (L - len(q)) + list(q)
    r = [a - b for a, b in zip(p, q)]
    while len(r) > 1 and r[0] == 0: r.pop(0)
    return r

def poly_eval(p, x):
    v = F(0)
    for a in p: v = v * x + a
    return v

def poly_shift(p, s):
    """return coefficients of p(x + s), high->low (Horner/Taylor)."""
    n = len(p) - 1
    a = [F(x) for x in p]
    # repeated synthetic division
    for i in range(n):
        for j in range(1, n - i + 1):
            a[j] += s * a[j - 1]
    return a

def hurwitz(p):
    """p high->low with p[0] > 0.  Returns ('stable'|'unstable'|'undecided', minors)."""
    p = [F(x) for x in p]
    assert p[0] != 0
    if p[0] < 0: p = [-x for x in p]
    n = len(p) - 1
    a = p  # a[0] x^n + a[1] x^{n-1} + ...
    def coef(k):
        return a[k] if 0 <= k <= n else F(0)
    H = [[coef(2 * (j + 1) - (i + 1)) for j in range(n)] for i in range(n)]
    minors = [det([row[:k] for row in H[:k]]) for k in range(1, n + 1)]
    if any(m < 0 for m in minors):
        # a negative minor while all prior positive is unstable; any negative minor rules out Hurwitz anyway
        first = next(i for i, m in enumerate(minors) if m <= 0)
        if minors[first] == 0: return 'undecided', minors
        return 'unstable', minors
    if any(m == 0 for m in minors): return 'undecided', minors
    return 'stable', minors
