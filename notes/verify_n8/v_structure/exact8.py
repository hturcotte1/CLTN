"""Independent exact primitives for the n = 8 structural verifier (v_structure).
Fractions / ints only.  No project module (ctln, symbolic, fastregion, graphs, n7/) is imported.
"""
from fractions import Fraction as F
from itertools import permutations

# ---------------------------------------------------------------- digraph6
def d6_decode(s):
    """'&' + N(n) + n*n bits row by row (bit (i,j) = 1 iff arc i->j), 6 bits/char, char = value+63, MSB first."""
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
    assert not any(bits[need:]), ("nonzero padding", s)
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

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

# ---------------------------------------------------------------- graph basics
def n_arcs(A):
    return sum(A[i][j] for i in range(len(A)) for j in range(len(A)) if i != j)

def missing_arcs(A):
    n = len(A)
    return sum(1 for i in range(n) for j in range(n) if i != j and A[i][j] == 0)

def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def graphically_dominates(A, k, j):
    """Curto arXiv:2510.05098 Definition 1.1, implemented verbatim:
    k graphically dominates j (k > j) iff
      (i)  for each vertex i in [n] \\ {j, k}: if i -> j then i -> k;
      (ii) j -> k and k -/-> j.
    A[i][j] == 1 iff arc i -> j."""
    n = len(A)
    if k == j: return False
    if not (A[j][k] == 1 and A[k][j] == 0): return False          # (ii)
    for i in range(n):                                            # (i)
        if i == j or i == k: continue
        if A[i][j] == 1 and A[i][k] == 0: return False
    return True

def dominated_pairs(A):
    """list of (k, j) with k graphically dominating j."""
    n = len(A)
    return [(k, j) for k in range(n) for j in range(n) if graphically_dominates(A, k, j)]

def is_domination_free(A):
    """no dominated node, i.e. no pair (k, j) with k > j."""
    return not dominated_pairs(A)

# ---------------------------------------------------------------- canonical form (individualisation / refinement)
def _refine(A, cells):
    """Equitable refinement of the ordered partition `cells` (list of lists).  Isomorphism-invariant:
    a vertex's signature is (its cell index, for every cell c: (#out-neighbours in c, #in-neighbours in c));
    each cell is split by sorting signatures; repeat until no cell splits."""
    n = len(A)
    while True:
        cell_of = [0] * n
        for ci, c in enumerate(cells):
            for v in c: cell_of[v] = ci
        new = []
        changed = False
        for ci, c in enumerate(cells):
            if len(c) == 1:
                new.append(c); continue
            sig = {}
            for v in c:
                out = [0] * len(cells); inn = [0] * len(cells)
                for w in range(n):
                    if w == v: continue
                    if A[v][w]: out[cell_of[w]] += 1
                    if A[w][v]: inn[cell_of[w]] += 1
                sig[v] = tuple(zip(out, inn))
            keys = sorted(set(sig.values()))
            if len(keys) == 1:
                new.append(c); continue
            changed = True
            for k in keys:
                new.append([v for v in c if sig[v] == k])
        cells = new
        if not changed:
            return cells

def _leaf_code(A, order):
    """adjacency of the relabelled graph (new vertex a = old vertex order[a]) as a tuple of bits."""
    n = len(A)
    return tuple(A[order[a]][order[b]] for a in range(n) for b in range(n))

def canon(A):
    """Canonical form: the minimum leaf code over the complete individualisation-refinement search tree
    (no automorphism pruning).  Equal canon <=> isomorphic: '=>' because each leaf code is a relabelling of
    the graph; '<=' because refinement and target-cell choice (first non-singleton cell of smallest size)
    depend only on the ordered partition structure, which is invariant under relabelling."""
    n = len(A)
    best = [None]
    def search(cells):
        cells = _refine(A, cells)
        if all(len(c) == 1 for c in cells):
            code = _leaf_code(A, [c[0] for c in cells])
            if best[0] is None or code < best[0]: best[0] = code
            return
        # target cell: first cell of minimum size among non-singleton cells
        sizes = [len(c) for c in cells]
        m = min(s for s in sizes if s > 1)
        ti = next(i for i, s in enumerate(sizes) if s == m)
        for v in cells[ti]:
            rest = [w for w in cells[ti] if w != v]
            search(cells[:ti] + [[v], rest] + cells[ti + 1:])
    search([list(range(n))])
    return best[0]

def canon_bruteforce(A):
    """minimum code over all n! relabellings (only for self-tests / small samples)."""
    n = len(A)
    best = None
    for p in permutations(range(n)):
        s = _leaf_code(A, p)
        if best is None or s < best: best = s
    return best

def code_to_matrix(code, n):
    return [list(code[i * n:(i + 1) * n]) for i in range(n)]

# ---------------------------------------------------------------- exact linear algebra
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
    c = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        T = [[Mk[i][j] + (c[-1] if i == j else 0) for j in range(n)] for i in range(n)]
        Mk = matmul(M, T)
        c.append(-sum(Mk[i][i] for i in range(n)) / k)
    return c

def hurwitz(p):
    """p high->low, p[0] != 0.  Routh-Hurwitz by the leading principal minors of the Hurwitz matrix.
    Returns ('stable' | 'unstable' | 'undecided', minors).  A vanishing minor (before any negative one) is UNDECIDED."""
    p = [F(x) for x in p]
    assert p[0] != 0
    if p[0] < 0: p = [-x for x in p]
    n = len(p) - 1
    def coef(k): return p[k] if 0 <= k <= n else F(0)
    H = [[coef(2 * (j + 1) - (i + 1)) for j in range(n)] for i in range(n)]
    minors = [det([row[:k] for row in H[:k]]) for k in range(1, n + 1)]
    for m in minors:
        if m == 0: return 'undecided', minors
        if m < 0: return 'unstable', minors
    return 'stable', minors

# ---------------------------------------------------------------- CTLN at a finite point
def W_of(A, eps, delta):
    n = len(A)
    return [[F(0) if i == j else ((-1 + eps) if A[j][i] == 1 else (-1 - delta)) for j in range(n)] for i in range(n)]

def ctln_check(A, eps, delta):
    """full-support check at exact (eps, delta): legality, positivity of the fixed point, Hurwitz verdict."""
    eps = F(eps); delta = F(delta); n = len(A)
    q = delta * (1 - eps) / eps
    legal = (0 < eps < 1) and delta > 0 and q > 1
    W = W_of(A, eps, delta)
    IW = [[F(int(i == j)) - W[i][j] for j in range(n)] for i in range(n)]
    # closed form cross-check  I - W = -delta I + (1+delta) J - (eps+delta) A^T
    CF = [[-delta * int(i == j) + (1 + delta) - (eps + delta) * A[j][i] for j in range(n)] for i in range(n)]
    assert IW == CF
    x = solve(IW, [1] * n)
    positive = x is not None and all(v > 0 for v in x)
    Jac = [[W[i][j] - int(i == j) for j in range(n)] for i in range(n)]
    st, minors = hurwitz(charpoly(Jac))
    return {'eps': str(eps), 'delta': str(delta), 'q': str(q), 'legal': legal, 'singular': x is None,
            'positive': positive, 'stable': st, 'min_x': str(min(x)) if x else None}

# ---------------------------------------------------------------- polynomials (high -> low)
def poly_eval(p, x):
    v = F(0)
    for a in p: v = v * x + a
    return v

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[0] == 0: p.pop(0)
    return p

def poly_rem(p, d):
    p = [F(x) for x in p]; d = [F(x) for x in d]
    while len(p) >= len(d) and any(p):
        c = p[0] / d[0]
        for k in range(len(d)): p[k] -= c * d[k]
        p.pop(0)
    return poly_trim(p) if p else [F(0)]

def poly_deriv(p):
    n = len(p) - 1
    return [F(p[i]) * (n - i) for i in range(n)] if n > 0 else [F(0)]

def sturm_sequence(p):
    seq = [poly_trim([F(x) for x in p]), poly_trim(poly_deriv(p))]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        r = poly_rem(seq[-2], seq[-1])
        r = [-x for x in r]
        if len(r) == 1 and r[0] == 0: break
        seq.append(r)
    return seq

def _sign_changes(seq, x):
    vals = [poly_eval(s, x) for s in seq]
    vals = [v for v in vals if v != 0]
    return sum(1 for a, b in zip(vals, vals[1:]) if (a < 0) != (b < 0))

def sturm_count(seq, a, b):
    """number of distinct real roots in (a, b]  (a < b)."""
    return _sign_changes(seq, F(a)) - _sign_changes(seq, F(b))

# ---------------------------------------------------------------- self tests
if __name__ == '__main__':
    import random
    assert hurwitz([1, 3, 3, 1])[0] == 'stable'
    assert hurwitz([1, 0, 1])[0] == 'undecided'
    assert hurwitz([1, -1, 1])[0] == 'unstable'
    assert hurwitz([1, 1, 1, 6])[0] == 'unstable'       # (x+2)(x^2-x+3)
    assert hurwitz([1, 2, 2, 1])[0] == 'stable'          # (x+1)(x^2+x+1)
    assert hurwitz([1, 1, 1, 1])[0] == 'undecided'       # (x+1)(x^2+1): imaginary pair
    assert charpoly([[0, 1], [1, 0]]) == [1, 0, -1]
    for s in ['&EZjtwWE', '&FNFL~[udo?', '&GM?WUZeEie^O']:
        assert d6_encode(d6_decode(s)) == s
    G = d6_decode('&EZjtwWE')
    assert len(G) == 6 and n_arcs(G) == 19
    # canonical form: invariance under random relabelling + agreement with brute force (n = 5, 6, 7)
    random.seed(1)
    for trial in range(60):
        n = random.choice([5, 6, 7]); dens = random.random()
        A = [[0 if i == j else int(random.random() < dens) for j in range(n)] for i in range(n)]
        p = list(range(n)); random.shuffle(p)
        B = [[A[p[a]][p[b]] for b in range(n)] for a in range(n)]
        assert canon(A) == canon(B)
        assert canon(A) == canon(code_to_matrix(canon(A), n))
        if n <= 6:
            # equality of my codes must agree with equality of brute-force codes (iso test), on a second random graph
            C = [[0 if i == j else int(random.random() < dens) for j in range(n)] for i in range(n)]
            assert (canon(A) == canon(C)) == (canon_bruteforce(A) == canon_bruteforce(C))
            assert (canon(B) == canon(C)) == (canon_bruteforce(B) == canon_bruteforce(C))
    # domination definition on a hand example: 3 nodes, 0->1, 1->2, 0->2 ... k=2 dominates j=1? (ii) 1->2 and 2-/->1 ok; (i) 0->1 => 0->2 ok
    A = [[0, 1, 1], [0, 0, 1], [0, 0, 0]]
    assert graphically_dominates(A, 2, 1) and graphically_dominates(A, 1, 0) and graphically_dominates(A, 2, 0)
    A = [[0, 1, 1], [0, 0, 1], [0, 1, 0]]   # now 2 -> 1 so (ii) fails for (2,1)
    assert not graphically_dominates(A, 2, 1)
    assert is_domination_free(d6_decode('&EZjtwWE'))      # Geneson's graph is domination free (known)
    # Sturm: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    seq = sturm_sequence([1, -6, 11, -6])
    assert sturm_count(seq, 0, 4) == 3 and sturm_count(seq, F(3, 2), F(5, 2)) == 1 and sturm_count(seq, 1, 3) == 2
    # Geneson at the paper's parameters is a counterexample (eps = 1/4? no: use CLAUDE.md q ~ 29/25 with eps=1/2 -> delta = 29/25*... )
    print('self tests ok')
