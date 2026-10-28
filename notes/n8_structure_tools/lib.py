"""Exact structural helpers for the n = 8 structure analysis (notes/n8_structure.md).

Conventions: CLAUDE.md section 2.  A[i][j] = 1 iff i -> j.  D = J - I - A is the
missing-arc ("non-edge") digraph: D[i][j] = 1 iff i != j and i -/-> j.
Everything here is integer / Fraction arithmetic.  No floating point anywhere.
"""
import os, sys, json, subprocess, tempfile
from fractions import Fraction as F
from itertools import combinations

REPO = '/home/user/CLTN'
LABELG = REPO + '/n7/labelg'


# ------------------------------------------------------------------ digraph6
def from_d6(s):
    """Own digraph6 decoder ('&' + N(n) + n*n bits row-major, bit (i,j) = arc i->j)."""
    assert s[0] == '&'
    if ord(s[1]) != 126:
        # standard for n <= 62; ord > 126 is the single-byte header written by notes/research_wf3/theory/tlib.to_d6
        # for n >= 64 (non-standard, e.g. the n = 73 Theorem Y witness): accepted, length-checked below
        n = ord(s[1]) - 63; p = s[2:]
    else:                                  # N(n) for 63 <= n <= 258047: 126 + three 6-bit bytes
        assert ord(s[2]) < 126
        n = ((ord(s[2]) - 63) << 12) | ((ord(s[3]) - 63) << 6) | (ord(s[4]) - 63); p = s[5:]
    assert n <= 1000 and len(p) == (n * n + 5) // 6, 'bad digraph6 length'
    A = [[0] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63
            A[i][j] = (c >> (5 - bit % 6)) & 1
            bit += 1
    for i in range(n):
        assert A[i][i] == 0, 'loop in ' + s
    return A


def to_d6(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    if n <= 62:
        s = '&' + chr(n + 63)
    else:
        s = '&' + chr(126) + chr((n >> 12) + 63) + chr(((n >> 6) & 63) + 63) + chr((n & 63) + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + b
        s += chr(v + 63)
    return s


def labelg_canon(d6_list):
    """nauty labelg canonical labels (exact). Returns list aligned with input."""
    if not d6_list:
        return []
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, 'in.d6'), os.path.join(td, 'out.d6')
        with open(fi, 'w') as f:
            f.write('\n'.join(d6_list) + '\n')
        subprocess.run([LABELG, '-q', fi, fo], check=True)
        out = open(fo).read().split()
    assert len(out) == len(d6_list)
    return out


def induced(A, nodes):
    return [[A[i][j] for j in nodes] for i in nodes]


# ------------------------------------------------------------------ basic counts
def arcs(A):
    n = len(A)
    return sum(A[i][j] for i in range(n) for j in range(n) if i != j)


def pair_types(A):
    """(bidirected pairs, single-arc pairs, non-adjacent pairs)."""
    n = len(A)
    b = s = z = 0
    for i, j in combinations(range(n), 2):
        t = A[i][j] + A[j][i]
        if t == 2: b += 1
        elif t == 1: s += 1
        else: z += 1
    return b, s, z


def nonedge(A):
    n = len(A)
    return [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]


def outdeg(M, i):
    return sum(M[i][j] for j in range(len(M)) if j != i)


def indeg(M, j):
    return sum(M[i][j] for i in range(len(M)) if i != j)


def has_directed_cycle(M):
    n = len(M)
    indg = [indeg(M, j) for j in range(n)]
    q = [j for j in range(n) if indg[j] == 0]
    seen = 0
    while q:
        v = q.pop(); seen += 1
        for w in range(n):
            if w != v and M[v][w]:
                indg[w] -= 1
                if indg[w] == 0:
                    q.append(w)
    return seen < n


def strongly_connected(M):
    n = len(M)
    def reach(adj):
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for w in range(n):
                if w != v and adj(v, w) and w not in seen:
                    seen.add(w); st.append(w)
        return len(seen) == n
    return reach(lambda v, w: M[v][w]) and reach(lambda v, w: M[w][v])


def underlying_edges(M):
    n = len(M)
    return [(i, j) for i, j in combinations(range(n), 2) if M[i][j] or M[j][i]]


def components(n, edges):
    par = list(range(n))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i, j in edges:
        par[f(i)] = f(j)
    comp = {}
    for v in range(n):
        comp.setdefault(f(v), []).append(v)
    return list(comp.values())


def undirected_shape(n, edges):
    """Classify a simple undirected graph (given by its edge list on vertices 0..n-1)."""
    deg = [0] * n
    for i, j in edges:
        deg[i] += 1; deg[j] += 1
    nontriv = [c for c in components(n, edges) if len(c) > 1]
    is_matching = max(deg) <= 1
    is_star = len(nontriv) == 1 and any(deg[v] == len(edges) for v in nontriv[0]) if edges else True
    def comp_is_star(c):
        ce = [(i, j) for i, j in edges if i in c]
        return any(deg[v] == len(ce) for v in c) and len(ce) == len(c) - 1
    is_star_forest = all(comp_is_star(c) for c in nontriv)
    is_forest = len(edges) == sum(len(c) - 1 for c in nontriv)
    return dict(matching=is_matching, star=is_star, star_forest=is_star_forest, forest=is_forest,
                n_components_nontrivial=len(nontriv), max_degree=max(deg), covered=sum(1 for d in deg if d))


def bidirected_clique_number(A):
    n = len(A)
    best = 0; arg = []
    for mask in range(1 << n):
        nodes = [i for i in range(n) if mask >> i & 1]
        if len(nodes) <= best:
            continue
        if all(A[i][j] and A[j][i] for i, j in combinations(nodes, 2)):
            best = len(nodes); arg = nodes
    return best, arg


def in_set(A, v, excl=()):
    return frozenset(i for i in range(len(A)) if i != v and i not in excl and A[i][v])


def out_set(A, v, excl=()):
    return frozenset(j for j in range(len(A)) if j != v and j not in excl and A[v][j])


def twin_pairs(A):
    """Counts of bidirected in-twins, closed twins (true twins), non-adjacent in-twins, bidirected out-twins."""
    n = len(A)
    res = dict(bi_in=[], closed=[], nonadj_in=[], bi_out=[])
    for u, v in combinations(range(n), 2):
        bi = A[u][v] and A[v][u]
        non = not A[u][v] and not A[v][u]
        same_in = in_set(A, u, (v,)) == in_set(A, v, (u,))
        same_out = out_set(A, u, (v,)) == out_set(A, v, (u,))
        if bi and same_in: res['bi_in'].append((u, v))
        if bi and same_in and same_out: res['closed'].append((u, v))
        if non and same_in: res['nonadj_in'].append((u, v))
        if bi and same_out: res['bi_out'].append((u, v))
    return res


def closed_twin_classes(A):
    n = len(A)
    cl = closed_twin_pairs = twin_pairs(A)['closed']
    return sorted(components(n, cl), key=lambda c: (-len(c), c))


# ------------------------------------------------------------------ domination (Curto Def 1.1, as used by the pipeline)
def dominated_nodes(A):
    """Graphical domination (Curto et al.): k dominates j (k != j) iff
    in(j) \\ {k} is contained in in(k) \\ {j}... the pipeline's Lemma D form is taken from graphs.py;
    here we import it so the statement tested is the production one."""
    sys.path.insert(0, REPO)
    import graphs
    return graphs.dominated_nodes(A)


# ------------------------------------------------------------------ exact linear algebra
def charpoly_int(M):
    """det(zI - M), hi -> lo, exact (Faddeev-LeVerrier with Fractions)."""
    n = len(M)
    Mf = [[F(x) for x in row] for row in M]
    c = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    for k in range(1, n + 1):
        # Mk = M (Mk-1 + c_{k-1} I)
        B = [[Mk[i][j] + (c[-1] if i == j else 0) for j in range(n)] for i in range(n)]
        Mk = [[sum(Mf[i][l] * B[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        ck = -sum(Mk[i][i] for i in range(n)) / k
        c.append(ck)
    assert all(x.denominator == 1 for x in c)
    return [int(x) for x in c]


def inwalk_vectors(A):
    """For every vertex i the tuple ((1^T A^k)_i for k = 1..n-1)."""
    n = len(A)
    row = [1] * n
    out = [[] for _ in range(n)]
    for k in range(1, n):
        row = [sum(row[l] * A[l][j] for l in range(n)) for j in range(n)]
        for j in range(n):
            out[j].append(row[j])
    return [tuple(v) for v in out]


def D_class_key(A):
    """Lemma J.2: D(A) = D(B) iff after relabelling cospectral with equal per-vertex in-walk counts."""
    return (tuple(charpoly_int(A)), tuple(sorted(inwalk_vectors(A))))


def solve(M, b):
    """Exact Gaussian elimination over Fractions; returns None if singular."""
    n = len(M)
    a = [[F(x) for x in row] + [F(y)] for row, y in zip(M, b)]
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None:
            return None
        a[c], a[p] = a[p], a[c]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c] / a[c][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return [a[i][n] / a[i][i] for i in range(n)]


def limit_fixed_point(A, beta):
    """w = (A^T + gamma I)^{-1} 1, gamma = 1 - 1/beta (Lemma D2: the eps -> 0 positivity vector, t-free)."""
    n = len(A)
    g = 1 - F(1) / F(beta)
    M = [[F(A[j][i]) + (g if i == j else 0) for j in range(n)] for i in range(n)]
    return solve(M, [1] * n)


def automorphism_count(A):
    """|Aut(A)| exactly: backtracking over vertex maps respecting (outdeg, indeg) classes and adjacency."""
    n = len(A)
    sig = [(outdeg(A, i), indeg(A, i), sorted((outdeg(A, j), indeg(A, j)) for j in range(n) if A[i][j]),
            sorted((outdeg(A, j), indeg(A, j)) for j in range(n) if A[j][i])) for i in range(n)]
    cnt = 0
    img = [-1] * n
    used = [False] * n
    def bt(i):
        nonlocal cnt
        if i == n:
            cnt += 1; return
        for c in range(n):
            if used[c] or sig[c] != sig[i]:
                continue
            ok = all(A[i][k] == A[c][img[k]] and A[k][i] == A[img[k]][c] for k in range(i))
            if ok:
                img[i] = c; used[c] = True
                bt(i + 1)
                used[c] = False
        img[i] = -1
    bt(0)
    return cnt


def load_hits(path=REPO + '/n8/n8_limit_HITS.jsonl'):
    return [json.loads(l) for l in open(path)]
