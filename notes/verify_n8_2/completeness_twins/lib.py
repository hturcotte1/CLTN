"""Own, self-contained helpers for agent completeness_twins (verify_n8_2).
Conventions: A[i][j] == 1 iff arc i -> j (CLAUDE.md section 2).
digraph6: '&' + chr(n+63) + n*n bits row by row (bit (i,j) = 1 iff arc i->j), 6 bits/char, char = value+63, MSB first.
Decision engine: exact Fractions only. Hurwitz minors three-valued.
"""
from fractions import Fraction as F
import itertools


# ---------------------------------------------------------------- digraph6
def d6_decode(s):
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    assert len(bits) >= n * n, (s, len(bits), n)
    A = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    return A


def d6_encode(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    out = ['&', chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return ''.join(out)


def arcs(A):
    return sum(map(sum, A))


def underlying_edges(A):
    n = len(A)
    return sum(1 for i in range(n) for j in range(i + 1, n) if A[i][j] or A[j][i])


def permute(A, p):
    """new node i = old node p[i]."""
    n = len(A)
    return [[A[p[i]][p[j]] for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------- own canonical form
def _refine(A, cells):
    """Iteratively refine an ordered partition (list of lists of vertices) by counting, for each vertex,
    its out- and in-neighbours in every cell. Cells are ordered by sorted invariant, deterministic."""
    n = len(A)
    while True:
        cell_of = [0] * n
        for ci, c in enumerate(cells):
            for v in c:
                cell_of[v] = ci
        newcells = []
        changed = False
        for c in cells:
            sig = {}
            for v in c:
                key = tuple((sum(A[v][w] for w in cc), sum(A[w][v] for w in cc)) for cc in cells)
                sig.setdefault(key, []).append(v)
            keys = sorted(sig)
            if len(keys) > 1:
                changed = True
            for k in keys:
                newcells.append(sorted(sig[k]))
        cells = newcells
        if not changed:
            return cells


def canon(A):
    """Canonical digraph6 string: lexicographically smallest adjacency bit-string over all relabellings
    compatible with the refined (ordered) partition; the partition is a relabelling-invariant so the
    minimum over compatible relabellings is a canonical form. Exhaustive within cells (exact)."""
    n = len(A)
    # initial partition: by (outdeg, indeg, bidirected degree)
    inv = {}
    for v in range(n):
        od = sum(A[v]); idg = sum(A[w][v] for w in range(n)); bd = sum(A[v][w] & A[w][v] for w in range(n))
        inv.setdefault((od, idg, bd), []).append(v)
    cells = [sorted(inv[k]) for k in sorted(inv)]
    cells = _refine(A, cells)
    best = None
    # enumerate relabellings: position block for each cell, product of permutations within cells
    perms_per_cell = [list(itertools.permutations(c)) for c in cells]
    for combo in itertools.product(*perms_per_cell):
        p = [v for block in combo for v in block]
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    bits = list(best)
    B = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    return d6_encode(B)


def canon_cost(A):
    n = len(A)
    inv = {}
    for v in range(n):
        od = sum(A[v]); idg = sum(A[w][v] for w in range(n)); bd = sum(A[v][w] & A[w][v] for w in range(n))
        inv.setdefault((od, idg, bd), []).append(v)
    cells = _refine(A, [sorted(inv[k]) for k in sorted(inv)])
    c = 1
    for cl in cells:
        for k in range(2, len(cl) + 1):
            c *= k
    return c


def find_isomorphism(A, B):
    """Explicit permutation p with permute(A, p) == B, or None. Exact (uses refinement of both, then search)."""
    n = len(A)
    if n != len(B):
        return None
    def start(X):
        inv = {}
        for v in range(n):
            od = sum(X[v]); idg = sum(X[w][v] for w in range(n)); bd = sum(X[v][w] & X[w][v] for w in range(n))
            inv.setdefault((od, idg, bd), []).append(v)
        return _refine(X, [sorted(inv[k]) for k in sorted(inv)])
    ca, cb = start(A), start(B)
    if [len(c) for c in ca] != [len(c) for c in cb]:
        return None
    # the refinement signature (which keys) must agree; cheap check via canon
    for combo in itertools.product(*[itertools.permutations(c) for c in ca]):
        pa = [v for block in combo for v in block]
        pb = [v for block in cb for v in block]
        # new node i = old node pa[i] in A, pb[i] in B. Want permute(A, p) == B: B[i][j] = A[p[i]][p[j]].
        # permute(A, pa) == permute(B, pb)  =>  A[pa[i]][pa[j]] = B[pb[i]][pb[j]]  => with p[pb[i]] = pa[i]
        ok = all(A[pa[i]][pa[j]] == B[pb[i]][pb[j]] for i in range(n) for j in range(n))
        if ok:
            p = [0] * n
            for i in range(n):
                p[pb[i]] = pa[i]
            assert permute(A, p) == B
            return p
    return None


# ---------------------------------------------------------------- twins, switching, blow-up
def bidirected_in_twins(A):
    """All pairs u < v with A[u][v] == A[v][u] == 1 and A[k][u] == A[k][v] for every k not in {u, v}."""
    n = len(A)
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            if A[u][v] == 1 and A[v][u] == 1 and all(A[k][u] == A[k][v] for k in range(n) if k not in (u, v)):
                out.append((u, v))
    return out


def in_twin_switches(A, u, v):
    """All single switches: for each w not in {u,v} with A[u][w] + A[v][w] == 1, move the arc to the other twin."""
    n = len(A)
    res = []
    for w in range(n):
        if w in (u, v):
            continue
        if A[u][w] + A[v][w] == 1:
            B = [row[:] for row in A]
            B[u][w], B[v][w] = A[v][w], A[u][w]
            res.append((w, B))
    return res


def blow_up(A, v):
    """Closed-twin blow-up: new vertex n (= v') with the in- and out-neighbourhood of v, plus v <-> v'."""
    n = len(A)
    B = [row[:] + [0] for row in A] + [[0] * (n + 1)]
    for k in range(n):
        B[n][k] = A[v][k]
        B[k][n] = A[k][v]
    B[n][v] = 1; B[v][n] = 1
    B[n][n] = 0
    return B


# ---------------------------------------------------------------- exact CTLN decision engine
def W_matrix(A, eps, delta):
    n = len(A)
    W = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            W[i][j] = F(-1) + eps if A[j][i] == 1 else F(-1) - delta
    return W


def solve(M, b):
    """Exact Gaussian elimination; returns x with M x = b, or None if singular."""
    n = len(M)
    a = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        piv = None
        for r in range(c, n):
            if a[r][c] != 0:
                piv = r; break
        if piv is None:
            return None
        a[c], a[piv] = a[piv], a[c]
        pv = a[c][c]
        a[c] = [x / pv for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return [a[i][n] for i in range(n)]


def charpoly(M):
    """Faddeev-LeVerrier: coefficients of det(lam I - M) = lam^n + c1 lam^(n-1) + ... + cn, exact."""
    n = len(M)
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    Mk = [row[:] for row in I]
    coeffs = [F(1)]
    for k in range(1, n + 1):
        # Mk = M * Mk_prev
        MM = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(MM[i][i] for i in range(n))
        ck = -tr / k
        coeffs.append(ck)
        Mk = [[MM[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
    return coeffs


def det(M):
    n = len(M)
    a = [row[:] for row in M]
    d = F(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if a[r][c] != 0:
                piv = r; break
        if piv is None:
            return F(0)
        if piv != c:
            a[c], a[piv] = a[piv], a[c]; d = -d
        d *= a[c][c]
        for r in range(c + 1, n):
            if a[r][c] != 0:
                f = a[r][c] / a[c][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return d


def hurwitz_minors(coeffs):
    """coeffs = [a0=1, a1, ..., an] of a0 lam^n + ... + an. Returns the n Hurwitz leading principal minors."""
    n = len(coeffs) - 1
    def a(k):
        return coeffs[k] if 0 <= k <= n else F(0)
    H = [[a(2 * j - i + 1) for j in range(n)] for i in range(n)]  # H[i][j] = a_{2j - i + 1}, 0-indexed
    # Standard Hurwitz matrix: row i (1-indexed), col j: a_{2j - i}. 0-indexed: a_{2(j+1) - (i+1)} = a_{2j - i + 1}.
    return [det([row[:k] for row in H[:k]]) for k in range(1, n + 1)]


def decide(A, eps, delta):
    """Three-valued exact decision for the full support. Returns dict with keys
    positive (True/False/None if singular), stability ('stable'|'unstable'|'undecided'), minors, x."""
    n = len(A)
    W = W_matrix(A, eps, delta)
    IW = [[(F(1) if i == j else F(0)) - W[i][j] for j in range(n)] for i in range(n)]
    x = solve(IW, [F(1)] * n)
    positive = None if x is None else all(xi > 0 for xi in x)
    Jm = [[(F(-1) if i == j else F(0)) + W[i][j] for j in range(n)] for i in range(n)]
    cp = charpoly(Jm)
    mins = hurwitz_minors(cp)
    if any(m == 0 for m in mins):
        stab = 'undecided'
    elif all(m > 0 for m in mins):
        stab = 'stable'
    else:
        stab = 'unstable'
    return {'positive': positive, 'stability': stab, 'minors': mins, 'x': x, 'charpoly': cp}


def is_complete(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)


# ---------------------------------------------------------------- exact algebraic endpoints
def poly_eval(coeffs_hi_to_lo, x):
    r = F(0)
    for c in coeffs_hi_to_lo:
        r = r * x + c
    return r


def refine_root(coeffs_hi_to_lo, lo, hi, steps):
    """Bisection on a sign-changing isolating interval [lo, hi] (rationals). Returns (lo, hi) after `steps` halvings.
    If the polynomial has a root exactly at an endpoint (value 0) that endpoint is returned as a point interval."""
    lo, hi = F(lo), F(hi)
    flo, fhi = poly_eval(coeffs_hi_to_lo, lo), poly_eval(coeffs_hi_to_lo, hi)
    if flo == 0:
        return lo, lo
    if fhi == 0:
        return hi, hi
    assert (flo < 0) != (fhi < 0), ('no sign change', lo, hi)
    for _ in range(steps):
        mid = (lo + hi) / 2
        fm = poly_eval(coeffs_hi_to_lo, mid)
        if fm == 0:
            return mid, mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
    return lo, hi


def endpoint_bounds(ep):
    """Return (lower, upper) rational bounds of an interval endpoint record (lo or hi)."""
    if 'beta' in ep:
        b = F(ep['beta']); return b, b
    lo, hi = ep['beta_interval']
    return refine_root(ep['beta_minpoly_hi_to_lo'], F(lo), F(hi), 0)


def interior_point(iv, extra_steps=0):
    """A rational beta strictly inside the stable interval iv (dict with lo, hi), by bisection until
    the endpoint bounds separate. Returns (beta, lo_bounds, hi_bounds, steps)."""
    lo, hi = iv['lo'], iv['hi']
    lob = endpoint_bounds(lo); hib = endpoint_bounds(hi)
    steps = 0
    while not (lob[1] < hib[0]):
        steps += 1
        assert steps < 200
        if 'beta' not in lo:
            lob = refine_root(lo['beta_minpoly_hi_to_lo'], lob[0], lob[1], 1)
        if 'beta' not in hi:
            hib = refine_root(hi['beta_minpoly_hi_to_lo'], hib[0], hib[1], 1)
    for _ in range(extra_steps):
        if 'beta' not in lo:
            lob = refine_root(lo['beta_minpoly_hi_to_lo'], lob[0], lob[1], 1)
        if 'beta' not in hi:
            hib = refine_root(hi['beta_minpoly_hi_to_lo'], hib[0], hib[1], 1)
    beta = (lob[1] + hib[0]) / 2
    assert lob[1] < beta < hib[0]
    return beta, lob, hib, steps


def legal(eps, delta):
    return eps > 0 and delta > 0 and delta * (1 - eps) / eps > 1
