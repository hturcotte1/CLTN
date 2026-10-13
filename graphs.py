"""
graphs.py -- digraph enumeration up to isomorphism, with canonical forms.

VALIDATION ANCHOR: the number of directed graphs on n nodes up to isomorphism
is OEIS A000273:  n = 0..6  ->  1, 1, 3, 16, 218, 9608, 1540944

The n=5 entry (9608) is the same count that appears in Curto's published
`n5-graphs-package`, so matching it is a genuine cross-check against the
literature, not just against ourselves. If enumerate_digraphs(5) does not
return 9608 canonical graphs, STOP -- every downstream count is wrong.

Production note: for n=6 (1.5M graphs) pure Python augmentation is slow but
feasible. If `nauty` is available (apt-get install nauty, or build from
https://pallini.di.uniroma1.it/), prefer `directg`/`geng` and use the
canonical_form() here only to cross-check a sample. Never swap in a faster
enumerator without re-running the A000273 check.
"""

from itertools import permutations, combinations

# OEIS A000273 -- digraphs on n nodes up to isomorphism
A000273 = [1, 1, 3, 16, 218, 9608, 1540944]


# ------------------------------------------------------------------ encoding
# A digraph on n nodes is a tuple of n tuples (row-major adjacency), where
# A[i][j] == 1 iff  i -> j.  Same convention as ctln.py.  Do not change it.


def to_key(A):
    """Flatten to an immutable comparable key (off-diagonal entries only)."""
    n = len(A)
    return tuple(A[i][j] for i in range(n) for j in range(n) if i != j)


def permute(A, p):
    """Relabel: new node i is old node p[i]."""
    n = len(A)
    return tuple(tuple(A[p[i]][p[j]] for j in range(n)) for i in range(n))


# ------------------------------------------------------------------ refinement


def _refine(A):
    """1-dimensional Weisfeiler-Leman colour refinement.

    Returns a list of cells (lists of vertices) partitioning 0..n-1, ordered
    canonically by colour. Vertices in different cells can never map to each
    other under an isomorphism, which collapses the permutation search.
    """
    n = len(A)
    colour = [(sum(A[i]), sum(A[j][i] for j in range(n))) for i in range(n)]
    while True:
        sig = []
        for i in range(n):
            out = sorted(colour[j] for j in range(n) if i != j and A[i][j])
            inn = sorted(colour[j] for j in range(n) if i != j and A[j][i])
            sig.append((colour[i], tuple(out), tuple(inn)))
        uniq = sorted(set(sig))
        newc = [uniq.index(s) for s in sig]
        if newc == colour or all(isinstance(c, int) for c in colour) and newc == colour:
            break
        if [sorted(set(newc)).index(c) for c in newc] == \
           [sorted(set(colour)).index(c) for c in colour] and len(set(newc)) == len(set(colour)):
            colour = newc
            break
        colour = newc
    cells = {}
    for i, c in enumerate(colour):
        cells.setdefault(c, []).append(i)
    return [cells[c] for c in sorted(cells)]


def canonical_form(A):
    """Lexicographically least adjacency key over all isomorphic relabelings.

    Uses colour refinement to restrict the search to permutations that respect
    the refined partition. Exact -- never heuristic. Returns (key, n).
    """
    n = len(A)
    cells = _refine(A)
    best = None
    # Only permutations that keep each refined cell setwise fixed can be
    # canonical candidates, because refinement is isomorphism-invariant.
    def gen(idx, prefix):
        nonlocal best
        if idx == len(cells):
            cand = to_key(permute(A, prefix))
            if best is None or cand < best:
                best = cand
            return
        for perm in permutations(cells[idx]):
            gen(idx + 1, prefix + list(perm))

    gen(0, [])
    return best, n


def from_key(key, n):
    """Rebuild an adjacency matrix from a canonical key."""
    A = [[0] * n for _ in range(n)]
    it = iter(key)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i][j] = next(it)
    return A


# ------------------------------------------------------------------ enumeration


def enumerate_digraphs(n, progress=None):
    """All digraphs on n nodes up to isomorphism, as canonical keys.

    Brute force over 2^(n(n-1)) labelled graphs for n <= 4, then vertex
    augmentation from the (n-1)-node representatives for larger n.
    """
    if n == 0:
        return {()}
    if n <= 4:
        seen = set()
        m = n * (n - 1)
        positions = [(i, j) for i in range(n) for j in range(n) if i != j]
        for mask in range(1 << m):
            A = [[0] * n for _ in range(n)]
            for b, (i, j) in enumerate(positions):
                if mask >> b & 1:
                    A[i][j] = 1
            seen.add(canonical_form(A)[0])
        return seen

    smaller = enumerate_digraphs(n - 1)
    seen = set()
    slots = 2 * (n - 1)          # edges between the new vertex and the old ones
    for count, key in enumerate(smaller):
        base = from_key(key, n - 1)
        for mask in range(1 << slots):
            A = [row[:] + [0] for row in base] + [[0] * n]
            for k in range(n - 1):
                if mask >> (2 * k) & 1:
                    A[k][n - 1] = 1
                if mask >> (2 * k + 1) & 1:
                    A[n - 1][k] = 1
            seen.add(canonical_form(A)[0])
        if progress and count % progress == 0:
            print(f"  augmented {count}/{len(smaller)} -> {len(seen)} so far",
                  flush=True)
    return seen


# ------------------------------------------------------------------ classes


def is_symmetric(A):
    n = len(A)
    return all(A[i][j] == A[j][i] for i in range(n) for j in range(n))


def is_oriented(A):
    """No bidirectional edge anywhere."""
    n = len(A)
    return all(not (A[i][j] and A[j][i])
               for i in range(n) for j in range(n) if i != j)


def is_mixed(A):
    """Neither symmetric nor oriented -- the only class still open."""
    return not is_symmetric(A) and not is_oriented(A)


def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)


# ------------------------------------------------------------------ domination
# Curto, arXiv:2510.05098 (Oct 2025), Definition 1.1, VERBATIM:
#   "Let j, k in [n] be vertices of G. We say that k graphically dominates j in G,
#    and write k > j, if the following two conditions hold:
#    (i)  For each vertex i in [n] \ {j, k}, if i -> j then i -> k.
#    (ii) j -> k and k -/-> j.
#    If there exists a k that graphically dominates j, we say that j is a
#    dominated node of G. If G has no dominated nodes, we say that it is
#    domination free."
# Corollary 3.6 there (any eps_i, delta_i > 0, no nondegeneracy needed): if k
# input-dominates j then no fixed point has x_j > 0; Lemma 3.8: graphical
# domination implies input domination for every CTLN. Hence a graph with a
# dominated node has NO full-support fixed point at any legal parameters.
# Theorem 1 there: FP(G) = FP(G minus j) for a dominated node j.
#
# The earlier placeholder `dominates()` (written from a secondhand summary)
# omitted condition (ii) and added an out-neighbour condition that is not in the
# paper; it was never used in a sweep and has been removed. See notes/lemmas.md
# for the proof that stable supports are also preserved (needed for Stage 4).


def graphically_dominates(A, k, j):
    """True iff k graphically dominates j in G (Definition 1.1 above)."""
    n = len(A)
    if j == k or A[j][k] != 1 or A[k][j] != 0:          # (ii)  j -> k and k -/-> j
        return False
    for i in range(n):                                   # (i)
        if i != j and i != k and A[i][j] == 1 and A[i][k] == 0:
            return False
    return True


def dominated_nodes(A):
    n = len(A)
    return [j for j in range(n) if any(graphically_dominates(A, k, j) for k in range(n) if k != j)]


def is_domination_free(A):
    return not dominated_nodes(A)
