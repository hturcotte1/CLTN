"""
forest_exact.py -- EXACT tools for the in-tree / arms-forest construction (no floating point in any decision).

Conventions: CLAUDE.md section 2.  A[i][j] = 1 iff i -> j.  D := J - I - A (D[i][j] = 1 iff i != j and i does NOT send to j).
beta = 1 + delta/eps, r := 1/beta, gamma = 1 - r.  (notes/lemma_I_spectral.md, notes/research_wf3/theory/notes_theory.md.)

Arms-forest core for an integer polynomial p(z) = z^d + z^(d-1) + c_2 z^(d-2) + ... + c_d  (c_2 >= 1, c_j >= 0):
  nodes c0, c1, c2 with D-arcs c0 -> c1 -> c2 -> c0  (c0 does not send to c1, etc.);
  for every 'arm' a: D-arcs c2 -> a, a -> z_a, z_a -> target(a)   (a, z_a are new nodes);
  arm lengths len(a) = 2 depth(a) + ell(exit), ell(c2, c1, c0) = (0, 1, 2);
  c_2 - 1 arms of length 2 (exit c2, depth 1); c_j arms of length j for j >= 3:
     odd j -> exit c1, depth (j-1)/2;  even j >= 4 -> exit c0, depth j/2 - 1;
  an arm of depth k >= 2 targets the first arm of depth k-1 with the same exit ('spine').
  Every node except c2 has exactly one D-out-neighbour, and all D-paths lead to c2 (an in-tree rooted at c2).

Lemma T (proved in the notes):  chi_{D - e_{c2} 1^T}(z) = z^(m-d) p(z).
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from fastregion import icharpoly

def build_core(c):
    """c = [1, 1, c2, ..., cd] integers.  Returns (D, names) with D an m x m 0/1 integer matrix (lists)."""
    c = [int(x) for x in c]
    assert c[0] == 1 and c[1] == 1 and len(c) >= 3 and c[2] >= 1 and all(x >= 0 for x in c)
    d = len(c) - 1
    names = ['c0', 'c1', 'c2']; arcs = [('c0', 'c1'), ('c1', 'c2'), ('c2', 'c0')]
    first = {}
    spec = [(2, 'c2', 1, c[2] - 1)]
    for j in range(3, d + 1):
        if j % 2 == 1: spec.append((j, 'c1', (j - 1) // 2, c[j]))
        else: spec.append((j, 'c0', j // 2 - 1, c[j]))
    spec.sort(key=lambda s: s[2])                         # by depth, so parents exist first
    for j, ex, dep, cnt in spec:
        for i in range(cnt):
            a = 'a%d_%d' % (j, i); z = 'z%d_%d' % (j, i)
            if dep == 1: tgt = ex
            else:
                if (ex, dep - 1) not in first:
                    raise ValueError('forest constraint: no arm of depth %d with exit %s' % (dep - 1, ex))
                tgt = first[(ex, dep - 1)]
            names += [a, z]; arcs += [('c2', a), (a, z), (z, tgt)]
            first.setdefault((ex, dep), a)
    idx = {v: i for i, v in enumerate(names)}; m = len(names)
    D = [[0] * m for _ in range(m)]
    for u, v in arcs: D[idx[u]][idx[v]] = 1
    return D, names

def A_of_D(D):
    m = len(D); return [[0 if i == j else 1 - D[i][j] for j in range(m)] for i in range(m)]

def L_matrix(D, names, root='c2'):
    """integer matrix L = D - e_root 1^T  (the sink-class limit matrix, Lemma X1 with b = e_root, c = 0)"""
    m = len(D); k = names.index(root)
    return [[D[i][j] - (1 if i == k else 0) for j in range(m)] for i in range(m)]

def check_lemma_T(c):
    """exact: integer charpoly of L equals z^(m-d) p."""
    D, names = build_core(c); m = len(D); d = len(c) - 1
    chi = icharpoly(L_matrix(D, names))
    return chi == [int(x) for x in c] + [0] * (m - d), m

def roots_right_of(p_hi_lo, s):
    """exact three-valued test: every root of p has Re > -s  (s rational > 0).
    Q(w) := sum_j c_j (-1)^j (w + s)^(d-j)  has roots w = -z - s; Hurwitz(Q) <=> all Re z > -s."""
    from_tlib = __import__('tlib')
    d = len(p_hi_lo) - 1
    ptil = [F(cj) * (-1) ** j for j, cj in enumerate(p_hi_lo)]     # (-1)^d p(-y), hi->lo in y (monic)
    Q = from_tlib.taylor_shift(ptil, F(s))
    return ctln.hurwitz_stable([F(x) for x in Q])

def positivity(D, names, beta, root='c2'):
    """exact u = (D^T + I/beta)^{-1} 1 and the admissibility facts:
    u > 0, sum u > 1 (=> w = (A^T + gamma I)^{-1} 1 = u/(sum u - 1) > 0, Lemma D2) and u_root < 1 (sink class
    received from all but root is activated: sigma_S > 1, Lemma X2)."""
    m = len(D); r = 1 / F(beta)
    M = [[F(D[j][i]) + (r if i == j else 0) for j in range(m)] for i in range(m)]
    u = ctln.solve_exact(M, [F(1)] * m)
    if u is None: return None
    k = names.index(root)
    return {'u': u, 'u_min': min(u), 'sum': sum(u), 'u_root': u[k],
            'ok': min(u) > 0 and sum(u) > 1 and u[k] < 1}

def closed_form_lambda(c, beta):
    """closed form of lambda = 1 - eps_{c2}/r  (notes, Lemma P): lambda = r^3 (1 - r + r^2) / Den,
    Den = F0 - r F1 + r^2 F2 + r^2 + r^5."""
    r = 1 / F(beta); d = len(c) - 1
    F0 = sum(F(c[j]) * r ** (4 - j) for j in range(4, d + 1, 2))
    F1 = sum(F(c[j]) * r ** (3 - j) for j in range(3, d + 1, 2))
    F2 = F(c[2] - 1)
    Den = F0 - r * F1 + r * r * F2 + r ** 2 + r ** 5
    return r ** 3 * (1 - r + r * r) / Den, (F0, F1, F2, Den)
