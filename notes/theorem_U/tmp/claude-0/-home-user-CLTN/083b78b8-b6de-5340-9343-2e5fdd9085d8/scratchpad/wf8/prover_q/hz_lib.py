"""H/Z designs (exact construction + float exploration helpers).
H nodes 0..h-1 with a parent map par[v] in H (functional graph, D-arc par[v] -> v); T = h x h nonneg integer matrix:
row t of T is realised by max(T[t]) Z leaves with parent t, the k-th leaf having D-arcs to {t' : T[t][t'] >= k}.
D = non-edge matrix (D[i][j] = 1 iff i does not send to j).  A = J - I - D."""
import numpy as np
def build(par, T):
    h = len(par); T = np.asarray(T, int)
    leaves = []
    for t in range(h):
        for k in range(1, int(T[t].max(initial=0)) + 1):
            leaves.append((t, [tp for tp in range(h) if T[t][tp] >= k]))
    m = h + len(leaves); D = np.zeros((m, m), int)
    for v in range(h):
        if par[v] is not None: D[par[v], v] = 1
    for i, (t, S) in enumerate(leaves):
        z = h + i; D[t, z] = 1
        for tp in S: D[z, tp] = 1
    assert all(D[i, i] == 0 for i in range(m))
    return D
def PH(par):
    h = len(par); P = np.zeros((h, h))
    for v in range(h): P[par[v], v] = 1
    return P
def psi_roots(par, T, th):
    """roots of 1_H^T adj(z^2 I - z P_H - T) theta_H (float, exploration)"""
    import numpy.polynomial.polynomial as npp
    h = len(par); P = PH(par); T = np.asarray(T, float)
    # evaluate polynomial via interpolation on 2h-1 points
    deg = 2 * h - 2
    xs = np.exp(2j * np.pi * np.arange(deg + 1) / (deg + 1)) * 1.3
    vals = []
    for x in xs:
        M = x * x * np.eye(h) - x * P - T
        vals.append(np.linalg.det(M) * (np.ones(h) @ np.linalg.solve(M, th)))
    c = np.linalg.solve(np.vander(xs, deg + 1, increasing=True), np.array(vals))
    return np.roots(c[::-1].real)
