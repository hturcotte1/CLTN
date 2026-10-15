"""HEURISTIC core design for the many-sink-class limit.  For a core digraph (m nodes) and beta:
light single nodes  = {i : u_i < 1},  u = (D^T + I/beta)^{-1} 1 > 0 (and sum u > 1);
limit moving polynomial p_theta = det(zI - D + bbar 1^T), bbar = sum theta_i e_i (theta in the simplex on light nodes).
score(core) = largest beta on a grid for which positivity holds and some theta (pairs of light nodes, theta grid) gives
min Re p_theta > -1/beta.  Hill-climb over single-arc flips.  Floats steer only; any find must be certified exactly."""
import numpy as np, sys, json, time
import tlib
BGRID = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
TH = np.linspace(0, 1, 21)
def best_minre(D, light):
    m = len(D); best = -9; arg = None
    for a in range(len(light)):
        for b in range(a, len(light)):
            for th in (TH if a != b else [1.0]):
                bb = np.zeros(m); bb[light[a]] += th; bb[light[b]] += 1 - th
                v = np.linalg.eigvals(D - np.outer(bb, np.ones(m))).real.min()
                if v > best: best, arg = v, (light[a], light[b], th)
    return best, arg
def score(A):
    A = np.asarray(A, float); m = len(A); D = 1 - np.eye(m) - A
    if D.sum() == 0: return 0.0, None
    sc, info = 0.0, None
    for beta in BGRID:
        try: u = np.linalg.solve(D.T + np.eye(m) / beta, np.ones(m))
        except np.linalg.LinAlgError: break
        if u.min() <= 1e-9 or u.sum() <= 1: break
        light = [i for i in range(m) if u[i] < 1 - 1e-9]
        if not light: break
        v, arg = best_minre(D, light)
        if v > -1 / beta: sc, info = beta, (v, arg)
        else: break
    return sc, info
def climb(A0, iters, rng):
    A = np.array(A0, int); s, inf = score(A); best = (s, A.copy(), inf)
    for it in range(iters):
        m = len(A); i, j = rng.integers(0, m, 2)
        if i == j: continue
        B = A.copy(); B[i, j] ^= 1
        sb, ib = score(B)
        if sb >= s:
            A, s, inf = B, sb, ib
            if s > best[0]: best = (s, A.copy(), inf); print('   improved', s, tlib.to_d6(A.tolist()), inf, flush=True)
    return best
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    seeds = sys.argv[1].split(',')
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    for s in seeds:
        A = tlib.from_d6(s)
        print('seed', s, 'score', score(A), flush=True)
        b = climb(A, iters, rng)
        print('BEST from', s, ':', b[0], tlib.to_d6(b[1].tolist()), b[2], flush=True)
