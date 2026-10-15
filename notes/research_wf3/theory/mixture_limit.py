"""HEURISTIC: multi-sink-class limit.  For a core D (base graph G) and beta, light sets B (sum_{i in B} u_i < 1,
u = (D^T + I/beta)^{-1} 1 > 0).  Many sink classes with in-set complements B_j and sizes k_j -> inf (ratios theta_j)
have limiting moving spectrum spec(D - bbar 1^T), bbar = sum theta_j 1_{B_j} (Lemma X1 with B = sum k_j b_j).
Search: maximise beta such that  Pos(G, beta)  and  min Re spec(D - bbar 1^T) > -1/beta  for some mixture bbar of light sets."""
import numpy as np, json, sys, itertools
import tlib
def light_sets(A, beta):
    m = len(A); D = 1 - np.eye(m) - np.array(A)
    u = np.linalg.solve(D.T + np.eye(m) / beta, np.ones(m))
    if u.min() <= 1e-12 or u.sum() <= 1: return None, None
    sets = [np.array([(bits >> i) & 1 for i in range(m)], float) for bits in range(1, 1 << m)]
    return D, [b for b in sets if b @ u < 1 - 1e-9]
def best_mix(D, sets, rng, tries=300):
    m = len(D); best = (-9, None)
    S = np.array(sets)
    for t in range(tries):
        if t < len(sets): th = np.zeros(len(sets)); th[t] = 1
        else: th = rng.dirichlet(np.ones(len(sets)) * 0.2)
        bb = th @ S
        v = np.linalg.eigvals(D - np.outer(bb, np.ones(m))).real.min()
        if v > best[0]: best = (v, th)
    # local refinement
    th = best[1]; v0 = best[0]
    for it in range(300):
        th2 = np.abs(th + rng.normal(0, 0.05, len(th))); th2 /= th2.sum()
        v = np.linalg.eigvals(D - np.outer(th2 @ S, np.ones(m))).real.min()
        if v > v0: th, v0 = th2, v
    return v0, th
if __name__ == '__main__':
    d6s = json.load(open(sys.argv[1]))
    rng = np.random.default_rng(0)
    out = []
    for s in d6s:
        A = tlib.from_d6(s); bestb = None
        for beta in np.linspace(2.05, 6.0, 80):
            D, sets = light_sets(A, beta)
            if not sets: continue
            v, th = best_mix(D, sets, rng, tries=120)
            if v > -1 / beta: bestb = (beta, v)
        out.append((bestb[0] if bestb else 0, s, bestb))
        print(s, 'max beta with a stabilising light mixture (limit):', bestb, flush=True)
    out.sort(key=lambda x: -x[0])
    json.dump(out, open('mixture_limit_out.json', 'w'))
