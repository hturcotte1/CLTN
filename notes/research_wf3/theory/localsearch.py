"""HEURISTIC local search (floats steer only).  Maximise q_lim(G) = sup{beta - 1 : 2 < beta < beta_c, w(beta) > 0}
over digraphs on n nodes by single-arc flips (best improvement + random kicks).  Every reported graph must be
certified exactly afterwards (certify_q.py)."""
import numpy as np, json, sys, time
import heur_extend as he, tlib
BETAS = np.concatenate([np.linspace(2.0005, 2.6, 200), np.linspace(2.603, 6.0, 500)])
def scores(Bs):
    Bs = np.asarray(Bs, dtype=np.int8)
    n = Bs.shape[1]; U = he.ortho_perp(n)
    C = np.einsum('ia,kij,jb->kab', U, Bs.astype(float), U)
    a = np.linalg.eigvals(C).real.max(axis=1)
    out = np.empty(len(Bs))
    BT = np.transpose(Bs, (0, 2, 1)).astype(float)
    for k in range(len(Bs)):
        if a[k] >= -0.5:
            out[k] = -3.0 - a[k]; continue
        bc = 1 / (1 + a[k])
        bs = BETAS[BETAS < bc]
        if len(bs) == 0: out[k] = -2.5; continue
        g = 1 - 1 / bs
        M = BT[k][None] + g[:, None, None] * np.eye(n)[None]
        try:
            w = np.linalg.solve(M, np.ones((len(bs), n, 1)))[..., 0]
        except np.linalg.LinAlgError:
            out[k] = -2.5; continue
        pm = w.min(axis=1) / np.abs(w).max(axis=1)
        ok = pm > 1e-9
        out[k] = (bs[ok].max() - 1) if ok.any() else (-2.0 + pm.max())
    return out
def neighbours(A):
    n = len(A); res = []
    for i in range(n):
        for j in range(n):
            if i != j:
                B = A.copy(); B[i, j] ^= 1; res.append(B)
    return np.array(res)
def search(A0, iters, rng, kick=3):
    A = np.array(A0, dtype=np.int8); s = scores([A])[0]
    best = (s, A.copy()); stall = 0
    for it in range(iters):
        nb = neighbours(A); sc = scores(nb)
        k = int(np.argmax(sc))
        if sc[k] > s + 1e-9:
            A, s = nb[k], sc[k]; stall = 0
            if s > best[0]: best = (s, A.copy())
        else:
            stall += 1
            # random kick from the best
            A = best[1].copy()
            for _ in range(kick):
                i, j = rng.integers(0, len(A), 2)
                if i != j: A[i, j] ^= 1
            s = scores([A])[0]
    return best
if __name__ == '__main__':
    seeds = json.load(open(sys.argv[1])); iters = int(sys.argv[2]); tag = sys.argv[3]
    rng = np.random.default_rng(1)
    res = []
    t0 = time.time()
    for d6 in seeds:
        b = search(tlib.from_d6(d6), iters, rng)
        res.append((float(b[0]), tlib.to_d6(b[1].tolist()), d6))
        print('%.4f' % b[0], tlib.to_d6(b[1].tolist()), 'from', d6, '%.0fs' % (time.time() - t0), flush=True)
    res.sort(key=lambda x: -x[0])
    json.dump(res, open('ls_%s.json' % tag, 'w'), indent=0)
