"""
HEURISTIC ONLY (floating point). Stream digraph6 positivity survivors and flag any graph that is
(positive at beta) & (stable at some legal finite t) on a grid, then compare with the limit.
Output lines: d6  flag  details    flag in {FIN_AND_LIM, FIN_NOT_LIM}
"""
import sys, numpy as np
from scan_float import gP2, g_inf
from scan_pos import positivity_grid

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = np.zeros((n, n)); bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i, j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A

taus = np.logspace(-2.5, 3.5, 90)
epss = np.linspace(0.004, 0.496, 124)
gp = gP2(taus)

def process(batch, out):
    As = np.array([from_d6(s) for s in batch])
    m = As.shape[1]
    Ns = np.ones((m, m)) - np.eye(m) - np.swapaxes(As, 1, 2)
    J = np.ones((m, m))
    G = np.empty((len(As), len(taus)))
    for k, tau in enumerate(taus):
        G[:, k] = np.linalg.eigvals(Ns + tau * J).real.min(axis=1)
    cand = np.nonzero((G > gp[None, :] + 1e-9).any(axis=1))[0]
    if len(cand) == 0:
        return 0
    Nc = Ns[cand]; Gc = G[cand]
    pos = positivity_grid(Nc, epss)
    gi = g_inf(Nc)
    L = -epss
    leg = taus[:, None] > (epss**2 / (1 - 2 * epss))[None, :] * (1 + 1e-9)
    stab = (Gc[:, :, None] > L[None, None, :] + 1e-9) & leg[None] & pos[:, None, :]
    n = 0
    for b in range(len(cand)):
        if stab[b].any():
            ti, ei = np.nonzero(stab[b])
            limit_ok = gi[b] >= L[ei] - 1e-9
            flag = 'FIN_AND_LIM' if limit_ok.all() else 'FIN_NOT_LIM'
            k = np.argmin(limit_ok) if not limit_ok.all() else 0
            out.write(f"{batch[cand[b]]} {flag} ginf={gi[b]:.6f} npts={len(ti)} ex_tau={taus[ti[k]]:.5g} ex_eps'={epss[ei[k]]:.5g} g={Gc[b, ti[k]]:.6f}\n")
            n += 1
    out.flush()
    return n

if __name__ == '__main__':
    path, outp = sys.argv[1], sys.argv[2]
    start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    stop = int(sys.argv[4]) if len(sys.argv) > 4 else 10**12
    batch = []; tot = 0; hits = 0
    with open(path) as fh, open(outp, 'a') as out:
        for i, line in enumerate(fh):
            if i < start: continue
            if i >= stop: break
            batch.append(line.split()[0])
            if len(batch) == 4000:
                hits += process(batch, out); tot += len(batch); batch = []
                print(tot, hits, flush=True)
        if batch:
            hits += process(batch, out); tot += len(batch)
    print('done', tot, hits, flush=True)
