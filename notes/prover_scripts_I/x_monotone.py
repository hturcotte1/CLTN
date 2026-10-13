"""EXPLORATION ONLY (floats; not a decision path): is phi(tau) := alpha(A - tau J) nonincreasing for tau >= 1?
By Lemma I.3(g), (3) at (beta,t) <=> phi(1 + t/beta) < -gamma, so monotone phi would make the stable t-set at fixed
beta an up-interval. Scan all n=4,5 digraphs (and an n=6 sample) on a log grid of tau."""
import sys, random
sys.path.insert(0, '/home/user/CLTN')
import numpy as np, graphs, enum_cache
taus = np.concatenate([np.linspace(1, 3, 81), np.logspace(np.log10(3.05), 4, 120)])
def scan(n, keys):
    bad = []; J = np.ones((n, n))
    for v in keys:
        A = np.array(graphs.from_key(enum_cache.key_from_int(v, n), n), float)
        phi = np.array([np.linalg.eigvals(A - t*J).real.max() for t in taus])
        inc = np.diff(phi)
        if inc.max() > 1e-7:
            k = int(inc.argmax()); bad.append((inc.max(), taus[k], taus[k+1], phi[k], phi[k+1], A.astype(int).tolist()))
    return bad
for n in [3, 4, 5]:
    keys = enum_cache.load_keys(n, f'/home/user/CLTN/certificates/digraphs_n{n}.keys') if n >= 5 else None
    if keys is None:
        keys = sorted(enum_cache.key_to_int(k) for k in graphs.enumerate_digraphs(n))
    bad = scan(n, keys)
    print(f"n={n}: graphs={len(keys)}  with phi increasing somewhere on grid: {len(bad)}")
    for b in sorted(bad, reverse=True)[:3]: print("   max rise %.4g between tau=%.4g and %.4g (phi %.6f -> %.6f)  A=%s" % b)
