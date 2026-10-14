"""HEURISTIC ONLY. Graph-level pure-stability form over all n=6 digraphs (coarse grid, then refine)."""
import numpy as np, sys
sys.path.insert(0, '/home/user/CLTN')
from enum_cache import key_from_int
from graphs import from_key
from scan_graph_level import nu_inf
cs = 1 + np.logspace(-3, 3.5, 120)
nuP2 = cs - np.sqrt(cs * cs - cs)
keys = [int(l) for l in open('/home/user/CLTN/certificates/digraphs_n6.keys')]
lo, hi = int(sys.argv[1]), int(sys.argv[2])
keys = keys[lo:hi]
J = np.ones((6, 6))
nst = 0; nviol = 0; close = []
for b0 in range(0, len(keys), 20000):
    As = np.array([from_key(key_from_int(v, 6), 6) for v in keys[b0:b0+20000]], dtype=float)
    best = np.full(len(As), -np.inf)
    for k, c in enumerate(cs):
        nu = np.linalg.eigvals(c * J - As).real.min(axis=1)
        best = np.maximum(best, nu - nuP2[k])
    ni = nu_inf(As)
    st = best > 1e-9
    nst += st.sum()
    v = st & (ni < 0.5 - 1e-9); nviol += v.sum()
    # near misses: graphs with nu_inf < 1/2 and margin > -0.02
    for b in np.nonzero((ni < 0.5 - 1e-9) & (best > -0.02))[0]:
        close.append((keys[b0+b], float(best[b]), float(ni[b])))
    print(lo + b0 + len(As), nst, nviol, len(close), flush=True)
close.sort(key=lambda r: -r[1])
print('TOP near-misses (key, margin, nu_inf):')
for r in close[:30]: print(r)
