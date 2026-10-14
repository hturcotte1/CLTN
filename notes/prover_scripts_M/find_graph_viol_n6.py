import numpy as np, sys
sys.path.insert(0, '/home/user/CLTN')
from enum_cache import key_from_int
from graphs import from_key, is_domination_free
from scan_graph_level import nu_inf
from to_d6 import to_d6
cs = 1 + np.logspace(-3, 3.5, 120)
nuP2 = cs - np.sqrt(cs * cs - cs)
keys = [int(l) for l in open('/home/user/CLTN/certificates/digraphs_n6.keys')]
lo, hi = int(sys.argv[1]), int(sys.argv[2])
J = np.ones((6, 6))
for b0 in range(lo, hi, 20000):
    ks = keys[b0:min(hi, b0 + 20000)]
    As = np.array([from_key(key_from_int(v, 6), 6) for v in ks], dtype=float)
    best = np.full(len(As), -np.inf); argc = np.zeros(len(As))
    for k, c in enumerate(cs):
        nu = np.linalg.eigvals(c * J - As).real.min(axis=1)
        d = nu - nuP2[k]; u = d > best; best[u] = d[u]; argc[u] = c
    ni = nu_inf(As)
    for b in np.nonzero((best > 1e-9) & (ni < 0.5 - 1e-9))[0]:
        A = As[b].astype(int).tolist()
        print(b0 + b, ks[b], to_d6(A), 'margin', best[b], 'at c', argc[b], 'nu_inf', ni[b], 'domfree', is_domination_free(A), A, flush=True)
