"""HEURISTIC. Is P2 extremal for the finite-c gain nu_G(c) - nu_G(inf)?"""
import numpy as np, sys
sys.path.insert(0, '/home/user/CLTN')
from run_scan_small import labeled, iso
from scan_graph_level import nu_inf
cs = np.array([1.0001, 1.01, 1.1, 1.5, 2, 3, 5, 10, 30, 100, 1000])
P2gain = cs - np.sqrt(cs * cs - cs) - 0.5
for m in [2, 3, 4, 5]:
    As = np.array(list(labeled(m)) if m <= 4 else list(iso(m)))
    ni = nu_inf(As)
    J = np.ones((m, m))
    row = []
    for k, c in enumerate(cs):
        nu = np.linalg.eigvals(c * J - As).real.min(axis=1)
        gain = nu - ni
        b = np.argmax(gain)
        row.append((c, round(gain[b], 5), round(P2gain[k], 5)))
        if gain[b] > P2gain[k] + 1e-9 and c >= 1.01:
            print('   m', m, 'c', c, 'gain', gain[b], '> P2', P2gain[k], As[b].astype(int).tolist(), 'nu', nu[b], 'nu_inf', ni[b])
    print(m, row)
