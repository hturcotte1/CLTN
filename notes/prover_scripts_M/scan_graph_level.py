"""
HEURISTIC ONLY. Graph-level form (pure stability): nu_G(c) := min Re spec(cJ - A).
legal&stable somewhere at c  <=>  nu_G(c) > nuP2(c) := c - sqrt(c^2 - c)   (gamma in (nuP2(c), nu_G(c)))
M_graph(pure): if nu_G(c) > nuP2(c) for some c>1 then nu_G(inf) >= 1/2.
Records max over c of nu_G(c) - nuP2(c) and nu_G(inf).
"""
import numpy as np, sys
sys.path.insert(0, '/home/user/CLTN')
from run_scan_small import labeled, iso
cs = 1 + np.logspace(-4, 4, 500)
nuP2 = cs - np.sqrt(cs * cs - cs)
def nu_inf(As):
    B, m, _ = As.shape
    V = np.zeros((m, m - 1))
    for j in range(m - 1): V[j, j] = 1; V[m - 1, j] = -1
    U, _ = np.linalg.qr(V)
    P = np.einsum('ji,bjk,kl->bil', U, -As, U)
    return np.linalg.eigvals(P).real.min(axis=1)
def run(As):
    B, m, _ = As.shape
    J = np.ones((m, m))
    best = np.full(B, -np.inf); argc = np.zeros(B)
    for k, c in enumerate(cs):
        nu = np.linalg.eigvals(c * J - As).real.min(axis=1)
        d = nu - nuP2[k]
        upd = d > best; best[upd] = d[upd]; argc[upd] = c
    return best, argc, nu_inf(As)
if __name__ == '__main__':
    for m in [2, 3, 4, 5]:
        As = np.array(list(labeled(m)) if m <= 4 else list(iso(m)))
        best, argc, ni = run(As)
        st = best > 1e-9
        viol = st & (ni < 0.5 - 1e-9)
        print(m, len(As), 'stable-legal-somewhere', st.sum(), 'with nu_inf<1/2:', viol.sum())
        for b in np.nonzero(viol)[0][:5]:
            print('   ', As[b].astype(int).tolist(), 'margin', best[b], 'at c', argc[b], 'nu_inf', ni[b])
