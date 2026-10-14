"""HEURISTIC. At beta = 2 (level L = -1/2): T(2) = {tau>0 : g(tau) > -1/2}, g = min Re spec(N + tau J).
M2: T(2) upward closed.   Mw(2): T(2) nonempty => g(inf) >= -1/2."""
import numpy as np, sys
sys.path.insert(0, '/home/user/CLTN')
from run_scan_small import labeled, iso
from scan_float import g_inf, N_from_A
taus = np.logspace(-4, 5, 600)
def run(As):
    Ns = N_from_A(As); B, m, _ = Ns.shape; J = np.ones((m, m))
    G = np.empty((B, len(taus)))
    for k, tau in enumerate(taus):
        G[:, k] = np.linalg.eigvals(Ns + tau * J).real.min(axis=1)
    gi = g_inf(Ns)
    st = G > -0.5 + 1e-9
    anyst = st.any(axis=1)
    # upward closure: after first stable index, all stable (and limit >= -1/2)
    first = np.where(anyst, st.argmax(axis=1), len(taus))
    viol_M2 = np.array([anyst[b] and (not st[b, first[b]:].all() or gi[b] < -0.5 - 1e-9) for b in range(B)])
    viol_Mw2 = anyst & (gi < -0.5 - 1e-9)
    return anyst, viol_M2, viol_Mw2, G, gi
if __name__ == '__main__':
    for m in [2, 3, 4, 5]:
        As = np.array(list(labeled(m)) if m <= 4 else list(iso(m)))
        anyst, v2, vw, G, gi = run(As)
        print(m, len(As), 'T(2) nonempty:', anyst.sum(), 'M2 viol:', v2.sum(), 'Mw(2) viol:', vw.sum())
        for b in np.nonzero(v2)[0][:5]:
            st = G[b] > -0.5
            print('    ', As[b].astype(int).tolist(), 'gi', gi[b], 'stable taus range', taus[st].min(), taus[st].max(), 'min g after first', G[b, st.argmax():].min())
