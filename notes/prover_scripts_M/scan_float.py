"""
HEURISTIC ONLY (floating point). Search for violations of Lemma M (pure stability form).

Formulation (exact identity, see notes): with eps' = 1/beta, tau = t/beta,
  K/beta = eps' I + N + tau J,
so K positive stable  <=>  g(tau) := min Re spec(N + tau J) > -1/beta =: L.
Legal  t(beta-2) > 1  <=>  tau > L^2/(1+2L)  <=>  L > gP2(tau) := tau - sqrt(tau^2+tau).
M_s violation: tau1 < tau2 with g(tau1) > max(gP2(tau1), L), g(tau2) <= L, L in (-1/2,0).
  <=> exists tau1<tau2: max(gP2(tau1), g(tau2)) < min(g(tau1), 0)   (then pick L in between).
M_w violation: g(tau1) > gP2(tau1) and g(inf) < min(g(tau1),0).
"""
import sys, itertools, numpy as np

def g_curve(Ns, taus):
    # Ns: (B, m, m) float arrays. returns (B, len(taus))
    B, m, _ = Ns.shape
    J = np.ones((m, m))
    out = np.empty((B, len(taus)))
    for k, tau in enumerate(taus):
        ev = np.linalg.eigvals(Ns + tau * J)
        out[:, k] = ev.real.min(axis=1)
    return out

def g_inf(Ns):
    B, m, _ = Ns.shape
    V = np.zeros((m, m - 1))
    for j in range(m - 1):
        V[j, j] = 1; V[m - 1, j] = -1
    Q, _ = np.linalg.qr(V)
    P = np.einsum('ji,bjk,kl->bil', Q, Ns, Q)
    ev = np.linalg.eigvals(P)
    return ev.real.min(axis=1)

def gP2(tau):
    return tau - np.sqrt(tau * tau + tau)

def check(Ns, taus, tol=1e-9):
    g = g_curve(Ns, taus)
    gi = g_inf(Ns)
    gp = gP2(taus)
    # suffix min over tau2 > tau1 (including g_inf as tau -> inf)
    suf = np.minimum.accumulate(np.concatenate([g[:, ::-1], ], axis=1), axis=1)[:, ::-1]
    sufnext = np.concatenate([suf[:, 1:], np.full((len(Ns), 1), np.inf)], axis=1)
    sufnext = np.minimum(sufnext, gi[:, None])
    upper = np.minimum(g, 0.0)
    lower = np.maximum(gp[None, :], sufnext)
    viol_s = (lower < upper - tol) & (g > gp[None, :] + tol)
    lower_w = np.maximum(gp[None, :], gi[:, None])
    viol_w = (lower_w < upper - tol) & (g > gp[None, :] + tol)
    return g, gi, viol_s, viol_w

def N_from_A(A):
    m = A.shape[-1]
    return np.ones((m, m)) - np.eye(m) - np.swapaxes(A, -1, -2)

if __name__ == '__main__':
    pass
