"""
HEURISTIC ONLY (float). Large-t first-order correction at the critical limit eigenvalue.
X_c = cJ - Q, Q = A^T + gamma I, gamma = 1 - 1/beta.  Eigenvalues (moving) solve S(k) = 1/c,
S(k) = 1^T (Q + k I)^{-1} 1.  Near a zero k0 of S: k ~ k0 + 1/(c S'(k0)).
At beta_c (alpha = -gamma) the critical k0 = i y. M near the limit boundary needs Re S'(k0) <= 0
(finite c less stable than c = inf).  Report sign of Re S'(k0) and positivity at beta_c.
"""
import numpy as np, sys

def compress(A):
    m = len(A)
    V = np.zeros((m, m - 1))
    for j in range(m - 1): V[j, j] = 1; V[m - 1, j] = -1
    U, _ = np.linalg.qr(V)
    return U

def analyse(A):
    m = len(A)
    U = compress(A)
    Ap = U.T @ A @ U
    ev = np.linalg.eigvals(Ap)
    alpha = ev.real.max()
    if not alpha < -0.5:
        return None
    gamma = -alpha
    beta_c = 1 / (1 - gamma)
    Q = A.T + gamma * np.eye(m)
    # critical eigenvalue of X_inf: k0 = -(gamma + lambda) for lambda in spec(A|perp) with Re = alpha
    lam = ev[np.argmax(ev.real)]
    k0 = -(gamma + lam)
    R = np.linalg.inv(Q + k0 * np.eye(m))
    one = np.ones(m)
    S = one @ R @ one
    Sp = -(one @ R @ R @ one)
    # positivity at beta_c: w = Q^{-1} 1 > 0 and sign
    w = np.linalg.solve(Q, one)
    C = np.eye(m) + beta_c * (np.ones((m, m)) - np.eye(m) - A.T)
    a = np.linalg.det(C) * np.linalg.solve(C, one)
    return dict(alpha=alpha, beta_c=beta_c, k0=k0, S=S, ReSp=Sp.real, Sp=Sp, pos=bool((a > 1e-9).all()), wpos=bool((w > 0).all()))

if __name__ == '__main__':
    rng = np.random.default_rng(0)
    st = {}
    ex = []
    for trial in range(200000):
        m = int(rng.integers(3, 10))
        p = rng.uniform(0.5, 1.0)
        A = (rng.random((m, m)) < p).astype(float); np.fill_diagonal(A, 0)
        r = analyse(A)
        if r is None: continue
        if abs(r['k0'].imag) < 1e-9:  # real critical eigenvalue: then k0 = 0 crossing; skip
            key = ('real', r['pos'])
        else:
            key = ('cplx', r['pos'], r['ReSp'] <= 1e-12)
        st[key] = st.get(key, 0) + 1
        if key == ('cplx', True, False) and len(ex) < 5: ex.append((A.astype(int).tolist(), r))
    for k in sorted(st, key=str): print(k, st[k])
    for e in ex: print(e)
