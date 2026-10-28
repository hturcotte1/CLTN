"""dynamics_float: float ODE sanity lens for CTLN full-support fixed points.
NOT a decision engine. Own implementations from the definitions in CLAUDE.md s.2.
"""
import json, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp

# ---------- digraph6 ----------
def parse_d6(s):
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        bits.extend(((v >> k) & 1) for k in range(5, -1, -1))
    assert len(bits) >= n * n, (s, len(bits), n)
    A = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        assert A[i][i] == 0, s
    return n, A

def is_complete(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

# ---------- CTLN ----------
def weight_matrix(A, eps, delta):
    """W[i][j] = 0 if i==j; -1+eps if A[j][i]==1 (j sends to i); -1-delta otherwise."""
    n = len(A)
    W = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            W[i][j] = F(-1) + eps if A[j][i] == 1 else F(-1) - delta
    return W

def solve_exact(M, b):
    """Gaussian elimination over Fractions. Returns x or None if singular."""
    n = len(M)
    aug = [list(M[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = None
        for r in range(c, n):
            if aug[r][c] != 0:
                p = r; break
        if p is None:
            return None
        aug[c], aug[p] = aug[p], aug[c]
        piv = aug[c][c]
        aug[c] = [v / piv for v in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [a - f * bb for a, bb in zip(aug[r], aug[c])]
    return [aug[i][n] for i in range(n)]

def fixed_point_exact(W, theta=F(1)):
    n = len(W)
    M = [[(F(1) if i == j else F(0)) - W[i][j] for j in range(n)] for i in range(n)]
    return solve_exact(M, [theta] * n)

def jac_eigs(W):
    """float eigenvalues of -I + W (Jacobian in the all-active linear regime)."""
    Wf = np.array([[float(v) for v in row] for row in W])
    n = len(W)
    return np.linalg.eigvals(-np.eye(n) + Wf), Wf

# ---------- simulation ----------
def simulate(Wf, x0, T, theta=1.0, method='LSODA', rtol=1e-10, atol=1e-13, n_eval=201):
    n = Wf.shape[0]
    def rhs(t, x):
        return -x + np.maximum(Wf @ x + theta, 0.0)
    def jac(t, x):
        act = (Wf @ x + theta) > 0
        return -np.eye(n) + Wf * act[:, None]
    t_eval = np.linspace(0, T, n_eval)
    kw = dict(method=method, rtol=rtol, atol=atol, t_eval=t_eval)
    if method in ('LSODA', 'BDF', 'Radau'):
        kw['jac'] = jac
    sol = solve_ivp(rhs, (0.0, T), x0, **kw)
    return sol

def initial_conditions(xf, n, seed):
    """10 ICs: 5 perturbed (x*(1+0.3u)), 5 uniform [0,2]^n. xf may contain negatives -> clip at 0."""
    rng = np.random.default_rng(seed)
    ics = []
    for k in range(5):
        u = rng.uniform(-1, 1, n)
        ics.append(('perturb', np.maximum(xf * (1 + 0.3 * u), 0.0)))
    for k in range(5):
        ics.append(('uniform', rng.uniform(0, 2, n)))
    return ics

def analyze_trajectory(sol, xf):
    """relative max error vs exact xf along t_eval; convergence flags."""
    X = sol.y.T  # (m, n)
    scale = np.max(np.abs(xf)) if np.max(np.abs(xf)) > 0 else 1.0
    err = np.max(np.abs(X - xf[None, :]), axis=1) / scale
    m = len(err)
    half = err[m // 2:]
    strict_noninc = bool(np.all(half[1:] <= half[:-1] * (1 + 1e-6) + 1e-15))
    viol = int(np.sum(half[1:] > half[:-1] * (1 + 1e-6) + 1e-15))
    xe = X[-1]
    return dict(
        ok=bool(sol.success), nfev=int(sol.nfev), nsteps=int(len(sol.t)),
        err_end=float(err[-1]), err_half=float(half[0]), err_max_lasthalf=float(np.max(half)),
        converged=bool(err[-1] < 1e-7 and strict_noninc),
        err_end_lt_1e7=bool(err[-1] < 1e-7), noninc_lasthalf=strict_noninc, noninc_violations=viol,
        all_pos_end=bool(np.all(xe > 0)), min_x_end=float(np.min(xe)),
        x_end=[float(v) for v in xe],
    )

def frac(s):
    return F(s)

# ---------- end-state classification ----------
def perturbed_ics(xf, n, seed, scales=(0.3, 0.03, 0.003), k=5):
    rng = np.random.default_rng(seed)
    ics = []
    for s in scales:
        for _ in range(k):
            u = rng.uniform(-1, 1, n)
            ics.append((f'perturb{s}', np.maximum(xf * (1 + s * u), 0.0)))
    for _ in range(k):
        ics.append(('uniform', rng.uniform(0, 2, n)))
    return ics

def classify_end(W, Wf, xe, theta=1.0):
    """support of the end state, residual of dx/dt, and whether it matches the exact fixed point on that support
    (conditions 1 and 2 checked exactly on the support; stability not decided)."""
    n = len(W)
    inp = Wf @ xe + theta
    S = [i for i in range(n) if inp[i] > 0]
    dx = -xe + np.maximum(inp, 0.0)
    res = float(np.max(np.abs(dx)))
    info = dict(support=S, size=len(S), residual=res)
    if not S:
        info['exact_fp_match'] = False; return info
    M = [[(F(1) if i == j else F(0)) - W[i][j] for j in S] for i in S]
    xs = solve_exact(M, [F(1)] * len(S))
    if xs is None:
        info['exact_fp_match'] = False; info['note'] = 'singular on support'; return info
    xfull = [F(0)] * n
    for a, i in enumerate(S): xfull[i] = xs[a]
    pos = all(v > 0 for v in xs)
    off_ok = all(sum(W[k][j] * xfull[j] for j in S) + 1 < 0 for k in range(n) if k not in S)
    xff = np.array([float(v) for v in xfull])
    err = float(np.max(np.abs(xe - xff)) / max(np.max(np.abs(xff)), 1e-300))
    info.update(exact_fp_positive=pos, exact_fp_offsupport_ok=off_ok, rel_err_to_support_fp=err,
                exact_fp_match=bool(pos and off_ok and err < 1e-5))
    return info

def linear_solution(Wf, xf, x0, T, theta=1.0):
    """analytic solution of the all-active linear system dx/dt = -x + Wx + theta at time T: xf + expm(J T)(x0-xf)."""
    from scipy.linalg import expm
    J = -np.eye(Wf.shape[0]) + Wf
    return xf + expm(J * T) @ (x0 - xf)
