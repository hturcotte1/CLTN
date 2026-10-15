"""r3_instances.py -- referee: independent recomputation of the Theorem U instances (N = 3, 4, 5), the 65-node core
thresholds, and the Lemma S / positivity spot checks.  Writes r3_instances.json.  Exact only."""
import sys, json, time
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q')
sys.path.insert(0, '/home/user/CLTN')
from ref_lib import *
import ctln
from fastregion import idet

OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q/r3_instances.json'
which = sys.argv[1:] or ['U3', 'core65', 'S', 'U4', 'U5']
try: out = json.load(open(OUT))
except Exception: out = {}


def target_coeffs(N, W):
    s = sp.Rational(1, 2 * N)
    P = sp.Poly(sp.prod([(z + s) ** 2 + W ** 2 * k ** 2 for k in range(1, N + 1)]), z)
    return P.all_coeffs()


def nint(q):
    q = F(q.p, q.q)
    fl = q.numerator // q.denominator
    return fl + 1 if 2 * (q - fl) >= 1 else fl


def rounded(N, W):
    co = target_coeffs(N, W)
    assert co[0] == 1 and co[1] == 1
    return [1, 1] + [nint(a) for a in co[2:]]


def positivity_tree(c, beta):
    """referee's structured exact solve of (D^T + rI) u = 1 on the arms forest, without the closed form:
    unknown lam; arm f-values by spine recursion; cycle rows solved as a linear system in (u_c0,u_c1,u_c2,lam)."""
    r = 1 / F(beta); d = len(c) - 1
    groups = {('c2', 1): c[2] - 1}
    for j in range(3, d + 1):
        key = ('c1', (j - 1) // 2) if j % 2 else ('c0', j // 2 - 1)
        groups[key] = groups.get(key, 0) + c[j]
    # f of the spine arm at depth k for each exit: f = 1 + r^-2 * (sum of f over arms at depth k+1); non-spine arms are leaves
    fsp = {}; Fsum = {}; fmax = F(1)
    for ex in ['c0', 'c1', 'c2']:
        deps = sorted(dd for (e, dd) in groups if e == ex and groups[(e, dd)] > 0)
        if not deps: Fsum[ex] = F(0); continue
        maxd = max(deps)
        f_next_sum = F(0)
        for dd in range(maxd, 0, -1):
            cnt = groups.get((ex, dd), 0)
            assert cnt > 0, 'forest constraint'
            fs = 1 + f_next_sum / (r * r)
            fmax = max(fmax, fs)
            f_next_sum = fs + (cnt - 1)        # sum of f over arms at depth dd (spine + leaves)
        Fsum[ex] = f_next_sum                  # sum of f over depth-1 arms of this exit
    # unknowns u0,u1,u2,lam ; lam := 1 - (1-u2)/r
    u0, u1, u2, lam = sp.symbols('u0 u1 u2 lam')
    R = sp.Rational(r.numerator, r.denominator)
    Fs = {e: sp.Rational(v.numerator, v.denominator) for e, v in Fsum.items()}
    eqs = [sp.Eq(lam, 1 - (1 - u2) / R),
           sp.Eq(u2 + lam * Fs['c0'] / R + R * u0, 1),     # row c0: D-in = c2 + z's of exit c0
           sp.Eq(u0 + lam * Fs['c1'] / R + R * u1, 1),     # row c1: D-in = c0 + z's of exit c1
           sp.Eq(u1 + lam * Fs['c2'] / R + R * u2, 1)]     # row c2: D-in = c1 + z's of exit c2
    sol = sp.solve(eqs, [u0, u1, u2, lam], dict=True)[0]
    L_ = F(str(sol[lam])); U0 = F(str(sol[u0])); U1 = F(str(sol[u1])); U2 = F(str(sol[u2]))
    ok = L_ > 0 and L_ * fmax < 1 and U0 > 0 and U1 > 0 and 0 < U2 < 1
    return {'lam': str(L_), 'lam_fmax': str(L_ * fmax), 'u_c0': str(U0), 'u_c1': str(U1), 'u_c2': str(U2), 'ok': ok,
            'Fsum': {e: str(v) for e, v in Fsum.items()}}


def instance(c, beta):
    p, G, m = reduced_from_c(c)
    rr = sp.Rational(F(beta).denominator, F(beta).numerator)
    C1 = strip_right_of(p, rr)
    pos = positivity_tree(c, beta)
    t0 = time.time()
    k0, passes, det = k_threshold(p, G, rr)
    return {'c': c, 'beta': str(beta), 'm': m, 'C1_limit_poly_right_of_-1/beta': C1, 'positivity': pos, 'k0': k0,
            'fails_at_k0_minus_1': (not passes(k0 - 1)) if k0 > 1 else None, 'n0': m + k0, 'secs': round(time.time() - t0, 1),
            'minor_degrees_and_top_roots': det}


if 'U3' in which:
    res = []
    for N in [3]:
        c = rounded(N, 1)
        rec = instance(c, N); rec['N'] = N; rec['W'] = 1
        res.append(rec); print(json.dumps({k: v for k, v in rec.items() if k != 'minor_degrees_and_top_roots'}, default=str), flush=True)
    out['U3'] = res
    # direct exact positivity on the explicit 257-node core and one-point checks of Lemma T / psi identity
    c = rounded(3, 1); D, idx = arms_core(c); m = len(D)
    r = F(1, 3)
    M = [[F(D[j][i]) + (r if i == j else 0) for j in range(m)] for i in range(m)]
    t0 = time.time(); u = ctln.solve_exact(M, [F(1)] * m)
    U = sum(u)
    direct = {'m': m, 'u_min_pos': min(u) > 0, 'U_gt_1': U > 1, 'u_c2': str(u[idx['c2']]), 'u_c2_lt_1': u[idx['c2']] < 1,
              'lam_direct': str(1 - (1 - u[idx['c2']]) / r), 'secs': round(time.time() - t0, 1)}
    p, G, mm = reduced_from_c(c); d = len(c) - 1
    L = [[D[i][j] - (1 if i == 0 else 0) for j in range(m)] for i in range(m)]
    pts = {}
    for z0 in [2, -3]:
        chiL = idet([[(z0 if i == j else 0) - L[i][j] for j in range(m)] for i in range(m)])
        psi0 = idet([[(z0 if i == j else 0) - D[i][j] for j in range(m)] for i in range(m)]) - \
               idet([[(z0 if i == j else 0) - D[i][j] - 1 for j in range(m)] for i in range(m)])
        pts[z0] = {'chiL_ok': chiL == z0 ** (m - d) * p.eval(z0), 'psi_ok': psi0 == z0 ** (m - 1 - d) * G.eval(z0)}
    direct['pointwise_LemmaT_Q_on_257_core'] = pts
    out['U3_direct'] = direct
    print(direct, flush=True)
    json.dump(out, open(OUT, 'w'), indent=1, default=str)

if 'core65' in which:
    c65 = [1, 1, 11, 7, 10, 3, 1]
    res = []
    for beta in [F(4), F(5), F(11, 2), F(28, 5)]:
        rec = instance(c65, beta); res.append(rec)
        print(json.dumps({k: v for k, v in rec.items() if k != 'minor_degrees_and_top_roots'}, default=str), flush=True)
    # ceiling: at beta = 6 the limit polynomial must fail (sum of roots = -1, d = 6)
    p, G, m = reduced_from_c(c65)
    out['core65_beta6_limit_poly'] = strip_right_of(p, sp.Rational(1, 6))
    out['core65'] = res
    json.dump(out, open(OUT, 'w'), indent=1, default=str)

if 'S' in which:
    res = []
    for N in range(3, 8):
        W0 = 2 * N * 4 ** N * (N + 1) ** (2 * N - 2)
        c = rounded(N, W0)
        p = sp.Poly(sum(c[j] * z ** (2 * N - j) for j in range(2 * N + 1)), z)
        right = strip_right_of(p, sp.Rational(3, 4 * N))                 # all Re z > -3/(4N)
        # all Re z < -1/(4N): roots of p(w - 1/(4N)) in Re w < 0, i.e. Hurwitz of p(w - 1/(4N))
        w_ = sp.Symbol('w')
        P2 = sp.Poly(sp.expand(p.as_expr().subs(z, w_ - sp.Rational(1, 4 * N))), w_)
        co = P2.all_coeffs(); left = None
        if all(a > 0 for a in co):
            mins = hurwitz_minors(co); left = all(mm > 0 for mm in mins)
        else: left = False
        pos = positivity_tree(c, N)
        cmin = min(c[2:]) >= 1
        res.append({'N': N, 'W0': W0, 'strip_right': right, 'strip_left': left, 'pos_ok': pos['ok'], 'c_j>=1': cmin})
        print(res[-1], flush=True)
    # c~_j >= 1 for W = 1..3, N = 3..10 (Step 1 claim)
    cj = {}
    for N in range(3, 11):
        for W in [1, 2, 3]:
            c = rounded(N, W); cj['%d,%d' % (N, W)] = min(c[2:]) >= 1
    res.append({'c_j_ge_1_small_W': cj, 'all': all(cj.values())})
    print(res[-1], flush=True)
    out['S'] = res
    json.dump(out, open(OUT, 'w'), indent=1, default=str)

for N in [4, 5]:
    if 'U%d' % N in which:
        W = 2
        c = rounded(N, W)
        rec = instance(c, N); rec['N'] = N; rec['W'] = W
        # also W = 1 must fail C1 or positivity (prover chose the least W)
        c1 = rounded(N, 1); p1, G1, m1 = reduced_from_c(c1)
        rec['W1_C1'] = strip_right_of(p1, sp.Rational(1, N)); rec['W1_pos'] = positivity_tree(c1, N)['ok']
        out['U%d' % N] = rec
        print(json.dumps({k: v for k, v in rec.items() if k != 'minor_degrees_and_top_roots'}, default=str), flush=True)
        json.dump(out, open(OUT, 'w'), indent=1, default=str)
