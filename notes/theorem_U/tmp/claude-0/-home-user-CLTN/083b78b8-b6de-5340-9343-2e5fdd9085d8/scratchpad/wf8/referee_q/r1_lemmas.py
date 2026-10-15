"""r1_lemmas.py -- referee: exact checks of Lemma T, Lemma Q (random in-trees + arms forests) and Lemma P (closed form vs
direct exact solve), independent implementation. Writes r1_lemmas.json."""
import sys, json, random, time
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q')
sys.path.insert(0, '/home/user/CLTN')
from ref_lib import *
import ctln
from fastregion import icharpoly

out = {'T_Q_random': [], 'T_Q_arms': [], 'P': []}
rng = random.Random(20261015)


def check_TQ(D, root, use_sympy=True):
    m = len(D)
    Lam, Gam, N = Lam_Gam(D, root)
    L = [[D[i][j] - (1 if i == root else 0) for j in range(m)] for i in range(m)]
    if use_sympy:
        chiL = charpoly_sym(L); chiD = charpoly_sym(D); ps = psi_sym(D)
    else:
        chiL = sp.Poly(icharpoly(L), z); chiD = sp.Poly(icharpoly(D), z)
        DJ = [[D[i][j] + 1 for j in range(m)] for i in range(m)]
        ps = sp.Poly(sp.Poly(icharpoly(D), z).as_expr() - sp.Poly(icharpoly(DJ), z).as_expr(), z)
    T1 = (chiL == to_z_poly(Lam, m))
    T2 = (chiD == to_z_poly(1 - Gam, m))
    num = sp.expand((m + 1 - len(N)) * (1 - Gam) + (len(N) - 1) * Lam)
    Qp, rem = sp.div(sp.Poly(num, x), sp.Poly(1 - x, x))
    Qok = rem.is_zero and (ps == to_z_poly(Qp.as_expr(), m - 1))
    return T1, T2, Qok


for trial in range(30):
    m = rng.randrange(4, 19)
    D, root = random_intree(m, rng)
    T1, T2, Qok = check_TQ(D, root)
    out['T_Q_random'].append({'m': m, 'T_L': T1, 'T_D': T2, 'Q': Qok})
print('random in-trees: all ok =', all(r['T_L'] and r['T_D'] and r['Q'] for r in out['T_Q_random']), flush=True)

arms_cases = [[1, 1, 1], [1, 1, 2, 1], [1, 1, 5, 3, 1], [1, 1, 4, 2, 1], [1, 1, 2, 1, 1, 1, 1], [1, 1, 3, 0, 2, 0, 1],
              [1, 1, 1, 2, 0, 1]]
for c in arms_cases:
    D, idx = arms_core(c)
    T1, T2, Qok = check_TQ(D, 0)
    p, G, m = reduced_from_c(c)
    # reduced formula vs explicit: chi_L = z^(m-d) p and psi = z^(m-1) Q(1/z) = z^(m-1-d) G
    d = len(c) - 1
    L = [[D[i][j] - (1 if i == 0 else 0) for j in range(len(D))] for i in range(len(D))]
    red_ok = (charpoly_sym(L) == sp.Poly(z ** (m - d) * p.as_expr(), z)) and (psi_sym(D) == sp.Poly(z ** (m - 1 - d) * G.as_expr(), z))
    out['T_Q_arms'].append({'c': c, 'm': len(D), 'T_L': T1, 'T_D': T2, 'Q': Qok, 'reduced_formula': red_ok})
    print(c, len(D), T1, T2, Qok, red_ok, flush=True)

# the 65-node core with the repo integer charpoly (second implementation), m = 65
t0 = time.time()
c65 = [1, 1, 11, 7, 10, 3, 1]
D, idx = arms_core(c65)
T1, T2, Qok = check_TQ(D, 0, use_sympy=False)
p, G, m = reduced_from_c(c65)
out['T_Q_arms'].append({'c': c65, 'm': len(D), 'T_L': T1, 'T_D': T2, 'Q': Qok, 'm_formula': m, 'secs': round(time.time() - t0, 1)})
print('65-core', len(D), m, T1, T2, Qok, flush=True)


# Lemma P: closed form vs direct exact solve; referee re-derives the cycle system with sympy
lam_, e0, e1, e2, r_, F0s, F1s, F2s = sp.symbols('lam e0 e1 e2 r F0 F1 F2')
sol = sp.solve([sp.Eq(e2, r_ * (1 - lam_)), sp.Eq(e1, r_ * (1 - e2) + lam_ * F2s / r_),
                sp.Eq(e0, r_ * (1 - e1) + lam_ * F1s / r_), sp.Eq(e2, r_ * (1 - e0) + lam_ * F0s / r_)],
               [lam_, e0, e1, e2], dict=True)[0]
closed = r_ ** 3 * (1 - r_ + r_ ** 2) / (F0s - r_ * F1s + r_ ** 2 * F2s + r_ ** 2 + r_ ** 5)
out['P_symbolic_closed_form_matches'] = bool(sp.simplify(sol[lam_] - closed) == 0)
print('Lemma P closed form symbolic:', out['P_symbolic_closed_form_matches'], flush=True)


def direct_u(D, beta):
    m = len(D); r = 1 / F(beta)
    M = [[F(D[j][i]) + (r if i == j else 0) for j in range(m)] for i in range(m)]
    return ctln.solve_exact(M, [F(1)] * m)


for c, beta in [([1, 1, 5, 3, 1], 3), ([1, 1, 4, 2, 1], F(13, 4)), (c65, 5), (c65, F(28, 5)), (c65, 4), (c65, F(11, 2)),
                ([1, 1, 2, 1, 1, 1, 1], 3)]:
    D, idx = arms_core(c)
    u = direct_u(D, beta)
    r = 1 / F(beta); d = len(c) - 1
    F0 = sum(F(c[j]) * r ** (4 - j) for j in range(4, d + 1, 2))
    F1 = sum(F(c[j]) * r ** (3 - j) for j in range(3, d + 1, 2))
    F2 = F(c[2] - 1)
    lam_cf = r ** 3 * (1 - r + r * r) / (F0 - r * F1 + r * r * F2 + r * r + r ** 5)
    lam_direct = 1 - (1 - u[idx['c2']]) / r
    U = sum(u)
    w_pos = min(u) > 0 and U > 1
    sigmaS_gt1 = u[idx['c2']] < 1  # equivalent to sum_{i != c2} w_i > 1 when U > 1
    # direct check of sigma_S > 1 with w = u/(U-1)
    w = [ui / (U - 1) for ui in u]
    sigS = sum(w) - w[idx['c2']]
    out['P'].append({'c': c, 'beta': str(beta), 'lam_closed': str(lam_cf), 'lam_direct': str(lam_direct),
                     'match': lam_cf == lam_direct, 'u_min': str(min(u)), 'U': str(U), 'w_pos': w_pos,
                     'sigma_S': str(sigS), 'sigma_S_gt_1': sigS > 1, 'u_c2_lt_1': sigmaS_gt1})
    print('P', c, beta, lam_cf == lam_direct, w_pos, sigS > 1, flush=True)

json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q/r1_lemmas.json', 'w'), indent=1, default=str)
