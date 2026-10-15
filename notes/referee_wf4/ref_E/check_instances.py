"""Referee: per-instance exact checks for E_a, a = 0..5.
(1) independent blow-up equals the certified adjacency; E_0 isomorphic to the FNFL candidate certificate graph;
(2) n chi_perp: (l^2+l+1) divides exactly once, cofactor roots all Re < -1/2 (exact Hurwitz, all minors > 0) => alpha = -1/2 exactly;
(3) psi_A = (z-1)^(n-7) Phi(z-1) with Phi the closed form;
(4) Lemma B hypotheses (i)-(iv) exactly per instance (R via exact Q(omega) arithmetic on explicit chi, psi);
(5) witnesses re-run with ctln.analyze_support AND an independent K-based check."""
import sys, json, itertools
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E')
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
from fractions import Fraction as F
from reflib import *
import ctln

W = json.load(open('/home/user/CLTN/notes/research_wf3/theory/family_E_witnesses.json'))
H = d6_decode('&FNDL~[uds?')
AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
res = []
for av in range(AMAX + 1):
    rec = {'a': av}
    A, cls = blowup(H, [1, 1, 2, av, 1, 1, 1]); n = len(A)
    wr = W[av]
    rec['adjacency_matches_json'] = [''.join(map(str, r)) for r in A] == wr['A_rows']
    rec['d6'] = d6_encode(A); rec['d6_matches_json'] = rec['d6'] == wr['d6']
    rec['not_clique'] = any(A[i][j] == 0 for i in range(n) for j in range(n) if i != j)
    # (2) alpha
    P = sp.Poly(nchi_perp(A), l)
    mult = 0; Pq = P
    while True:
        q_, r_ = sp.div(Pq, sp.Poly(l**2 + l + 1, l))
        if r_.is_zero: Pq = q_; mult += 1
        else: break
    rec['phi3_mult'] = mult
    v, mins = all_roots_left_of(Pq.as_expr(), l, sp.Rational(-1, 2))
    rec['cofactor_all_Re_lt_-1/2'] = v
    rec['cofactor_minors_positive'] = all(x > 0 for x in mins)
    rec['alpha_eq_-1/2'] = (mult >= 1 and v is True)
    # (3) psi identity
    chiA, psiA = chi_psi_A(A, z)
    Phi = z*(z**2+z+1)*((av+7)*z*(z**2+z+1) + 2)
    rec['psi_eq_(z-1)^(n-7)Phi(z-1)'] = sp.expand(psiA - (z-1)**(n-7)*Phi.subs(z, z-1)) == 0
    # also perp poly vs psi: nchi_perp(l) relation psi_A(z) = (-1)^(n-1) nchi_perp(-z)
    rec['psi_vs_perp'] = sp.expand(psiA - (-1)**(n-1)*P.as_expr().subs(l, -z)) == 0
    # (4) Lemma B hypotheses
    pr = red(psiA, z); rec['i_psi(z0)=0'] = (pr == (0, 0))
    cr = red(chiA, z); rec['i_chi(z0)!=0'] = (cr != (0, 0))
    dpr = red(sp.diff(psiA, z), z); rec['i_simple'] = (dpr != (0, 0))
    # (ii): other roots of psi Re > 1/2  <=>  cofactor of perp poly Re < -1/2 (already) ; check directly too
    psi_co = sp.quo(sp.Poly(psiA, z), sp.Poly(z**2 - z + 1, z))
    v2, _ = all_roots_left_of(sp.expand(-psi_co.as_expr().subs(z, -z)), z, sp.Rational(-1, 2))
    rec['ii_other_psi_roots_Re_gt_1/2'] = v2 and sp.rem(psi_co, sp.Poly(z**2 - z + 1, z)).is_zero is False
    val = qmul(cr, qinv(dpr)); R = qre(val)
    rec['iii_R'] = str(R); rec['iii_R_gt_1/8'] = bool(R > sp.Rational(1, 8))
    rec['inv_sprime'] = [str(val[0]), str(val[1])]
    # (iv) w(2) = (A^T + I/2)^{-1} 1 > 0 ; also W = 1^T w > 1
    Mt = [[F(A[j][i]) + (F(1, 2) if i == j else 0) for j in range(n)] for i in range(n)]
    w = ctln.solve_exact(Mt, [F(1)]*n)
    rec['iv_w2_positive'] = w is not None and all(x > 0 for x in w)
    rec['iv_W_total'] = str(sum(w)); rec['iv_W_gt_1'] = sum(w) > 1
    claimed_y = [F(k, 18*av+17) for k in (6, 2, 2, 18, 6, 2, 10)]
    rec['iv_w_equals_class_vector'] = all(w[i] == claimed_y[cls[i]] for i in range(n))
    # (5) witness
    wt = wr['witness']; eps = F(wt['eps']); delta = F(wt['delta'])
    r = ctln.analyze_support(A, list(range(n)), eps, delta)
    rec['w_legal'] = ctln.legal(eps, delta); rec['w_q'] = str(ctln.q_of(eps, delta))
    t = (1 - eps)/eps; beta = 1 + delta/eps
    rec['w_t_beta'] = [str(t), str(beta)]; rec['w_legal_curve'] = t*(beta - 2) > 1
    rec['w_positive'] = r['positive']; rec['w_stable'] = r['stable']; rec['w_detail'] = r['stability_detail']
    rec['w_sfp'] = r['stable_fixed_point']; rec['w_target_free'] = r['target_free']
    rec['w_x_matches_json'] = r['x'] == wt['x']
    # independent: K = I + beta N + tJ, N = J - I - A^T; I - W == eps K ; K positive stable; K^{-1} 1 > 0
    Wm = ctln.weight_matrix(A, eps, delta)
    K = [[(1 if i == j else 0) + beta*((1 if i != j else 0) - A[j][i]) + t for j in range(n)] for i in range(n)]
    rec['IminusW_eq_epsK'] = all((F(1) if i == j else F(0)) - Wm[i][j] == eps*K[i][j] for i in range(n) for j in range(n))
    xk = ctln.solve_exact(K, [F(1)]*n); rec['indep_K_x_pos'] = all(x > 0 for x in xk)
    cp = sp.Matrix([[sp.Rational(e.numerator, e.denominator) for e in row] for row in K]).charpoly(z).as_expr()
    # K positive stable <=> roots of charpoly(K)(-y) Re<0
    v3, mins3 = all_roots_left_of(sp.expand(cp.subs(z, -z)), z, 0)
    rec['indep_K_positive_stable'] = v3; rec['indep_K_minors_all_pos'] = all(x > 0 for x in mins3)
    res.append(rec); print(json.dumps(rec), flush=True)
json.dump(res, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E/check_instances.json', 'w'), indent=1)
