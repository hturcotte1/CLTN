"""Candidate runner. usage: run_cands.py <sample.json> <out.jsonl> [limit]
short: T=2000, 20 ICs (5 each perturb 0.3/0.03/0.003 + 5 uniform[0,2]^8).
long: T_long = 30/rate if <= 2e6 else 5e5 (horizon-limited), ICs perturb0.3#0, perturb0.003#0, uniform#0."""
import sys, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/dynamics_float')
from dyn import *
from multiprocessing import Pool

T_SHORT = 2000.0
import os
CAP_FULL = float(os.environ.get('CAP_FULL', 2e6))
CAP_LIMITED = float(os.environ.get('CAP_LIMITED', 5e5))
SCALES = tuple(float(v) for v in os.environ.get('SCALES', '0.3,0.03,0.003').split(','))
LONG_KINDS = os.environ.get('LONG_KINDS', 'perturb0.3,perturb0.003,uniform').split(',')

def traj_stats(W, Wf, sol, xf, T, kind):
    r = analyze_trajectory(sol, xf)
    X = sol.y.T
    scale = max(np.max(np.abs(xf)), 1e-300)
    err = np.max(np.abs(X - xf[None, :]), axis=1) / scale
    m = len(err); half = err[m // 2:]
    bl = np.array_split(half, 8); bm = [float(np.max(b)) for b in bl]
    r['block_noninc'] = bool(all(bm[i + 1] <= bm[i] * (1 + 1e-6) + 1e-15 for i in range(7)))
    r['err_start'] = float(err[0]); r['decreased_end_lt_half'] = bool(err[-1] < half[0])
    r['T'] = T; r['kind'] = kind
    inp = (X @ Wf.T) + 1.0
    r['min_input_along'] = float(np.min(inp))          # <0 somewhere => some unit was off at a sample time
    r['stayed_all_active'] = bool(np.all(inp > 0))
    r['end'] = classify_end(W, Wf, X[-1])
    r['ends_at_full_fp'] = bool(r['end']['size'] == len(xf) and r['end']['exact_fp_match'])
    if r['stayed_all_active']:
        xl = linear_solution(Wf, xf, X[0], T)
        r['ode_vs_linear_rel'] = float(np.max(np.abs(X[-1] - xl)) / scale)
    del r['x_end']
    return r

def one(c):
    t0 = time.time()
    n, A = parse_d6(c['d6'])
    eps, delta = F(c['eps']), F(c['delta'])
    q = delta * (1 - eps) / eps
    W = weight_matrix(A, eps, delta)
    x = fixed_point_exact(W)
    out = dict(c)
    out['q_check'] = str(q); out['legal'] = bool(eps > 0 and delta > 0 and q > 1)
    out['complete_check'] = is_complete(A)
    if x is None:
        out['exact'] = dict(singular=True); out['RED'] = True; out['red_reasons'] = ['singular']; return out
    out['exact'] = dict(singular=False, x=[str(v) for v in x], all_pos=bool(all(v > 0 for v in x)), min_x=str(min(x)))
    ev, Wf = jac_eigs(W)
    rate = float(-max(ev.real))
    out['float'] = dict(eig_re=[float(v) for v in ev.real], eig_im=[float(v) for v in ev.imag],
                        max_re=float(max(ev.real)), slow_rate=rate)
    xf = np.array([float(v) for v in x])
    ics = perturbed_ics(xf, n, seed=1000 + c.get('idx', c.get('never_idx', 0)), scales=SCALES)
    short = [traj_stats(W, Wf, simulate(Wf, x0, T_SHORT), xf, T_SHORT, kind) for kind, x0 in ics]
    if rate > 0 and 30.0 / rate <= CAP_FULL:
        T_long = max(30.0 / rate, T_SHORT); suff = True
    else:
        T_long = CAP_LIMITED; suff = False
    long = [traj_stats(W, Wf, simulate(Wf, x0, T_long), xf, T_long, kind) for kind, x0 in [next((k, x0) for k, x0 in ics if k == lk) for lk in LONG_KINDS]]
    out['short'] = short; out['long'] = long; out['T_long'] = T_long; out['horizon_sufficient'] = suff
    small = [r for r in short + long if r['kind'] == 'perturb0.003']
    reasons = []
    if not out['exact']['all_pos']: reasons.append('exact_x_not_positive')
    if out['float']['max_re'] >= 0: reasons.append('float_eig_nonneg')
    if not all(r['ok'] for r in short + long): reasons.append('solver_failure')
    if not all(r['stayed_all_active'] for r in small): reasons.append('small_perturb_left_linear_regime')
    if not all(r.get('ode_vs_linear_rel', 0.0) < 1e-6 for r in small): reasons.append('ode_disagrees_with_linear_solution')
    if not suff and rate * T_long >= 1.0 and not all(r['err_end'] < r['err_start'] for r in long if r['kind'] == 'perturb0.003'): reasons.append('limited_horizon_small_perturb_err_not_below_start')
    out['rate_x_Tlong'] = rate * T_long
    if suff and not all(r['err_end'] < 1e-7 for r in long if r['kind'] == 'perturb0.003'): reasons.append('expected_convergence_missing_small')
    if suff and not all(r['ends_at_full_fp'] for r in long if r['kind'] != 'uniform' and r['stayed_all_active']): reasons.append('stayed_linear_but_not_at_fp')
    out['RED'] = bool(reasons); out['red_reasons'] = reasons
    out['n_short_converged'] = sum(r['converged'] for r in short)
    out['n_short_ends_full_fp'] = sum(r['ends_at_full_fp'] for r in short)
    out['n_long_converged'] = sum(r['err_end'] < 1e-7 and r['block_noninc'] for r in long)
    out['n_long_ends_full_fp'] = sum(r['ends_at_full_fp'] for r in long)
    out['wall'] = time.time() - t0
    return out

if __name__ == '__main__':
    cands = json.load(open(sys.argv[1]))
    if len(sys.argv) > 3: cands = cands[:int(sys.argv[3])]
    t0 = time.time()
    with open(sys.argv[2], 'w') as f, Pool(2) as pool:
        for i, out in enumerate(pool.imap(one, cands, chunksize=1)):
            f.write(json.dumps(out) + '\n'); f.flush()
            print(i, out['d6'], 't', out.get('t'), 'rate', '%.3g' % out['float']['slow_rate'], 'T_long', out['T_long'],
                  'RED', out['red_reasons'], 'longconv', out['n_long_converged'], 'wall', round(out['wall'], 1), 'elapsed', round(time.time() - t0), flush=True)
    print('done', time.time() - t0)
