"""Never runner: usage run_never.py <never_sample.json> <out.jsonl>. betas 2.05, 2.5, 4 at t=10000."""
import sys, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/dynamics_float')
from dyn import *
from run_cands import traj_stats
from multiprocessing import Pool
T_SHORT = 2000.0; T_LONG = 1e5
BETAS = [F('41/20'), F('5/2'), F(4)]

def one(c):
    t0 = time.time(); n, A = parse_d6(c['d6']); out = dict(c); out['complete_check'] = is_complete(A); cases = []
    for beta in BETAS:
        t = F(10000); eps = 1 / (t + 1); delta = (beta - 1) * eps; q = delta * (1 - eps) / eps
        W = weight_matrix(A, eps, delta); x = fixed_point_exact(W)
        cs = dict(beta=str(beta), eps=str(eps), delta=str(delta), q=str(q), legal=bool(q > 1))
        if x is None:
            cs['exact'] = dict(singular=True); xf = None
        else:
            cs['exact'] = dict(singular=False, all_pos=bool(all(v > 0 for v in x)), min_x=str(min(x)), x=[str(v) for v in x])
            xf = np.array([float(v) for v in x])
        ev, Wf = jac_eigs(W)
        cs['float'] = dict(max_re=float(max(ev.real)), eig_re=[float(v) for v in ev.real], eig_im=[float(v) for v in ev.imag])
        ref = xf if xf is not None else np.full(n, 1.0 / n)
        ics = perturbed_ics(np.abs(ref), n, seed=5000 + c['never_idx'], scales=(0.3, 0.003), k=5)
        trajs = [traj_stats(W, Wf, simulate(Wf, x0, T_SHORT), np.abs(ref), T_SHORT, kind) for kind, x0 in ics]
        trajs.append(traj_stats(W, Wf, simulate(Wf, ics[5][1], T_LONG), np.abs(ref), T_LONG, 'perturb0.003'))
        for r in trajs:
            r.pop('ode_vs_linear_rel', None)
        cs['n_traj'] = len(trajs)
        cs['n_end_full_support'] = sum(r['end']['size'] == n for r in trajs)
        cs['n_settled_full_support'] = sum(r['end']['size'] == n and r['end']['residual'] < 1e-9 for r in trajs)
        cs['n_ends_at_full_fp'] = sum(r['ends_at_full_fp'] for r in trajs)
        cs['end_supports'] = [r['end']['size'] for r in trajs]
        cs['end_residuals'] = [r['end']['residual'] for r in trajs]
        cs['end_exact_fp_match'] = [r['end']['exact_fp_match'] for r in trajs]
        cs['stayed_all_active'] = [r['stayed_all_active'] for r in trajs]
        cs['solver_ok'] = all(r['ok'] for r in trajs)
        reasons = []
        if cs['n_ends_at_full_fp'] > 0: reasons.append('trajectory_settled_at_positive_full_fp')
        if cs['n_settled_full_support'] > 0: reasons.append('trajectory_settled_full_support')
        if xf is not None and cs['exact']['all_pos'] and cs['float']['max_re'] < 0: reasons.append('positive_and_float_stable')
        if not cs['solver_ok']: reasons.append('solver_failure')
        cs['RED'] = bool(reasons); cs['red_reasons'] = reasons
        cases.append(cs)
    out['cases'] = cases; out['RED'] = any(cs['RED'] for cs in cases); out['wall'] = time.time() - t0
    return out

if __name__ == '__main__':
    nev = json.load(open(sys.argv[1]))
    if len(sys.argv) > 3: nev = nev[:int(sys.argv[3])]
    t0 = time.time()
    with open(sys.argv[2], 'w') as f, Pool(2) as pool:
        for i, out in enumerate(pool.imap(one, nev, chunksize=1)):
            f.write(json.dumps(out) + '\n'); f.flush()
            print(i, out['d6'], out['reason'], 'RED', out['RED'], [cs['red_reasons'] for cs in out['cases']], 'pos', [cs['exact'].get('all_pos') for cs in out['cases']], 'maxre', ['%.2g' % cs['float']['max_re'] for cs in out['cases']], 'wall', round(out['wall'], 1), 'elapsed', round(time.time() - t0), flush=True)
    print('done', time.time() - t0)
