"""Follow-up on candidates flagged small_perturb_left_linear_regime: tinier perturbations and linear transient amplification."""
import sys, json, os
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/dynamics_float')
from dyn import *
from scipy.linalg import expm
os.chdir('/home/user/CLTN/notes/verify_n8_2/dynamics_float')
recs = []
for fn in ('candidates_out_A.jsonl', 'candidates_out_B.jsonl'):
    if os.path.exists(fn):
        recs += [json.loads(l) for l in open(fn)]
flag = [c for c in recs if 'small_perturb_left_linear_regime' in c['red_reasons'] or 'expected_convergence_missing_small' in c['red_reasons'] or 'stayed_linear_but_not_at_fp' in c['red_reasons'] or 'ode_disagrees_with_linear_solution' in c['red_reasons']]
print('flagged', len(flag))
out = []
for c in flag:
    n, A = parse_d6(c['d6']); W = weight_matrix(A, F(c['eps']), F(c['delta'])); x = fixed_point_exact(W)
    xf = np.array([float(v) for v in x]); ev, Wf = jac_eigs(W); rate = -max(ev.real)
    J = -np.eye(n) + Wf
    # worst-case linear amplification over a log grid of times
    ts = np.unique(np.concatenate([np.logspace(0, np.log10(max(10 / rate, 10)), 60)]))
    amp = max(np.linalg.norm(expm(J * t), np.inf) for t in ts)
    T = min(30 / rate, 6e5)
    res = dict(d6=c['d6'], idx=c['idx'], t=c['t'], rate=rate, linear_amplification_max=float(amp), T=T, rate_x_T=rate * T, probes=[])
    for s in (3e-4, 3e-5, 3e-6):
        rng = np.random.default_rng(1000 + c['idx'])
        for k in range(5):
            u = rng.uniform(-1, 1, n); x0 = xf * (1 + s * u)
            sol = simulate(Wf, x0, T)
            X = sol.y.T; inp = X @ Wf.T + 1.0
            err = np.max(np.abs(X - xf), axis=1) / np.max(np.abs(xf))
            res['probes'].append(dict(scale=s, k=k, stayed_all_active=bool(np.all(inp > 0)), err_start=float(err[0]), err_max=float(err.max()),
                                      err_end=float(err[-1]), amplification_observed=float(err.max() / err[0]), converged_1e7=bool(err[-1] < 1e-7)))
    pr = res['probes']
    res['n_probes'] = len(pr); res['n_stayed'] = sum(p['stayed_all_active'] for p in pr); res['n_converged'] = sum(p['converged_1e7'] for p in pr)
    res['n_err_end_lt_start'] = sum(p['err_end'] < p['err_start'] for p in pr)
    res['tiny_scale_ok'] = all(p['stayed_all_active'] and p['err_end'] < p['err_start'] for p in pr if p['scale'] == 3e-6)
    out.append(res)
    print(c['d6'], c['t'], 'rate %.3g' % rate, 'amp %.3g' % amp, 'rate*T %.2f' % (rate * T), 'stayed', res['n_stayed'], '/', res['n_probes'], 'conv', res['n_converged'], 'end<start', res['n_err_end_lt_start'], 'tiny_ok', res['tiny_scale_ok'], flush=True)
json.dump(out, open('followup.json', 'w'), indent=1)
