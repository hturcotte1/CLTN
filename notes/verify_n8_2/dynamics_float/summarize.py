import json, sys, os
from collections import Counter
os.chdir('/home/user/CLTN/notes/verify_n8_2/dynamics_float')
S = {}
# candidates
C = []
for fn in ('candidates_out_A.jsonl', 'candidates_out_B.jsonl'):
    if os.path.exists(fn): C += [json.loads(l) for l in open(fn)]
configA = {json.loads(l)['d6'] for l in open('candidates_out_A.jsonl')} if os.path.exists('candidates_out_A.jsonl') else set()
cs = dict(n=len(C), n_every23=sum('every23' in c['why'] for c in C), n_null_stored_legal=sum('null_stored_legal' in c['why'] for c in C))
cs['by_t'] = dict(Counter(c['t'] for c in C)); cs['n_configA'] = len(configA); cs['n_configB'] = len(C) - len(configA)
cs['all_legal'] = all(c['legal'] for c in C); cs['n_complete'] = sum(c['complete_check'] for c in C)
cs['exact_x_positive'] = sum(c['exact'].get('all_pos', False) for c in C)
cs['float_maxre_negative'] = sum(c['float']['max_re'] < 0 for c in C if 'float' in C[0])
cs['slow_rate_quantiles'] = None
import numpy as np
rates = np.array([c['float']['slow_rate'] for c in C]); 
if len(C): cs['slow_rate_quantiles'] = {k: float(v) for k, v in zip(['min','q25','med','q75','max'], np.quantile(rates, [0,.25,.5,.75,1]))}
cs['horizon_sufficient'] = sum(c['horizon_sufficient'] for c in C)
cs['n_RED'] = sum(c['RED'] for c in C); cs['red_reasons'] = dict(Counter(r for c in C for r in c['red_reasons']))
cs['RED_list'] = [dict(idx=c['idx'], d6=c['d6'], t=c['t'], reasons=c['red_reasons'], rate=c['float']['slow_rate'], rate_x_Tlong=c.get('rate_x_Tlong', c['float']['slow_rate']*c['T_long'])) for c in C if c['RED']]
cs['horizon_unresolvable_rate_x_Tlong_lt_1'] = sum(1 for c in C if c['float']['slow_rate']*c['T_long'] < 1)
cs['decay_measurable_rate_x_Tlong_ge_1'] = sum(1 for c in C if c['float']['slow_rate']*c['T_long'] >= 1)
Lm = [r for c in C if c['float']['slow_rate']*c['T_long'] >= 1 for r in c['long'] if r['kind']=='perturb0.003']
cs['measurable_small_long_err_end_lt_start'] = [sum(r['err_end']<r['err_start'] for r in Lm), len(Lm)]
cs['measurable_small_long_err_end_lt_half'] = [sum(r['err_end']<r['err_half'] for r in Lm), len(Lm)]
if os.path.exists('followup.json'):
    cs['followup'] = json.load(open('followup.json'))
    for f in cs['followup']: f.pop('probes', None)
# short horizon stats (T=2000) over the specified 10 ICs (perturb0.3 + uniform) and small perturbations
def kinds(c, ks, key): return [r for r in c[key] if r['kind'] in ks]
spec = [r for c in C for r in kinds(c, ('perturb0.3', 'uniform'), 'short')]
small = [r for c in C for r in kinds(c, ('perturb0.003',), 'short')]
cs['short_T2000'] = dict(
    n_traj_specified=len(spec), converged_1e7_strict=sum(r['converged'] for r in spec),
    err_end_lt_1e7=sum(r['err_end_lt_1e7'] for r in spec), all_pos_end=sum(r['all_pos_end'] for r in spec),
    ends_all_active=sum(r['end']['size'] == 8 for r in spec), stayed_all_active=sum(r['stayed_all_active'] for r in spec),
    n_small_perturb=len(small), small_stayed_all_active=sum(r['stayed_all_active'] for r in small),
    small_ode_vs_linear_max=max([r.get('ode_vs_linear_rel', 0) for r in small] or [0]))
L = [r for c in C for r in c['long']]
Ls = [r for c in C for r in c['long'] if r['kind'] == 'perturb0.003']
Lsuff = [r for c in C if c['horizon_sufficient'] for r in c['long']]
cs['long'] = dict(n_traj=len(L), n_sufficient_horizon_traj=len(Lsuff),
    sufficient_converged_1e7=sum(r['err_end'] < 1e-7 for r in Lsuff),
    sufficient_converged_1e7_blocknoninc=sum(r['err_end'] < 1e-7 and r['block_noninc'] for r in Lsuff),
    sufficient_by_kind={k: [sum(1 for r in Lsuff if r['kind'] == k), sum(1 for r in Lsuff if r['kind'] == k and r['err_end'] < 1e-7)] for k in ('perturb0.3','perturb0.003','uniform')},
    small_stayed_all_active=sum(r['stayed_all_active'] for r in Ls), small_err_end_lt_start=sum(r['err_end'] < r['err_start'] for r in Ls),
    small_err_end_lt_half=sum(r['err_end'] < r['err_half'] for r in Ls),
    small_ode_vs_linear_max=max([r.get('ode_vs_linear_rel', 0) for r in Ls] or [0]),
    nonuniform_left_linear_regime_to_other_support=sum(1 for r in L if r['kind'] != 'uniform' and not r['stayed_all_active']),
    other_attractor_supports=dict(Counter(r['end']['size'] for r in L if r['end']['size'] != 8)),
    other_support_exact_fp_match=sum(1 for r in L if r['end']['size'] != 8 and r['end']['exact_fp_match']))
cs['wall_total_s'] = sum(c['wall'] for c in C)
S['candidates'] = cs
# never
N = [json.loads(l) for l in open('never_out.jsonl')] if os.path.exists('never_out.jsonl') else []
ns = dict(n=len(N), reasons=dict(Counter(n['reason'] for n in N)), n_complete=sum(n['complete_check'] for n in N))
per_beta = {}
for b in ('41/20', '5/2', '4'):
    cases = [cs_ for n in N for cs_ in n['cases'] if cs_['beta'] == b]
    per_beta[b] = dict(n=len(cases), legal=sum(c['legal'] for c in cases), x_singular=sum(c['exact'].get('singular') for c in cases),
        x_positive=sum(c['exact'].get('all_pos', False) for c in cases), float_maxre_negative=sum(c['float']['max_re'] < 0 for c in cases),
        positive_and_float_stable=sum(c['exact'].get('all_pos', False) and c['float']['max_re'] < 0 for c in cases),
        n_traj=sum(c['n_traj'] for c in cases), traj_ending_full_support=sum(c['n_end_full_support'] for c in cases),
        traj_settled_full_support_res1em9=sum(c['n_settled_full_support'] for c in cases), traj_at_positive_full_fp=sum(c['n_ends_at_full_fp'] for c in cases),
        RED=sum(c['RED'] for c in cases))
ns['per_beta'] = per_beta; ns['n_RED'] = sum(n['RED'] for n in N)
ns['RED_list'] = [dict(d6=n['d6'], never_idx=n['never_idx'], reasons=[c['red_reasons'] for c in n['cases']]) for n in N if n['RED']]
ns['never_full_records'] = [dict(d6=n['d6'], never_idx=n['never_idx'], per_beta=[dict(beta=c['beta'], x_pos=c['exact'].get('all_pos'), max_re=c['float']['max_re'], end_supports=c['end_supports']) for c in n['cases']]) for n in N if n['has_cert']]
ns['end_support_hist'] = dict(Counter(s for n in N for c in n['cases'] for s in c['end_supports']))
S['never'] = ns
if os.path.exists('anchors.json'):
    A = json.load(open('anchors.json'))
    S['anchors'] = [dict(name=a['name'], d6=a['d6'], eps=a['eps'], delta=a['delta'], q=a['q_check'], legal=a['legal'], x_pos=a['exact']['all_pos'], max_re=a['float']['max_re'],
        T_long=a['T_long'], horizon_sufficient=a['horizon_sufficient'], RED=a['red_reasons'], n_long_converged=a['n_long_converged'], n_long=len(a['long']),
        n_short_converged=a['n_short_converged'], long_err_end=[r['err_end'] for r in a['long']], long_end_supports=[r['end']['size'] for r in a['long']],
        short_end_supports=[r['end']['size'] for r in a['short']], expect=a['expect']) for a in A]
json.dump(S, open('summary.json', 'w'), indent=1)
print(json.dumps(S, indent=1))
