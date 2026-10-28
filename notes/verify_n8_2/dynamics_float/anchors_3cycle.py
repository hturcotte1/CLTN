import sys, json
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/dynamics_float')
from dyn import *
from run_cands import traj_stats
C3 = [[int(j == (i + 1) % 3) for j in range(3)] for i in range(3)]
res = []
for eps, delta in ((F(1,4), F(1,2)), (F(1,100), F(29,2500))):
    W = weight_matrix(C3, eps, delta); x = fixed_point_exact(W); ev, Wf = jac_eigs(W)
    xf = np.array([float(v) for v in x]); q = delta*(1-eps)/eps
    r = dict(name='3cycle', eps=str(eps), delta=str(delta), q=str(q), legal=bool(q > 1), x=[str(v) for v in x], x_pos=bool(all(v > 0 for v in x)),
             eig_re=[float(v) for v in ev.real], max_re=float(max(ev.real)))
    ics = perturbed_ics(xf, 3, seed=7, scales=(0.3, 0.003), k=3)
    trajs = []
    for kind, x0 in ics:
        for T in (2000.0, 20000.0):
            s = traj_stats(W, Wf, simulate(Wf, x0, T, rtol=1e-8, atol=1e-11), xf, T, kind)
            trajs.append(dict(kind=kind, T=T, err_start=s['err_start'], err_end=s['err_end'], stayed_all_active=s['stayed_all_active'], end_support=s['end']['size'], end_residual=s['end']['residual'], min_x_end=s['min_x_end']))
    r['trajectories'] = trajs
    r['n_traj'] = len(trajs); r['n_converged'] = sum(t['err_end'] < 1e-7 for t in trajs); r['n_stayed_all_active'] = sum(t['stayed_all_active'] for t in trajs)
    r['verdict'] = 'does NOT converge to positive full-support fixed point (as expected)' if r['n_converged'] == 0 else 'UNEXPECTED convergence'
    print(json.dumps({k: v for k, v in r.items() if k != 'trajectories'}))
    print('  sample', trajs[1], trajs[-1])
    res.append(r)
json.dump(res, open('anchors_3cycle.json', 'w'), indent=1)
