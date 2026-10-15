"""
certify_members.py -- exact instance certificates (ctln.analyze_support on the explicit graphs) for fixed-parameter families.
Usage: python3 certify_members.py <job> ...      jobs defined in JOBS below; writes cert_<job>.json (one record per member).
Each record: legal, positive (x > 0), stable (all Hurwitz minors of the exact charpoly of -I + W positive), not clique,
plus the Psi-based Hurwitz verdict (fixlib.stable_at_k) as an independent cross-check, digraph6 of the explicit graph,
and sha256 of the exact characteristic polynomial (the polynomial itself is stored for n <= 30).
These members are CANDIDATE counterexamples pending human hand-verification (CLAUDE.md section 3).
"""
import json, sys, time, hashlib
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

JOBS = {
    # family, beta, t, ks
    'R_eps1_21': ('R', F(21, 10), F(20), [0, 1, 2, 5, 10, 20, 100]),
    'R_eps1_8': ('R', F(43, 20), F(7), [0, 1, 2, 5, 10, 20, 100]),
    'Q2_eps1_2': ('Q2', F(13, 4), F(1), [9, 10, 12]),   # pointwise threshold at this point: k0 = 9 (t4/pointwise_thresholds)
    'Y_eps1_2': ('Y', F(13, 4), F(1), [31]),
}


def run(job, ks_override=None):
    name, beta, t, ks = JOBS[job]
    base, classes, rhos = FAMILIES[name]
    A = fl.from_d6(base)
    if ks_override: ks = ks_override
    eps, delta = fl.eps_delta(beta, t)
    out = {'job': job, 'family': name, 'base': base, 'classes': [list(S) for S in classes], 'rhos': rhos,
           'beta': str(beta), 't': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(fl.ctln.q_of(eps, delta)),
           'legal': fl.ctln.legal(eps, delta), 'members': []}
    assert out['legal']
    for kk in ks:
        t0 = time.time()
        sizes = [rho * kk for rho in rhos]
        rec = fl.certify_member(A, classes, sizes, beta, t)
        psi_v = fl.stable_at_k(A, classes, rhos, beta, t, kk)
        cp = rec.pop('charpoly')
        h = hashlib.sha256(json.dumps(cp).encode()).hexdigest()
        r = {'k': kk, 'sizes': sizes, 'n': rec['n'], 'd6': rec['d6'], 'legal': rec['legal'], 'positive': rec['positive'],
             'stable': rec['stable'], 'stability_detail': rec['stability_detail'], 'clique': rec['clique'],
             'stable_fixed_point': rec['stable_fixed_point'], 'psi_hurwitz_verdict': psi_v[0], 'charpoly_sha256': h,
             'x': rec['x'], 'seconds': round(time.time() - t0, 1)}
        if rec['n'] <= 30: r['charpoly'] = cp
        ok = r['legal'] and r['positive'] and r['stable'] is True and not r['clique'] and r['stable_fixed_point'] and psi_v[0] is True
        r['CERTIFIED'] = bool(ok)
        out['members'].append(r)
        print(job, 'k', kk, 'n', r['n'], 'certified', ok, r['stability_detail'], r['seconds'], 's', flush=True)
        json.dump(out, open('cert_%s.json' % job, 'w'), indent=1)
    return out


if __name__ == '__main__':
    job = sys.argv[1]
    ks = [int(v) for v in sys.argv[2:]] or None
    run(job, ks)
