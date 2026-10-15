"""
crosscheck_region.py -- randomized EXACT cross-checks of the region theorems (tests of the proofs, not the proofs).
 (X1) R family, region 'beta <= 43/20, all k >= 0': 25 random legal rational points (beta in (2, 43/20], t legal) and random k in
      [0, 30]: ctln.analyze_support must report legal, x > 0, stable (all Hurwitz minors > 0), not a clique.
 (X2) Q2 family, region 'beta <= 5/2, k >= 2': 10 random legal points, random k in [2, 6] (n = 15..31): same.
 (X3) Negative controls: R at random legal points with beta in [11/5, 14/5] and random k in [0, 30]: the pointwise Psi verdict
      (fixlib.stable_at_k) must equal ctln's verdict (both stable and unstable cases occur).
Writes crosscheck_region.json.
"""
import json, random
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

random.seed(7)
out = {'X1': [], 'X2': [], 'X3': []}


def rand_point(bmin, bmax):
    while True:
        beta = F(random.randint(int(bmin * 1000) + 1, int(bmax * 1000)), 1000)
        tmin = 1 / (beta - 2)
        t = tmin * F(random.randint(1001, 5000), 1000) if random.random() < 0.7 else tmin + F(random.randint(1, 10 ** 5), 10)
        if t * (beta - 2) > 1: return beta, t


for job, nm, bmin, bmax, kr, n_pts in (('X1', 'R', F(2), F(43, 20), (0, 30), 25), ('X2', 'Q2', F(2), F(5, 2), (2, 6), 10)):
    base, classes, rhos = FAMILIES[nm]; A = fl.from_d6(base)
    for _ in range(n_pts):
        beta, t = rand_point(bmin, bmax); kk = random.randint(*kr)
        rec = fl.certify_member(A, classes, [r * kk for r in rhos], beta, t)
        ok = rec['legal'] and rec['positive'] and rec['stable'] is True and not rec['clique']
        out[job].append({'beta': str(beta), 't': str(t), 'k': kk, 'n': rec['n'], 'ok': ok})
        print(job, beta, t, kk, ok, flush=True)
        assert ok
base, classes, rhos = FAMILIES['R']; A = fl.from_d6(base)
nst = 0
for _ in range(30):
    beta, t = rand_point(F(11, 5), F(14, 5)); kk = random.randint(0, 30)
    rec = fl.certify_member(A, classes, [kk], beta, t)
    v = fl.stable_at_k(A, classes, rhos, beta, t, kk)
    out['X3'].append({'beta': str(beta), 't': str(t), 'k': kk, 'ctln_stable': rec['stable'], 'psi': v[0]})
    print('X3', beta, t, kk, rec['stable'], v[0], flush=True)
    assert rec['stable'] == v[0] and rec['positive']
    nst += (v[0] is True)
out['X3_counts'] = {'stable': nst, 'unstable': 30 - nst}
json.dump(out, open('crosscheck_region.json', 'w'), indent=1)
print('ALL CROSS-CHECKS PASSED', out['X3_counts'])
