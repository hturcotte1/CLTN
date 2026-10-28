"""Exact check of the smallest q upper end: for &GMBQMHtMH~]W, limit stability at beta requires every root of
psi_D(z) = chi_D - chi_{D+J} in Re z > -1/beta (Lemma D1 / I.3).  Hurwitz test of psi_D(-x - 1/beta) at
beta = 2 + 1/10000 (q = 1.0001) and beta = 2 + 2/10000 (q = 1.0002).  Output: check_min_q.json"""
import json, os, sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
from ctln import hurwitz_stable
import lib
HERE = os.path.dirname(os.path.abspath(__file__))

def shifted(p, s):
    res = [F(0)]
    for c in p:
        out = [F(0)] * (len(res) + 1)
        for i, r in enumerate(res):
            out[i] += -r; out[i + 1] += -s * r
        out[-1] += c
        res = out
    while len(res) > 1 and res[0] == 0: res = res[1:]
    if res[0] < 0: res = [-r for r in res]
    return res

d6 = '&GMBQMHtMH~]W'
A = lib.from_d6(d6); D = lib.nonedge(A)
psi = [a - b for a, b in zip(lib.charpoly_int(D), lib.charpoly_int([[x + 1 for x in r] for r in D]))][1:]
out = {'d6': d6, 'psi_D': psi}
for beta in (F(20001, 10000), F(10001, 5000)):
    v, det = hurwitz_stable(shifted(psi, 1 / beta))
    out[f'beta={beta}, q={beta - 1}'] = [v, det]
json.dump(out, open(os.path.join(HERE, 'check_min_q.json'), 'w'), indent=1, default=str)
print(out)
