import sys, json
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from fractions import Fraction as F
R = '/home/user/CLTN/'
n7 = {json.loads(l)['d6']: json.loads(l) for l in open(R + 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')}
hits = {json.loads(l)['d6']: json.loads(l) for l in open(R + 'n8/n8_limit_HITS.jsonl')}
P = n7['&F^v|^Fo\\B?']; H = hits['&G]O|FyZBz~\\w']; B = d6_decode('&G^z~Vwn?zBRw'); A = d6_decode(P['d6'])
assert find_isomorphism(B, d6_decode(H['d6'])) is not None
out = []
for name, iv in (('parent_mid', P['cert']['stable_intervals'][0]), ('hit_mid', H['cert']['stable_intervals'][0])):
    beta, lob, hib, steps = interior_point(iv, extra_steps=4)
    for t in (10000, 100000):
        eps = F(1, t + 1); delta = (beta - 1) * eps
        dB = decide(B, eps, delta); dA = decide(A, eps, delta)
        out.append({'point': name, 'beta': str(beta), 't': t, 'blowup': [dB['stability'], dB['positive']], 'parent': [dA['stability'], dA['positive']]})
print(json.dumps(out, indent=1)); json.dump(out, open('part5c_disjoint.json', 'w'), indent=1)
