"""
check_X2.py -- EXACT test of Lemma X2' on cases where positivity actually holds (check_lemmaP.py C3 drew only negative cases).
For 80 random draws from the 130 catalogue families (heur_catalogue_good.json: base, S) plus the multi-class families Q2, Y, a
random rational beta in (2, 3] and random class scale k in [1, 4]: prediction [w_G(beta) > 0 and sigma_j(beta) > 1] must equal
ctln.analyze_support's 'positive' on the explicit graph (t = 7; positivity is t-free).  Writes check_X2.json.
"""
import json, random
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

random.seed(11)
cat = json.load(open('heur_catalogue_good.json'))
pool = [(o['base'], [tuple(o['S'])], [1]) for o in cat] + [FAMILIES['Q2'], FAMILIES['Y']]
out = []; npos = 0
for i in range(80):
    base, classes, rhos = random.choice(pool)
    A = fl.from_d6(base)
    beta = F(random.randint(2001, 3000), 1000); kk = random.randint(1, 4)
    pos, w, sig = fl.positivity_core(A, classes, beta)
    G = fl.build(A, classes, [r * kk for r in rhos]); n = len(G)
    eps, delta = fl.eps_delta(beta, F(7))
    r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
    out.append({'base': base, 'classes': [list(S) for S in classes], 'rhos': rhos, 'k': kk, 'beta': str(beta), 'pred': pos,
                'ctln_positive': r['positive']})
    assert pos == r['positive'], out[-1]
    npos += pos
json.dump({'cases': out, 'n_positive': npos, 'n': len(out)}, open('check_X2.json', 'w'), indent=1)
print('X2 check passed:', npos, 'positive and', len(out) - npos, 'non-positive cases agree with ctln')
