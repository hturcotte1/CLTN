"""Positive controls: the engine must DETECT stable+positive points where they are known to exist."""
import json, sys
sys.path.insert(0, '.')
from fractions import Fraction as F
from engine import parse_d6, decide_point, decide_point_direct, is_domination_free
from alpha import classify_alpha
out = {}
# 1. FNFL (n = 7 finite-eps-only counterexample, alpha = -1/2 exactly, stable for t in [26, 53] at beta = 53/26)
A = parse_d6('&FNFL~[udo?')
r = decide_point(A, 40, 53, 26)
out['FNFL_t40_beta53/26'] = {'legal': r['legal'], 'stable': r['stable'], 'positive': r['positive'], 'alpha_class': classify_alpha(A)[0], 'domfree': is_domination_free(A)}
r2 = decide_point(A, 10000, 53, 26)
out['FNFL_t10000_beta53/26'] = {'legal': r2['legal'], 'stable': r2['stable'], 'detail': r2['stable_detail'], 'positive': r2['positive']}
# 2. 20 n = 8 limit candidates (every 450th HIT): certified interior witness from certificates/n8_candidates_legal_witnesses.json
W = {g['d6']: g for g in json.load(open('/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'))['graphs']}
hits = [json.loads(l)['d6'] for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')][::450]
res = []
for d6 in hits:
    w = W[d6]['witness_in_first_stable_cell']
    eps, delta = F(w['eps']), F(w['delta'])
    r = decide_point_direct(parse_d6(d6), eps, delta)
    res.append((d6, r['stable'], r['positive'], str(w.get('beta')), str(w.get('t'))))
out['hits_controls'] = res
out['hits_all_stable_positive'] = all(s is True and p is True for _, s, p, _, _ in res)
json.dump(out, open('controls.json', 'w'), indent=1); print(json.dumps(out)[:1500])
