"""Extra check on the 53 graphs whose first legal stored witness is (stored and re-verified) unstable:
interval width, my beta vs stored midpoint beta, and the stored midpoint beta at larger t."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run8 import *
recs = {json.loads(l)['d6']: json.loads(l) for l in open(HITS)}
out = []
for l in open(f'{OUT}/per_graph.jsonl'):
    r = json.loads(l)
    if r['verdict'] != 'REFUTED': continue
    rec = recs[r['d6']]; A = decode_d6(r['d6'])
    lo_ub = F(r['interval']['lo_bounds'][1]); hi_lb = F(r['interval']['hi_bounds'][0])
    beta_mid = F(rec['witnesses'][0]['beta']); beta_mine = F(r['point_b']['beta'])
    rows = []
    for t in (10000, 100000, 10**6, 10**7):
        eps, delta = ed_from_bt(beta_mid, t)
        row, _ = point(A, eps, delta, f'midpoint_t{t}')
        rows.append((t, row['q'], row['legal'], row['positive'], row['stable']))
    for t in (100000, 10**6):
        eps, delta = ed_from_bt(beta_mine, t)
        row, _ = point(A, eps, delta, f'mine_t{t}')
        rows.append(('mine', t, row['q'], row['legal'], row['positive'], row['stable']))
    out.append(dict(d6=r['d6'], edges=r['edges'], interval=[S(lo_ub), S(hi_lb)], width=S(hi_lb - lo_ub), width_float_display=float(hi_lb - lo_ub),
                    beta_mid=S(beta_mid), beta_mine=S(beta_mine), mid_inside=bool(lo_ub < beta_mid < hi_lb),
                    stored_witness_stable=[w['stable'] for w in rec['witnesses']], first_legal=r['first_legal_witness_index'], checks=rows))
json.dump(out, open(f'{OUT}/extra53.json', 'w'), indent=1)
for o in out:
    print(o['d6'], o['edges'], 'width~%.2e' % o['width_float_display'], 'mid_inside', o['mid_inside'], o['stored_witness_stable'], o['first_legal'], o['checks'])
