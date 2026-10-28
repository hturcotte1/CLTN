import json, sys, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8/v_structure')
from exact8 import *
from verify import refine_root_interval
OUT = '/home/user/CLTN/notes/verify_n8/v_structure/'
recs = [json.loads(l) for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
D = json.load(open(OUT + 'verify_out.json')); rows = {r['d6']: r for r in D['rows']}
t0 = time.time(); res = []
for rec in recs:
    hi = rec['cert']['stable_intervals'][0]['hi']; p = hi['beta_minpoly_hi_to_lo']
    a, b = map(F, hi['beta_interval'])
    a2, b2 = refine_root_interval(a, b, p, steps=20)
    res.append({'d6': rec['d6'], 'stored_q_upper_bound': str(b - 1), 'q_hi_refined_interval': [str(a2 - 1), str(b2 - 1)], 'q_hi_float': float((a2 + b2) / 2 - 1),
                'deg': len(p) - 1, 'lo_kind': rec['cert']['stable_intervals'][0]['lo']['kind'], 'in_37_no_legal_stored': rows[rec['d6']]['stored_legal_witness'] is None,
                '_lo': a2 - 1, '_hi': b2 - 1})
print('refined %d in %.1fs' % (len(res), time.time() - t0))
res.sort(key=lambda r: r['_hi'], reverse=True)
top = res[:10]
# ordering certified if consecutive refined intervals are disjoint
cert = all(top[i]['_lo'] > top[i + 1]['_hi'] for i in range(9)) and top[9]['_lo'] > res[10]['_hi']
largest = {'d6': top[0]['d6'], 'q_upper_end_refined_rational_bounds': top[0]['q_hi_refined_interval'], 'float': top[0]['q_hi_float'], 'deg': top[0]['deg']}
sec = [r for r in res if r['_hi'] >= F(3, 2)]
with open(OUT + 'q_upper_ends.txt', 'a') as f:
    f.write('\n--- after 20 exact bisection steps of every hi isolating interval (stored bounds are coarse, e.g. [2,3] -> [8495/4096, 531/256]) ---\n')
    f.write('largest stable-q upper end: %s in [%s, %s] (~%.6f), minpoly degree %d\n' % (top[0]['d6'], top[0]['q_hi_refined_interval'][0], top[0]['q_hi_refined_interval'][1], top[0]['q_hi_float'], top[0]['deg']))
    f.write('ten largest (ranking certified by disjoint refined intervals: %s):\n' % cert)
    for r in top: f.write('  %-16s ~%.6f  in [%s, %s]  deg %d  lo %s\n' % (r['d6'], r['q_hi_float'], r['q_hi_refined_interval'][0], r['q_hi_refined_interval'][1], r['deg'], r['lo_kind']))
    f.write('all refined upper ends < 2: %s; number with upper end >= 3/2: %d; max refined upper bound rational: %s\n' % (all(r['_hi'] < 2 for r in res), len(sec), str(top[0]['_hi'])))
for r in res: r.pop('_lo'); r.pop('_hi')
json.dump({'largest': largest, 'ten_largest': top, 'ordering_certified': cert, 'all_below_2': all(F(r['q_hi_refined_interval'][1]) < 2 for r in res), 'all': res}, open(OUT + 'q_upper_ends_refined.json', 'w'), indent=1)
print(open(OUT + 'q_upper_ends.txt').read())
print('top-10 among the 37 w/o legal stored witness:', [r['d6'] for r in top if r['in_37_no_legal_stored']])
print('witness betas of the two deg-12 graphs:', [(rec['d6'], rec['witnesses'][0]['beta'], rec['witnesses'][0]['stable']) for rec in recs if rec['d6'] in ('&GUru\\Aw^qz]w', '&GZr[^`\\VJ}Iw')])
