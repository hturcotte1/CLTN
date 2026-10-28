"""Finalize: rebuild table and summary from per_graph.jsonl (complete run of run8.py) and re-run the controls."""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run8 import *
from collections import Counter
t0 = time.time()
lines = open(HITS).read().splitlines()
results = [json.loads(l) for l in open(f'{OUT}/per_graph.jsonl')]
assert len(results) == len(lines) and [r['d6'] for r in results] == [json.loads(l)['d6'] for l in lines]
# re-apply the (corrected) verdict rule of run8.verify_one to the stored per-graph records
def ok(row): return row is not None and row['legal'] and row['positive'] and row['stable'] == 'stable'
def bad(row): return row is not None and row['legal'] and (not row['positive'] or row['stable'] == 'unstable')
for r in results:
    pa, pb = r.get('point_a'), r.get('point_b')
    r['point_a_verdict'] = None if pa is None else ('stable_positive' if ok(pa) else ('unstable_or_nonpositive' if bad(pa) else 'undecided'))
    contradiction = r['complete'] or any('disagreement' in f or 'edge count' in f for f in r['failures'])
    if contradiction or bad(pb): v = 'REFUTED'
    elif ok(pb) and not r['failures']: v = 'CONFIRMED'
    else: v = 'UNDECIDED'
    r['notes'] = [n for n in r['notes'] if 'interior point' not in n]
    if v == 'CONFIRMED' and pa is None: r['notes'].append('no legal stored witness; confirmed on my interior point')
    if v == 'CONFIRMED' and pa is not None and not ok(pa): r['notes'].append('first legal stored witness is decided unstable, in agreement with its stored flag (finite-t point); confirmed on my interior point')
    r['verdict'] = v
extra = {x['d6']: x for x in json.load(open(f'{OUT}/extra53.json'))}
ctl = controls()
verd = Counter(r['verdict'] for r in results)
table = {}
for r in results:
    pb = r.get('point_b'); pa = r.get('point_a')
    table[r['d6']] = dict(point=dict(beta=pb['beta'], t=pb['t'], eps=pb['eps'], delta=pb['delta'], q=pb['q']) if pb else None,
                          x_positive=pb['positive'] if pb else None, hurwitz=pb['minor_signs'] if pb else None,
                          hurwitz_verdict=pb['stable'] if pb else None,
                          first_legal_witness=dict(index=r['first_legal_witness_index'], q=pa['q'], x_positive=pa['positive'], hurwitz=pa['stable']) if pa else None,
                          first_witness_illegal=r['first_witness_illegal'], witness_mismatches=r['witness_mismatches'],
                          point_a_verdict=r['point_a_verdict'],
                          stored_witness_beta_outside_refined_interval=(not extra[r['d6']]['mid_inside']) if r['d6'] in extra else False,
                          verdict=r['verdict'], failures=r['failures'], notes=r['notes'])
json.dump(table, open(f'{OUT}/v_fraction_n8_table.json', 'w'), indent=0)
npts = sum(len(r.get('witnesses', [])) + (1 if r.get('point_b') else 0) for r in results)
summary = dict(verifier='v_fraction', n_candidates=len(results), distinct_d6=len({r['d6'] for r in results}),
               verdicts=dict(verd), points_checked=npts,
               witness_points=sum(len(r.get('witnesses', [])) for r in results), interior_points=sum(1 for r in results if r.get('point_b')),
               first_witness_illegal=sum(1 for r in results if r.get('first_witness_illegal')),
               no_legal_stored_witness=sum(1 for r in results if r.get('point_a') is None),
               first_legal_witness_index=dict(Counter(str(r['first_legal_witness_index']) for r in results)),
               witness_disagreements=sum(len(r.get('witness_mismatches', [])) for r in results),
               witness_agree_stable_both_decided=sum(1 for r in results for w in r.get('witnesses', []) if w['agree_stable'] is True),
               witness_agree_positive=sum(1 for r in results for w in r.get('witnesses', []) if w['agree_positive'] is True),
               witness_undecided_mine=sum(1 for r in results for w in r.get('witnesses', []) if w['stable'] == 'undecided'),
               stored_unstable_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if not w['stored_stable']),
               my_unstable_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if w['stable'] == 'unstable'),
               legal_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if w['legal']),
               illegal_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if not w['legal']),
               interior_t_used=dict(Counter(str(r['point_b']['t']) for r in results if r.get('point_b'))),
               interior_strictly_inside=sum(1 for r in results if r.get('point_b') and r['point_b']['strictly_inside']),
               interior_stable_positive_legal=sum(1 for r in results if r.get('point_b') and r['point_b']['legal'] and r['point_b']['positive'] and r['point_b']['stable'] == 'stable'),
               interior_undecided=sum(1 for r in results if r.get('point_b') and r['point_b']['stable'] == 'undecided'),
               endpoint_notes=dict(Counter(r['interval']['lo_note'] + '|' + r['interval']['hi_note'] for r in results if 'interval' in r)),
               graphs_with_failures=sum(1 for r in results if r['failures']),
               graphs_with_notes=dict(Counter(n for r in results for n in r['notes'])),
               not_confirmed=[r['d6'] for r in results if r['verdict'] != 'CONFIRMED'],
               point_a_verdicts=dict(Counter(str(r['point_a_verdict']) for r in results)),
               point_a_unstable_agreeing_with_stored=[r['d6'] for r in results if r['point_a_verdict'] == 'unstable_or_nonpositive'],
               point_a_unstable_because_stored_beta_outside_interval=[d for d, x in extra.items() if not x['mid_inside']],
               point_a_unstable_only_t100_finite_t_effect=[d for d, x in extra.items() if x['mid_inside']],
               seconds_finalize=round(time.time() - t0, 1))
json.dump(dict(summary=summary, controls=ctl), open(f'{OUT}/v_fraction_n8_results.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
for r in results:
    if r['failures'] or r['verdict'] != 'CONFIRMED': print('NOTCONF', r['d6'], r['verdict'], r['failures'], r['notes'])
