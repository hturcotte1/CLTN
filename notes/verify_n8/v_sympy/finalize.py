"""Merge the chunk tables, cross-check coverage, and write summary + table. Usage: python3 -I finalize.py <dir>"""
import sys, json
from collections import Counter
D = sys.argv[1]
rows = []
for c in range(5):
    rows += [json.loads(l) for l in open(f'{D}/table_part{c}.jsonl')]
rows.sort(key=lambda r: r['idx'])
inp = [json.loads(l)['d6'] for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
assert [r['d6'] for r in rows] == inp and [r['idx'] for r in rows] == list(range(len(inp))), 'coverage mismatch'
assert len(set(inp)) == len(inp), 'duplicate d6 in input'
S = {}
S['n_graphs'] = len(rows)
S['verdicts'] = dict(Counter(r['verdict'] for r in rows))
S['verdict_first_legal'] = dict(Counter(r.get('verdict_first_legal') for r in rows))
S['first_legal_idx_t'] = dict(Counter(str([100, 1000, 10000, 100000][r['first_legal_idx']]) for r in rows if r.get('first_legal_idx') is not None))
S['n_errors'] = sum(1 for r in rows if r['verdict'] == 'ERROR')
S['errors'] = [(r['d6'], r.get('error')) for r in rows if r['verdict'] == 'ERROR']
S['integrity_ok_all'] = all(r.get('integrity_ok') for r in rows)
S['integrity_failures'] = [r['d6'] for r in rows if not r.get('integrity_ok')]
S['n_complete'] = sum(1 for r in rows if r.get('complete'))
S['n_loops'] = sum(1 for r in rows if r.get('loops'))
S['sturm_count_hi'] = dict(Counter(str(r.get('sturm_count_hi')) for r in rows))
S['hi_minpoly_degree'] = dict(Counter(r.get('hi_minpoly_degree') for r in rows))
fl = [r['first_legal'] for r in rows if r.get('first_legal')]
S['first_legal_points'] = len(fl)
S['first_legal_agree_with_stored'] = sum(1 for b in fl if b['agrees_with_stored'])
S['first_legal_disagree'] = [(r['d6'], r['first_legal']['stable'], r['first_legal']['positive'], r['first_legal']['stored']) for r in rows if r.get('first_legal') and not r['first_legal']['agrees_with_stored']]
S['first_legal_crosschecks_all_ok'] = all(b['closed_form_ok'] and b['charpoly_crosscheck'] and b['routh_agrees'] is True for b in fl)
S['first_legal_undecided'] = sum(1 for b in fl if b['stable'] == 'UNDECIDED' or b['positive'] is None)
S['first_legal_stable_positive'] = sum(1 for b in fl if b['stable'] is True and b['positive'] is True)
# every-10th sample: all four witnesses
samp10 = [r for r in rows if r['idx'] % 10 == 0]
S['sample10_graphs'] = len(samp10)
S['sample10_points'] = sum(len(r['all_witnesses']) for r in samp10)
S['sample10_all_agree_with_stored'] = sum(1 for r in samp10 if r['all_witnesses_agree_with_stored'])
S['sample10_disagreements'] = [(r['d6'], [(w['t'], w['stable'], w['positive'], w['stored']) for w in r['all_witnesses'] if not w['agrees_with_stored']]) for r in samp10 if not r['all_witnesses_agree_with_stored']]
S['sample10_undecided_points'] = sum(1 for r in samp10 for w in r['all_witnesses'] if w['stable'] == 'UNDECIDED')
S['sample10_routh_all_agree'] = all(w['routh_agrees'] is True for r in samp10 for w in r['all_witnesses'] if w['stable'] != 'UNDECIDED')
# all graphs where all four were checked (sample + non-first-confirmed)
allw_rows = [r for r in rows if 'all_witnesses' in r]
S['all_four_checked_graphs'] = len(allw_rows)
S['all_four_checked_points'] = sum(len(r['all_witnesses']) for r in allw_rows)
S['all_four_checked_agree'] = sum(1 for r in allw_rows if r['all_witnesses_agree_with_stored'])
# every-50th sample: sharpness
samp50 = [r for r in rows if r['idx'] % 50 == 0]
S['sample50_graphs'] = len(samp50)
S['sample50_sharp'] = sum(1 for r in samp50 if r.get('sharpness', {}).get('sharp'))
S['sample50_not_sharp'] = [(r['d6'], r.get('sharpness')) for r in samp50 if not r.get('sharpness', {}).get('sharp')]
S['sample50_beyond_retries'] = sum(1 for r in samp50 if 'beyond_retries' in r.get('sharpness', {}))
S['sample50_inside_retries'] = sum(1 for r in samp50 if 'inside_retries' in r.get('sharpness', {}))
S['sample50_sturm_count_1'] = sum(1 for r in samp50 if r.get('sturm_count_hi') == 1)
S['sample50_d_in_values'] = dict(Counter(r['sharpness']['d_in'] for r in samp50 if 'sharpness' in r))
# own witnesses
own = [r for r in rows if 'own_witness' in r]
S['own_witness_graphs'] = len(own)
S['own_witness_confirmed'] = sum(1 for r in own if r['own_witness_confirms'])
S['own_witness_t_used'] = dict(Counter(r['own_witness'][-1]['t'] for r in own if r['own_witness_confirms']))
S['own_witness_failed'] = [(r['d6'], r['own_witness']) for r in own if not r['own_witness_confirms']]
S['later_stored_witness_graphs'] = sum(1 for r in rows if r['verdict'] == 'CONFIRMED_AT_LATER_STORED_WITNESS')
# points checked in total
def npts(r):
    k = 1 if r.get('first_legal') else 0
    k += len(r.get('all_witnesses', []))
    sh = r.get('sharpness', {})
    if sh:
        k += 2 + len(sh.get('beyond_retries', [])) + len(sh.get('inside_retries', []))
    k += len(r.get('own_witness', []))
    return k
S['points_checked_total'] = sum(npts(r) for r in rows)
S['undecided_graphs'] = [r['d6'] for r in rows if r['verdict'] == 'UNDECIDED']
S['not_confirmed_graphs'] = [r['d6'] for r in rows if r['verdict'] not in ('CONFIRMED_AT_STORED_FIRST_LEGAL_WITNESS', 'CONFIRMED_AT_LATER_STORED_WITNESS', 'CONFIRMED_AT_OWN_WITNESS')]
S['seconds_sum'] = round(sum(r['seconds'] for r in rows), 1)
S['seconds_max'] = max(r['seconds'] for r in rows)
json.dump({'summary': S, 'rows': rows}, open(f'{D}/v_sympy_n8_results.json', 'w'))
with open(f'{D}/v_sympy_n8_table.jsonl', 'w') as f:
    for r in rows:
        b = r.get('first_legal') or {}
        f.write(json.dumps({'idx': r['idx'], 'd6': r['d6'], 'arcs': r.get('arcs'), 'verdict': r['verdict'],
                            'first_legal_t': b.get('t'), 'first_legal_q': b.get('q'), 'first_legal_stable': b.get('stable'),
                            'first_legal_positive': b.get('positive'), 'first_legal_agrees_with_stored': b.get('agrees_with_stored'),
                            'sturm_count_hi': r.get('sturm_count_hi'), 'hi_minpoly_degree': r.get('hi_minpoly_degree'),
                            'all_witnesses_agree_with_stored': r.get('all_witnesses_agree_with_stored'),
                            'sharp': (r.get('sharpness') or {}).get('sharp'), 'hi_root_refined': r.get('hi_root_refined'),
                            'own_witness': [(o['beta'], o['t'], o['stable'], o['positive']) for o in r.get('own_witness', [])] or None,
                            'integrity_ok': r.get('integrity_ok'), 'seconds': r['seconds']}) + '\n')
print(json.dumps(S, indent=1, default=str)[:6000])
