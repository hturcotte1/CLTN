import sys, json, time
sys.path.insert(0, '.')
from fractions import Fraction as F
from vroots import d6_decode, is_complete_bidirected, check_point, canon
ROOT = '/home/user/CLTN/certificates/'
fin = [json.loads(l) for l in open(ROOT + 'COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
lim = {r['d6']: r for r in (json.loads(l) for l in open(ROOT + 'COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl'))}
fnfl = json.load(open(ROOT + 'COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))['graph']
out = {'decoder_check': {}, 'rows': [], 'sanity': []}

# decoder check
G = d6_decode('&EZjtwWE')
gen_file = open('/home/user/CLTN/n7/geneson.d6').read().strip()
out['decoder_check'] = {'d6': '&EZjtwWE', 'n': len(G), 'arcs': sum(map(sum, G)),
    'rows': [''.join(map(str, r)) for r in G], 'n7/geneson.d6': gen_file}
gp = check_point(G, F(1, 100), F(29, 2500))
out['sanity'].append({'what': 'Geneson 6-node graph at eps=1/100 delta=29/2500 (expect counterexample)',
    'clique': is_complete_bidirected(G), **{k: gp[k] for k in ('legal', 'positive', 'stable', 'rhp_routh', 'rhp_cauchy', 'rhp_sympy', 'det_IminusW')}})
Afn = d6_decode('&FNFL~[udo?')
out['decoder_check']['FNFL_rows_match_certificate'] = [''.join(map(str, r)) for r in Afn] == fnfl['adjacency_rows']
print(json.dumps(out['decoder_check'], indent=1)); print(out['sanity'][-1]); sys.stdout.flush()

def tdpt(beta, t):
    eps = 1 / (t + 1); delta = (beta - 1) * eps
    return eps, delta

def summarize(tag, p):
    keep = ('eps', 'delta', 'q', 'legal', 'det_IminusW', 'nondegenerate', 'positive', 'stable',
            'rhp_routh', 'rhp_cauchy', 'rhp_sympy', 'methods_agree')
    d = {'point': tag}; d.update({k: p.get(k) for k in keep})
    for k in ('routh_undecided_first_col', 'cauchy_note', 'sympy_note'):
        if k in p: d[k] = p[k]
    d['counterexample_point'] = bool(p['legal'] and p['positive'] is True and p['stable'] is True)
    return d

canon_seen = {}
t0 = time.time()
for idx, r in enumerate(fin):
    A = d6_decode(r['d6'])
    row = {'d6': r['d6'], 'edges_decoded': sum(map(sum, A)), 'edges_stored': r['edges'],
           'clique': is_complete_bidirected(A), 'clique_stored': r['clique'], 'points': [], 'issues': []}
    if row['edges_decoded'] != r['edges']: row['issues'].append('edge count mismatch')
    if row['clique'] != r['clique']: row['issues'].append('clique flag mismatch')
    ck = canon(A)
    if ck in canon_seen: row['issues'].append('isomorphic to ' + canon_seen[ck])
    canon_seen[ck] = r['d6']
    pts = []
    if r['d6'] in lim:
        L = lim[r['d6']]
        if [list(x) for x in L['A']] != A: row['issues'].append('limit-file A differs from decoded d6')
        for w in L['pointwise_witnesses']:
            eps, delta = F(w['eps']), F(w['delta'])
            e2, d2 = tdpt(F(w['beta']), F(w['t']))
            if (e2, d2) != (eps, delta): row['issues'].append('witness beta/t inconsistent with eps/delta')
            pts.append(('witness', eps, delta, w))
    elif r['d6'] == '&FNFL~[udo?':
        for w in fnfl['witnesses']['stable_positive_legal']:
            pts.append(('witness', F(w['eps']), F(w['delta']), w))
    else:
        row['issues'].append('no stored pointwise witness found')
    c0 = r['stable_cells'][0]
    blo, bhi, bs, tlo, thi = c0
    bs = F(bs); tlo = F(tlo)
    # stored CAD sample point of first stable cell
    samp = r['region2d']['cells']
    ts = None
    for c in samp:
        if c['beta_cell'] == [blo, bhi]:
            for tc in c['t_cells']:
                if tc['t_lo'] == c0[3] and tc['t_hi'] == c0[4]:
                    ts = F(tc['t_sample'])
    if ts is None: row['issues'].append('could not locate stored t_sample of first stable cell')
    t2 = tlo + (F(thi) - tlo) / 3 if thi is not None else 3 * tlo + 5
    if ts is not None and t2 == ts: t2 = tlo + (F(thi) - tlo) / 5 if thi is not None else 4 * tlo + 7
    e, d = tdpt(bs, t2); pts.append(('cell1_second_point beta=%s t=%s' % (bs, t2), e, d, None))
    if ts is not None:
        e, d = tdpt(bs, ts); pts.append(('cell1_stored_sample beta=%s t=%s' % (bs, ts), e, d, None))
    for tag, eps, delta, w in pts:
        p = check_point(A, eps, delta)
        s = summarize(tag, p)
        if w is not None and 'x' in w and 'x' in p and w['x'] != p['x']:
            s['x_mismatch_with_stored'] = True; row['issues'].append('x mismatch at ' + tag)
        if w is not None and 'charpoly' in w and w['charpoly'] != p['charpoly']:
            row['issues'].append('charpoly mismatch at ' + tag)
        if w is not None:
            stored = (w.get('legal'), w.get('is_fixed_point', w.get('positive')), w.get('stable'))
            mine = (p['legal'], p['positive'], p['stable'])
            s['stored_flags_legal_positive_stable'] = list(stored)
            if stored != mine:
                row['issues'].append('stored flags %s != recomputed %s at %s eps=%s' % (stored, mine, tag, eps))
        elif not s['counterexample_point']:
            row['issues'].append('point does not check: ' + tag)
        if not p['nondegenerate']:
            row['issues'].append('det(I-W)=0 at ' + tag)
        if not p['methods_agree']:
            row['issues'].append('methods disagree at ' + tag)
        row['points'].append(s)
    wit = [s for s in row['points'] if s['point'] == 'witness']
    row['all_witnesses_check'] = bool(wit) and all(s['counterexample_point'] for s in wit)
    row['n_witnesses'] = len(wit)
    row['n_witnesses_counterexample'] = sum(s['counterexample_point'] for s in wit)
    row['n_witnesses_stored_illegal'] = sum(1 for s in wit if s['stored_flags_legal_positive_stable'][0] is False)
    row['second_point_checks'] = row['points'][len(wit)]['counterexample_point']
    row['verdict'] = ('CONFIRMED' if (not row['clique'] and row['n_witnesses_counterexample'] >= 1 and row['second_point_checks']
                       and not row['issues']) else 'NOT CONFIRMED')
    out['rows'].append(row)
    print(idx, r['d6'], row['verdict'], row['issues'], '%.0fs' % (time.time() - t0)); sys.stdout.flush()

# sanity: FNFL stored unstable points must come out unstable
for w in fnfl['witnesses']['unstable_same_beta_larger_t']:
    p = check_point(Afn, F(w['eps']), F(w['delta']))
    out['sanity'].append({'what': 'FNFL stored UNSTABLE point eps=%s delta=%s (expect stable False)' % (w['eps'], w['delta']),
                          'legal': p['legal'], 'positive': p['positive'], 'stable': p['stable'],
                          'rhp': [p['rhp_routh'], p['rhp_cauchy'], p['rhp_sympy']]})
rows = out['rows']
out['summary'] = {'n_graphs': len(rows), 'n_confirmed': sum(r['verdict'] == 'CONFIRMED' for r in rows),
    'n_points_checked': sum(len(r['points']) for r in rows),
    'n_points_counterexample': sum(s['counterexample_point'] for r in rows for s in r['points']),
    'n_undecided_points': sum(s['stable'] is None for r in rows for s in r['points']),
    'n_method_disagreements': sum(not s['methods_agree'] for r in rows for s in r['points']),
    'n_degenerate_points': sum(not s['nondegenerate'] for r in rows for s in r['points']),
    'n_distinct_isomorphism_classes': len(canon_seen),
    'n_graphs_all_stored_witnesses_counterexample': sum(r['all_witnesses_check'] for r in rows),
    'n_stored_witnesses_flagged_illegal_in_file': sum(r['n_witnesses_stored_illegal'] for r in rows),
    'n_points_with_repeated_roots': sum(1 for r in rows for s in r['points'] if 'repeated' in str(s.get('sympy_note'))),
    'failures': [(r['d6'], r['issues']) for r in rows if r['verdict'] != 'CONFIRMED']}
json.dump(out, open('v_roots_results.json', 'w'), indent=1)
print(json.dumps(out['summary'], indent=1)); print(json.dumps(out['sanity'], indent=1))
