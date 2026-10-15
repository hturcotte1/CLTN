import json, sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vf import *

REPO = '/home/user/CLTN'
OUT = os.path.dirname(os.path.abspath(__file__))
fin = [json.loads(l) for l in open(f'{REPO}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
lim = {json.loads(l)['d6']: json.loads(l) for l in open(f'{REPO}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')}
fnfl = json.load(open(f'{REPO}/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
bulge = json.load(open(f'{REPO}/certificates/n7_limit_hits_finite_eps_bulge_witnesses.json'))
S = lambda f: str(f)
failures = []
log = {}

# ---------- decoder checks ----------
G = decode_d6('&EZjtwWE')
gen = open(f'{REPO}/n7/geneson.d6').read().split()[0]
dec = dict(geneson_string=gen, n=len(G), arcs=sum(map(sum, G)), rows=[''.join(map(str, r)) for r in G],
           node0_out=[j for j in range(6) if G[0][j]], loops=has_loops(G), roundtrip=encode_d6(G) == '&EZjtwWE')
discrepancies = []
if not (dec['n'] == 6 and dec['arcs'] == 19 and dec['roundtrip'] and gen == '&EZjtwWE' and not dec['loops']):
    failures.append(f'decoder check on Geneson failed: {dec}')
if dec['node0_out'] != [1, 2, 3, 4]:
    # bits of '&EZjtwWE': Z=27=011011 j=43=101011 t=53=110101 w=56=111000 W=24=011000 E=6=000110
    discrepancies.append('task text says Geneson node 0 -> 1,2,3,4; the specified decoding gives node 0 -> %s '
                         '(row 0 = 011011 from Z=27). Column-major reading would give 0 -> {1,2,3}. Neither matches; the '
                         'same decoder reproduces all 123 stored A matrices, FNFL rows and both bulge rows exactly, so the '
                         'hint in the task text is taken to be wrong, not the decoder.' % dec['node0_out'])
Af = decode_d6('&FNFL~[udo?')
dec['FNFL_rows_match'] = [''.join(map(str, r)) for r in Af] == fnfl['graph']['adjacency_rows']
if not dec['FNFL_rows_match']: failures.append('FNFL rows mismatch')
for d6, b in bulge['graphs'].items():
    ok = [''.join(map(str, r)) for r in decode_d6(d6)] == b['adjacency_rows']
    dec[f'bulge_rows_match_{d6}'] = ok
    if not ok: failures.append(f'bulge rows mismatch {d6}')
mism = [d for d, r in lim.items() if decode_d6(d) != r['A']]
dec['limit_A_mismatches'] = mism
if mism: failures.append(f'limit A mismatches: {mism}')
# controls
g_res = analyze(G, F(1, 100), F(29, 2500))
dec['geneson_control_eps1/100_delta29/2500'] = dict(legal=g_res['legal'], positive=g_res['positive'], stability=g_res['stability'],
                                                    q=S(F(29, 2500) * F(99, 100) * 100))
if not (g_res['legal'] and g_res['positive'] and g_res['stability'] == 'stable'):
    failures.append(f'Geneson positive control did not flag: {dec["geneson_control_eps1/100_delta29/2500"]}')
K7 = [[int(i != j) for j in range(7)] for i in range(7)]
k_res = analyze(K7, F(1, 4), F(1, 2))
dec['K7_control'] = dict(positive=k_res['positive'], stability=k_res['stability'])
if not (k_res['positive'] and k_res['stability'] == 'stable'): failures.append('K7 control failed')
C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
c_res = analyze(C3, F(1, 4), F(1, 2))
dec['C3_control_expect_unstable'] = dict(positive=c_res['positive'], stability=c_res['stability'])
if c_res['stability'] != 'unstable': failures.append('3-cycle control not unstable')
log['decoder_and_controls'] = dec

def ed_from_bt(beta, t):
    eps = 1 / (t + 1)
    return eps, (beta - 1) * eps

def point(A, eps, delta, label, stored=None):
    r = analyze(A, eps, delta)
    q = delta * (1 - eps) / eps
    row = dict(label=label, eps=S(eps), delta=S(delta), q=S(q), legal=r['legal'], positive=r['positive'],
               stable=r['stability'], minor_signs=[(d > 0) - (d < 0) for d in r['minors']],
               min_x=S(min(r['x'])) if r['x'] else None)
    if stored is not None:
        if 'x' in stored:
            row['x_matches_stored'] = [F(v) for v in stored['x']] == r['x']
        if 'charpoly' in stored:
            row['charpoly_matches_stored'] = [F(v) for v in stored['charpoly']] == r['charpoly']
    return row, r

table = []
illegal_first = []
for rec in fin:
    d6 = rec['d6']
    A = decode_d6(d6)
    clique = is_complete_bidirected(A)
    ent = dict(d6=d6, edges=sum(map(sum, A)), clique=clique, loops=has_loops(A), points=[])
    if ent['edges'] != rec['edges']: failures.append(f'{d6}: edge count {ent["edges"]} != stored {rec["edges"]}')
    if clique != rec['clique']: failures.append(f'{d6}: clique flag mismatch')
    # witness point
    if d6 in lim:
        w = lim[d6]['pointwise_witnesses'][0]
        eps, delta = F(w['eps']), F(w['delta'])
        assert (F(w['beta']), F(w['t'])) and ed_from_bt(F(w['beta']), F(w['t'])) == (eps, delta)
        row, _ = point(A, eps, delta, 'stored_witness_0', w)
        row['stored_legal_flag'] = w['legal']
        if not row['legal']:
            # the first stored witness is outside the legal region (the file itself says legal=false);
            # record it as an informational row and also check every stored witness.
            row['note'] = 'first stored witness is ILLEGAL (q<=1); stored flag legal=%s' % w['legal']
            ent['illegal_first_witness'] = row
            ent['all_stored_witnesses'] = []
            for k, ww in enumerate(lim[d6]['pointwise_witnesses']):
                rr, _ = point(A, F(ww['eps']), F(ww['delta']), f'stored_witness_{k}', ww)
                rr['stored_legal_flag'] = ww['legal']; ent['all_stored_witnesses'].append(rr)
            legal_ws = [rr for rr in ent['all_stored_witnesses'] if rr['legal']]
            if not legal_ws:
                failures.append(f'{d6}: no legal stored witness'); row = dict(row, label='NO_LEGAL_WITNESS')
            else:
                row = legal_ws[0]; row['label'] = row['label'] + '_first_legal'
    elif d6 == '&FNFL~[udo?':
        w = [x for x in fnfl['graph']['witnesses']['stable_positive_legal'] if x['eps'] == '1/36'][0]
        row, _ = point(A, F(1, 36), F(17, 588), 'stored_witness_eps1/36', w)
    else:
        failures.append(f'{d6}: no witness found'); continue
    ent['points'].append(row)
    # point inside first stable cell, and in every stable cell
    for ci, cell in enumerate(rec['stable_cells']):
        blo, bhi, bs, tlo, thi = cell
        beta = F(bs); tl = F(tlo)
        t = tl + 1 if thi is None else (tl + F(thi)) / 2
        inside = F(blo) < beta < F(bhi) and t > tl and (thi is None or t < F(thi))
        eps, delta = ed_from_bt(beta, t)
        row, _ = point(A, eps, delta, f'cell{ci}' + ('_first' if ci == 0 else ''))
        row['beta'] = S(beta); row['t'] = S(t); row['inside_cell'] = inside
        ent['points'].append(row)
    for row in ent['points']:
        ok = row['legal'] and row['positive'] and row['stable'] == 'stable' and not clique
        row['verdict'] = 'COUNTEREXAMPLE_CONFIRMED' if ok else 'NOT_CONFIRMED'
        if not ok: failures.append(f"{d6} {row['label']}: legal={row['legal']} positive={row['positive']} stable={row['stable']} clique={clique}")
        for k in ('x_matches_stored', 'charpoly_matches_stored', 'inside_cell'):
            if row.get(k) is False: failures.append(f"{d6} {row['label']}: {k} False")
    if 'illegal_first_witness' in ent:
        r0 = ent['illegal_first_witness']
        if r0['legal'] != r0['stored_legal_flag']: failures.append(f'{d6}: my legality disagrees with stored flag')
        illegal_first.append(dict(d6=d6, eps=r0['eps'], delta=r0['delta'], q=r0['q'], stable_there=r0['stable'],
                                  positive_there=r0['positive'], stored_legal_flag=r0['stored_legal_flag'],
                                  substitute=ent['points'][0]['label'], substitute_eps=ent['points'][0]['eps']))
        for rr in ent['all_stored_witnesses']:
            if rr['legal'] != rr['stored_legal_flag']: failures.append(f'{d6} {rr["label"]}: legality disagrees with stored flag')
            if rr.get('x_matches_stored') is False or rr.get('charpoly_matches_stored') is False:
                failures.append(f'{d6} {rr["label"]}: stored x/charpoly mismatch')
    ent['verdict'] = 'COUNTEREXAMPLE_CONFIRMED' if all(r['verdict'] == 'COUNTEREXAMPLE_CONFIRMED' for r in ent['points'][:2]) else 'NOT_CONFIRMED'
    ent['all_cells_confirmed'] = all(r['verdict'] == 'COUNTEREXAMPLE_CONFIRMED' for r in ent['points'])
    table.append(ent)

# ---------- unstable points ----------
unst = []
def chk_unstable(d6, w):
    A = decode_d6(d6)
    row, r = point(A, F(w['eps']), F(w['delta']), 'stored_unstable', w)
    neg = [k + 1 for k, d in enumerate(r['minors']) if d < 0]
    zero = [k + 1 for k, d in enumerate(r['minors']) if d == 0]
    row.update(d6=d6, negative_minors=neg, zero_minors=zero,
               decided_unstable_strict=(r['stability'] == 'unstable' and len(neg) > 0),
               stored_detail=w.get('stable_detail'))
    if not row['decided_unstable_strict']: failures.append(f'{d6} eps={w["eps"]}: not strictly decided unstable ({row})')
    for k in ('x_matches_stored', 'charpoly_matches_stored'):
        if row.get(k) is False: failures.append(f"{d6} unstable eps={w['eps']}: {k} False")
    unst.append(row)
for w in fnfl['graph']['witnesses']['unstable_same_beta_larger_t']:
    chk_unstable('&FNFL~[udo?', w)
for d6, b in bulge['graphs'].items():
    for w in b['witnesses']['unstable_same_beta_larger_t']:
        chk_unstable(d6, w)
    for w in b['witnesses']['stable_positive_legal_beyond_limit_interval']:
        row, r = point(decode_d6(d6), F(w['eps']), F(w['delta']), 'stored_bulge_stable', w)
        row['d6'] = d6; unst.append(row)
        if not (row['legal'] and row['positive'] and row['stable'] == 'stable'): failures.append(f'{d6} bulge stable witness failed')
# also FNFL's second stored stable witness eps=1/49
w = [x for x in fnfl['graph']['witnesses']['stable_positive_legal'] if x['eps'] == '1/49'][0]
row, _ = point(decode_d6('&FNFL~[udo?'), F(w['eps']), F(w['delta']), 'stored_witness_eps1/49', w)
row['d6'] = '&FNFL~[udo?'; unst.append(row)
if not (row['legal'] and row['positive'] and row['stable'] == 'stable'): failures.append('FNFL eps=1/49 witness failed')

trans = {}
for ent in table:
    A = decode_d6(ent['d6']); AT = [list(r) for r in zip(*A)]
    p0 = ent['points'][0]; r = analyze(AT, F(p0['eps']), F(p0['delta']))
    trans[ent['d6']] = (r['positive'], r['stability'])
n_trans_cex = sum(1 for v in trans.values() if v[0] and v[1] == 'stable')
summary = dict(transposed_convention_still_cex_at_primary_point=n_trans_cex, n_graphs=len(table), n_confirmed=sum(e['verdict'] == 'COUNTEREXAMPLE_CONFIRMED' for e in table),
               n_all_cells_confirmed=sum(e['all_cells_confirmed'] for e in table),
               n_points_checked=sum(len(e['points']) for e in table), distinct_d6=len({e['d6'] for e in table}),
               n_failures=len(failures), n_first_witness_illegal=len(illegal_first))
json.dump(dict(summary=summary, failures=failures, discrepancies=discrepancies, first_stored_witness_illegal=illegal_first, decoder_and_controls=log['decoder_and_controls'],
               table=table, unstable_and_extra_points=unst), open(f'{OUT}/v_fraction_results.json', 'w'), indent=1)
print(json.dumps(summary, indent=1)); print('FAILURES:', failures[:50]); print('DISCREPANCIES:', discrepancies); print(json.dumps(illegal_first, indent=0)); print(json.dumps(log, indent=1))
for r in unst: print(r['d6'], r['label'], r['eps'], r['stable'], r.get('negative_minors'), r.get('zero_minors'), r['minor_signs'], r.get('x_matches_stored'), r.get('charpoly_matches_stored'))
