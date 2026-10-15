"""Join exact R* (rstar_raw.json) with the two-parameter certificates and classify agreement.
Exact comparisons only (Fractions); decimals in the output are labelled display-only."""
import json, sys, itertools
from fractions import Fraction as F
D = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf5/bulge/'
L = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
R = {r['d6']: r for r in json.load(open(D + 'rstar_raw.json'))}
H = json.load(open('/home/user/CLTN/notes/prover_scripts_M/n7_first_order_table.json'))
PLASTIC = ["&F^UXk\\Esb?", "&F^U|eNUTw?", "&F^VLOSekx?", "&F^u\\eNVoJ?", "&F^uk\\^_TU?", "&F^vlKUe[R?"]

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A
def canon(A):
    n = len(A); best = None
    for P in itertools.permutations(range(n)):
        key = tuple(A[P[i]][P[j]] for i in range(n) for j in range(n))
        if best is None or key > best: best = key
    return best
cmap = {}
for r in L: cmap[canon(from_d6(r['d6']))] = r['d6']
Hmap = {}
for h in H: Hmap[cmap.get(canon(from_d6(h['d6'])))] = h
plast = {cmap.get(canon(from_d6(g))) for g in PLASTIC}

rows = []; dis = []
for r in L:
    q = R[r['d6']]
    bc = [F(q['beta_c']['interval'][0]), F(q['beta_c']['interval'][1])] if q['beta_c']['interval'][0] != 'None' else None
    bounded = [c for c in r['stable_cells'] if c[4] is not None]
    unb = [c for c in r['stable_cells'] if c[4] is None]
    loc = []
    for c in bounded:
        lo, hi = F(c[0]), F(c[1])
        if bc and lo >= bc[0] - F(1, 10**3) and lo <= bc[1] + F(1, 10**2):
            loc.append('adjacent right of beta_c')
        elif bc and lo >= bc[1]:
            loc.append('right of beta_c, not adjacent (beta_lo - beta_c ~ %.4f)' % float(lo - bc[1]))
        else:
            loc.append('left of beta_c (inside limit-stable range)')
    unb_right = [c for c in unb if bc and F(c[0]) >= bc[1]]
    sign = q.get('Rstar', {}).get('sign') if q['status'] == 'ok' else None
    nu_half = q['nu_vs_half']
    # first-order prediction of a bulge cell at (beta_c, inf)
    if q['status'] != 'ok':
        pred = None; why = q['status']
    elif nu_half == 0:
        pred = (q['Rstar_vs_1/8'] == 1) and q['positive_at_beta_c']; why = 'beta_c = 2: compare R* with 1/8 (Lemma B)'
    else:
        pred = (sign == 1) and q['positive_at_beta_c']; why = 'Theorem AS(a)/(b)'
    row = {'d6': r['d6'], 'alpha_minpoly_hi_to_lo': q['alpha']['minpoly_hi_to_lo'], 'alpha_interval': q['alpha']['interval'],
           'nu_inf_interval': q['nu_inf']['interval'],
           'beta_c_minpoly_hi_to_lo': q['beta_c']['minpoly_hi_to_lo'], 'beta_c_interval': q['beta_c']['interval'],
           'z0': [{'factor_hi_lo': c['factor_hi_lo'], 'real': c['real'], 'kind': c['kind'],
                   'isolating_rect_re_im': c['rect'], 'mult_in_psi': c['mult_psi'], 'mult_in_chi': c['mult_chi']} for c in q['critical_roots']],
           'tie_between_nonconjugate_roots': q['tie_between_distinct_conjugate_classes'],
           'gcd_chi_psi_hi_lo': q['gcd_chi_psi_hi_lo'],
           'Rstar_minpoly_hi_to_lo': q.get('Rstar', {}).get('minpoly_hi_to_lo'),
           'Rstar_interval': q.get('Rstar', {}).get('interval'),
           'Rstar_exact_rational': q.get('Rstar', {}).get('exact'),
           'Rstar_sign': sign, 'Rstar_vs_1_8': q.get('Rstar_vs_1/8'),
           'beta_c3_Rstar_display_only': q.get('beta_c3_Rstar_interval_rounded', [None])[0],
           'positive_at_beta_c_large_t': q['positive_at_beta_c'], 'positivity_signs_at_beta_c': q['positivity_signs_at_beta_c'],
           'bounded_in_t_somewhere': r['bounded_in_t_somewhere'],
           'bounded_cells': [[c[0], c[1], c[3], c[4]] for c in bounded], 'bounded_cell_location': loc,
           'unbounded_stable_cells_right_of_beta_c': len(unb_right),
           'limit_intervals_lo_kinds': [li['lo']['kind'] for li in (r['limit_intervals'] or [])],
           'limit_intervals_hi_kinds': [li['hi'].get('kind') for li in (r['limit_intervals'] or [])],
           'limit_intervals_hi_vanishing': [li['hi'].get('vanishing') for li in (r['limit_intervals'] or [])],
           'plastic_graph': r['d6'] in plast,
           'first_order_prediction_bulge_at_beta_c': pred, 'prediction_rule': why,
           'Rstar_pos': sign == 1,
           'agree_Rpos_vs_bounded': (sign == 1) == r['bounded_in_t_somewhere'],
           'agree_prediction_vs_bounded': pred == r['bounded_in_t_somewhere'],
           'heuristic_float_R_display_only': Hmap.get(r['d6'], {}).get('R'),
           'heuristic_prediction': Hmap.get(r['d6'], {}).get('prediction')}
    if Hmap.get(r['d6']) and Hmap[r['d6']].get('R') is not None and sign is not None:
        hs = 1 if Hmap[r['d6']]['R'] > 0 else -1
        row['heuristic_sign_agrees'] = hs == sign
    rows.append(row)
json.dump(rows, open(D + 'rstar_table.json', 'w'), indent=1)
from collections import Counter
print('n', len(rows))
print('status', Counter(R[r['d6']]['status'] for r in L))
print('sign', Counter(x['Rstar_sign'] for x in rows))
print('Rpos x bounded', Counter((x['Rstar_sign'], x['bounded_in_t_somewhere']) for x in rows))
print('pred x bounded', Counter((x['first_order_prediction_bulge_at_beta_c'], x['bounded_in_t_somewhere']) for x in rows))
print('pos at beta_c', Counter((x['Rstar_sign'], x['positive_at_beta_c_large_t']) for x in rows))
print('loc', Counter(l for x in rows for l in x['bounded_cell_location']))
print('heur agree', Counter(x.get('heuristic_sign_agrees') for x in rows))
print('unb right', sum(x['unbounded_stable_cells_right_of_beta_c'] for x in rows))
print('ties', [x['d6'] for x in rows if x['tie_between_nonconjugate_roots']])
for x in rows:
    if not x['agree_Rpos_vs_bounded'] or not x['agree_prediction_vs_bounded']:
        print('DIS', x['d6'], 'R*sign', x['Rstar_sign'], 'pos', x['positive_at_beta_c_large_t'], 'bdd', x['bounded_in_t_somewhere'],
              x['bounded_cell_location'], 'plastic', x['plastic_graph'], x['limit_intervals_lo_kinds'], x['limit_intervals_hi_kinds'], x['limit_intervals_hi_vanishing'], 'bc', [float(F(v)) for v in x['beta_c_interval']])
