import json, collections, sys
from fractions import Fraction as F
OUT = '/home/user/CLTN/notes/verify_n8/v_structure/'
D = json.load(open(OUT + 'verify_out.json')); rows = D['rows']
N = len(rows)
c = collections.Counter()
table = []; legal = []; notes = []
for r in rows:
    iv = r['intervals'][0]
    ok = (r['n'] == 8 and not r['complete_bidirected'] and r['domination_free'] and r['recommended_point'] is not None
          and r['recommended_point']['stable'] == 'stable' and r['recommended_point']['positive'] and F(r['recommended_point']['q']) > 1)
    c['confirmed'] += ok
    c['not_complete'] += not r['complete_bidirected']; c['domfree'] += r['domination_free']
    c['first_illegal'] += r['first_witness_illegal']
    c['stored_legal'] += r['stored_legal_witness'] is not None
    c['all4'] += all(w['stable'] == 'stable' and w['positive'] for w in r['witnesses'])
    c['wit_conf'] += r['n_witnesses_confirmed']; c['wit_legal'] += sum(w['legal'] for w in r['witnesses'])
    c['wit_stable_pos'] += sum(w['stable'] == 'stable' and w['positive'] for w in r['witnesses'])
    c['wit_undecided'] += sum(w['stable'] == 'undecided' for w in r['witnesses'])
    c['wit_fields'] += all(w['fields_consistent'] for w in r['witnesses'])
    c['wit_inside'] += r['n_witnesses_outside_interval'] == 0
    c['lo_root'] += iv['lo_kind'] == 'root'
    c['iso_ok'] += iv['hi_isolating_ok'] and iv.get('lo_isolating_ok', True)
    c['rec_tried'] += len(r['recommended_tried']); c['rec_undecided'] += sum(x[1] == 'undecided' for x in r['recommended_tried'])
    c['points'] += len(r['witnesses']) + len(r['recommended_tried'])
    row = {'d6': r['d6'], 'arcs': r['arcs'], 'missing_arcs': r['missing_arcs'], 'complete_bidirected': r['complete_bidirected'],
           'domination_free': r['domination_free'], 'lo': iv['lo_kind'], 'hi_beta_interval': iv['hi_beta_interval'], 'hi_deg': iv['hi_deg'],
           'q_upper_bound': iv['hi_q_upper_bound'], 'witnesses': [(w['t'], w['legal'], w['stable'], w['positive']) for w in r['witnesses']],
           'first_witness_illegal': r['first_witness_illegal'], 'recommended': r['recommended_point'],
           'verdict': 'CONFIRMED counterexample' if ok else 'NOT CONFIRMED'}
    table.append(row)
    legal.append({'d6': r['d6'],
                  'witness_stored_legal': ({k: r['stored_legal_witness'][k] for k in ('eps', 'delta', 'q')} if r['stored_legal_witness'] else None),
                  'witness_in_first_stable_cell': {k: r['recommended_point'][k] for k in ('eps', 'delta', 'beta', 't', 'q')},
                  'beta_certified_by': r['recommended_point']['beta_certified_by']})
    if r['n_witnesses_outside_interval']:
        notes.append({'d6': r['d6'], 'stored_witness_beta': r['witnesses'][0]['beta'], 'hi_beta_interval': iv['hi_beta_interval'],
                      'refined_hi': r['refined_intervals']['hi'], 'stored_stable_flags': [w['stored_stable'] for w in r['witnesses']],
                      'my_verdicts': [w['stable'] for w in r['witnesses']], 'legal_flags': [w['legal'] for w in r['witnesses']]})
# q upper ends
qs = sorted(((F(r['intervals'][0]['hi_q_upper_bound']), r['d6'], r['intervals'][0]['hi_beta_interval'], r['intervals'][0]['hi_deg']) for r in rows), reverse=True)
top10 = [{'d6': d, 'q_upper_bound_rational': str(q), 'q_upper_bound_float_approx': float(q), 'hi_beta_isolating_interval': bi, 'hi_minpoly_degree': dg} for q, d, bi, dg in qs[:10]]
qmin = qs[-1]
summary = {
 'verifier': 'v_structure (n = 8)', 'n_graphs': N, 'n_confirmed': c['confirmed'],
 'pairwise_non_isomorphic': {'own_individualisation_refinement_canonical_forms_distinct': D['own_canon_distinct'], 'duplicates': len(D['duplicates']), 'labelg': D['labelg']},
 'none_complete_bidirected': c['not_complete'] == N, 'missing_arcs_distribution': dict(sorted(collections.Counter(r['missing_arcs'] for r in rows).items())),
 'arcs_distribution': dict(sorted(collections.Counter(r['arcs'] for r in rows).items())),
 'domination_free_all': c['domfree'] == N, 'domination_free_count': c['domfree'],
 'edges_field_equals_arcs_all': not any('edges field' in p for d, p in D['problems']),
 'witness_fields_consistent_all': c['wit_fields'] == N,
 'witness_beta_strictly_inside_interval_graphs': c['wit_inside'], 'witness_beta_outside_interval_graphs': N - c['wit_inside'],
 'stored_witness_points_checked': 4 * N, 'stored_witness_points_stable_positive_mine': c['wit_stable_pos'], 'stored_witness_points_legal': c['wit_legal'],
 'stored_witness_points_legal_and_confirmed': c['wit_conf'], 'stored_witness_points_undecided': c['wit_undecided'],
 'stored_flags_agree_with_mine_all_points': not any(('stable:' in p or 'positive:' in p) for d, p in D['problems']),
 'graphs_all_four_stored_witnesses_stable_positive': c['all4'], 'graphs_with_stored_legal_witness': c['stored_legal'],
 'graphs_first_witness_illegal_q_le_1': c['first_illegal'],
 'recommended_points_checked': c['rec_tried'], 'recommended_points_undecided': c['rec_undecided'],
 'recommended_point_found_for_all': all(r['recommended_point'] for r in rows),
 'recommended_t_distribution': dict(collections.Counter(r['recommended_point']['t'] for r in rows)),
 'bisection_steps_distribution': dict(collections.Counter(r['refined_intervals']['bisection_steps'] for r in rows)),
 'isolating_intervals_certified_all': c['iso_ok'] == N,
 'stable_intervals_per_record': dict(collections.Counter(r['n_stable_intervals'] for r in rows)),
 'cells_inconsistent_with_intervals': sum(r['cells_inconsistent_with_intervals'] for r in rows), 'cells_total': sum(r['n_cells'] for r in rows),
 'lo_kind_root_count': c['lo_root'], 'lo_kind_beta_min_count': N - c['lo_root'],
 'hi_minpoly_degree_distribution': dict(sorted(collections.Counter(r['intervals'][0]['hi_deg'] for r in rows).items())),
 'lo_minpoly_degree_distribution': dict(sorted(collections.Counter(r['intervals'][0].get('lo_deg', 0) for r in rows).items())),
 'largest_q_upper_bound': {'rational': str(qs[0][0]), 'float_approx': float(qs[0][0]), 'd6': qs[0][1]},
 'smallest_q_upper_bound': {'rational': str(qmin[0]), 'float_approx': float(qmin[0]), 'd6': qmin[1]},
 'ten_largest_q_upper_bounds': top10,
 'q_upper_bound_below_2_all': all(q < 2 for q, *_ in qs),
 'total_exact_points_checked': c['points'],
 'notes_witness_beta_outside_interval': notes,
 'discrepancies': D['problems'], 'elapsed_s': D['elapsed_s']}
json.dump({'summary': summary, 'table': table}, open(OUT + 'v_structure_table.json', 'w'), indent=1)
json.dump({'purpose': 'one LEGAL exact witness per n = 8 candidate (q > 1) at which the full support {0..7} is a stable positive fixed point, re-decided by v_structure (independent exact implementation: Fractions, Faddeev-LeVerrier charpoly, Hurwitz minors), plus a recommended rational interior point for hand checks',
           'note': '%d of the 9,133 candidates have an ILLEGAL first stored witness (t = 100, q <= 1); %d candidates have NO stored legal stable witness (their stored beta = 2001/1000 lies above the stable interval; witness_stored_legal is null for them) -- for those use witness_in_first_stable_cell, whose beta is certified inside the stable interval by the isolating intervals (after exact bisection where needed) and whose stability/positivity was re-decided exactly at the stated t. eps = 1/(t+1), delta = (beta-1)*eps.' % (c['first_illegal'], N - c['stored_legal']),
           'graphs': legal}, open(OUT + 'n8_candidates_legal_witnesses.json', 'w'), indent=1)
# q table as text
with open(OUT + 'q_upper_ends.txt', 'w') as f:
    f.write('stable q-interval upper ends (rational upper bound b - 1 of the hi isolating interval [a, b] of beta)\n')
    f.write('largest: %s (~%.6f) at %s\n' % (str(qs[0][0]), float(qs[0][0]), qs[0][1]))
    f.write('ten largest:\n')
    for q, d, bi, dg in qs[:10]: f.write('  %-16s %-10.6f %s  beta in [%s, %s]  deg %d\n' % (str(q), float(q), d, bi[0], bi[1], dg))
    f.write('smallest: %s (~%.6f) at %s\n' % (str(qmin[0]), float(qmin[0]), qmin[1]))
    f.write('lo kind: root (algebraic beta > 2): %d; beta_min = 2: %d\n' % (c['lo_root'], N - c['lo_root']))
    hist = collections.Counter(min(int(float(q) * 10) / 10, 1.9) for q, *_ in qs)
    f.write('histogram of q upper bound (0.1 bins): %s\n' % dict(sorted(hist.items())))
print(json.dumps({k: v for k, v in summary.items() if k not in ('notes_witness_beta_outside_interval', 'ten_largest_q_upper_bounds', 'discrepancies')}, indent=1))
print(open(OUT + 'q_upper_ends.txt').read())
