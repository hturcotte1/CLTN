import json, collections
D = json.load(open('verify_out.json')); Fn = json.load(open('fnfl_out.json'))
rows = []
c = collections.Counter()
for r in D['table']:
    s0 = r['finite_cell_samples'][0]
    row = {'d6': r['d6'], 'arcs': r['arcs'], 'missing_arcs': r['missing_arcs'], 'clique': r['clique'],
           'domination_free': r['domination_free'], 'in_limit_candidates': r['in_limit_candidates_record'],
           'finite_point': {'beta': s0['beta'], 't': s0['t'], 'eps': s0['eps'], 'delta': s0['delta'], 'q': s0['q']},
           'legal': s0['legal'], 'positive': s0['positive'], 'stable': s0['stable'],
           'stable_cell_samples_confirmed': '%d/%d' % (r['n_cell_samples_confirmed'], len(r['finite_cell_samples']))}
    if 'limit_check' in r:
        L = r['limit_check']
        row['limit_point'] = {'beta': L['beta'], 'gamma': L['gamma']}
        row['limit_legal'] = L['legal_limit']; row['limit_w_positive'] = L['w_positive']
        row['limit_shifted_hurwitz'] = L['shifted_hurwitz']
        row['limit_file_witnesses_confirmed'] = '%d/%d' % (r['n_limit_witnesses_confirmed'], len(r['limit_file_witnesses']))
        c['limit_ok'] += L['legal_limit'] and L['w_positive'] and L['shifted_hurwitz'] == 'stable'
        c['lw_conf'] += r['n_limit_witnesses_confirmed']; c['lw_tot'] += len(r['limit_file_witnesses'])
    else:
        row['limit'] = 'never (alpha = -1/2 exactly; factor l^2+l+1)'
    ok = s0['legal'] and s0['positive'] and s0['stable'] == 'stable' and not r['clique'] and r['domination_free']
    row['verdict'] = 'CONFIRMED counterexample' if ok else 'NOT CONFIRMED'
    c['confirmed'] += ok
    c['cells_conf'] += r['n_cell_samples_confirmed']; c['cells_tot'] += len(r['finite_cell_samples'])
    c['domfree'] += r['domination_free']
    rows.append(row)
summary = {'n_graphs': len(rows), 'n_confirmed': c['confirmed'], 'domination_free': c['domfree'],
           'stable_cell_sample_points_confirmed': '%d/%d' % (c['cells_conf'], c['cells_tot']),
           'limit_checks_ok': '%d/123' % c['limit_ok'],
           'limit_file_pointwise_witnesses_confirmed': '%d/%d' % (c['lw_conf'], c['lw_tot']),
           'missing_arcs_distribution': dict(sorted(collections.Counter(r['missing_arcs'] for r in rows).items())),
           'pairwise_non_isomorphic': 'yes (own 5040-perm canonical form, and nauty labelg: 124 distinct)',
           'in_limit_candidates_match': '123 finite records == the 123 limit-file graphs up to isomorphism; FNFL not among them',
           'fnfl_factorization': Fn['sympy_factor_list_of_m_chi_perp'], 'fnfl_quotient_shifted_hurwitz': Fn['quotient_shifted_hurwitz'],
           'fnfl_witnesses': [(w['eps'], w['delta'], w['legal'], w['positive'], w['stable']) for w in Fn['finite_witnesses']],
           'discrepancies': D['problems'], 'geneson_decode': D['geneson']}
json.dump({'summary': summary, 'table': rows}, open('v_structure_table.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
