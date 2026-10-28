"""Aggregate the exact per-candidate features (features.jsonl) and the comparison lists (small_features.json)
into structure_summary.json: every distribution and every universal-statement count quoted in
notes/n8_structure.md section 1-2 comes from this file (or the other JSON files named there)."""
import json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
R = [json.loads(l) for l in open(os.path.join(HERE, 'features.jsonl'))]
S = json.load(open(os.path.join(HERE, 'small_features.json')))
mu = json.load(open(os.path.join(HERE, 'more_universals.json')))
sets = {'n6': S['n6'], 'n7_limit': S['n7'], 'n7_FNFL': [S['n7_FNFL']], 'n8': R}

def dist(L, f):
    return {str(k): v for k, v in sorted(Counter(f(r) for r in L).items(), key=lambda kv: str(kv[0]))}

F = {
    'arcs': lambda r: r['arcs'], 'missing_arcs': lambda r: r['missing'],
    'bidirected_pairs': lambda r: r['bi'], 'single_arc_pairs': lambda r: r['single'], 'nonadjacent_pairs': lambda r: r['nonadj'],
    'bi_minus_nonadj': lambda r: r['bi'] - r['nonadj'],
    'omega_bi': lambda r: r['omega_bi'],
    'nodes_sending_to_all': lambda r: len(r['Dsinks']), 'nodes_receiving_from_all': lambda r: len(r['Dsources']),
    'bidirected_universal_nodes': lambda r: len(r['Diso']),
    'A_sources_or_sinks': lambda r: len(r['Asources']) + len(r['Asinks']),
    'D_matching': lambda r: r['D_shape']['matching'], 'D_star': lambda r: r['D_shape']['star'],
    'D_star_forest': lambda r: r['D_shape']['star_forest'], 'D_forest': lambda r: r['D_shape']['forest'],
    'D_max_degree_underlying': lambda r: r['D_shape']['max_degree'], 'D_vertices_covered': lambda r: r['D_shape']['covered'],
    'D_has_directed_cycle': lambda r: r['D_has_cycle'], 'oneway_part_has_directed_cycle': lambda r: r['oneway_has_cycle'],
    'induced_oriented_3cycle': lambda r: r['oneway_3cycles'] > 0,
    'D_strongly_connected': lambda r: r['D_strong'], 'D_on_core_strongly_connected': lambda r: r['Dcore_strong'],
    'A_strongly_connected': lambda r: r['A_strong'],
    'bidirected_in_twin_pairs': lambda r: r['n_bi_in_twins'], 'closed_twin_pairs': lambda r: r['n_closed_twins'],
    'nonadjacent_in_twin_pairs': lambda r: r['n_nonadj_in_twins'], 'dominated_nodes': lambda r: len(r['dominated']),
    'closed_twin_classes_nontrivial': lambda r: str([x for x in r['closed_twin_class_sizes'] if x > 1]),
    'sink_class_sizes': lambda r: str([len(x[1]) for x in r['sink_classes']]),
    'core_size': lambda r: r['core_size'], 'peeling': lambda r: str((r['peel_layers'], r['peel_residual'])),
    'aut': lambda r: r['aut'], 'min_outdeg': lambda r: min(r['outdeg']), 'min_indeg': lambda r: min(r['indeg']),
    'max_outdeg': lambda r: max(r['outdeg']), 'max_indeg': lambda r: max(r['indeg']),
}
out = {nm: {k: dist(L, f) for k, f in F.items()} for nm, L in sets.items()}
n8 = sets['n8']
out['n8_extra'] = dict(
    lo_kind=dist(n8, lambda r: r['lo_kind']),
    hi_spectral=sum(1 for r in n8 if r['hi_vanishing'] == ['Delta']),
    D_strong_iff_no_sink=all(r['D_strong'] == (not r['Dsinks']) for r in n8),
    uncovered_vertices_are_exactly_bidirected_universal=all(8 - r['D_shape']['covered'] == len(r['Diso']) for r in n8),
)
out['more_universals'] = mu
json.dump(out, open(os.path.join(HERE, 'structure_summary.json'), 'w'), indent=1)
print(json.dumps(out['n8_extra']))
for k in ['D_max_degree_underlying', 'D_vertices_covered', 'missing_arcs', 'D_on_core_strongly_connected', 'nodes_receiving_from_all']:
    print(k, out['n8'][k], out['n7_limit'][k])
