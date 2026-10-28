"""(positive_for_all_large_beta / w_inf_positive = the exact positivity set of positivity.py is unbounded)
Exact descriptions of the special subsets of the 9,133 (output special_sets.json):
 the two 29-arc candidates, the 16 with q upper end > 3/2 (two extremal), the 37 narrow-interval graphs,
 the 61 lifts vs the 9,072 others; plus for every candidate the large-beta positivity vector
 w(oo) = (A^T + I)^{-1} 1 (w(beta) -> w(oo) as beta -> oo; w(oo) > 0 iff positive for all large beta)."""
import json, os
from fractions import Fraction as F
from collections import Counter
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def summary(rows, motif, lifts_idx, cat_idx, idxs):
    c = Counter()
    for i in idxs:
        f = rows[i]
        c['n'] += 1
        c['arcs=' + str(f['arcs'])] += 0
    def dist(fn):
        return dict(sorted(Counter(fn(rows[i]) for i in idxs).items(), key=lambda kv: str(kv[0])))
    return {
        'count': len(idxs),
        'arcs': dist(lambda f: f['arcs']),
        'nonadjacent_pairs': dist(lambda f: f['nonadj']),
        'bidirected_pairs': dist(lambda f: f['bi']),
        'omega_bi': dist(lambda f: f['omega_bi']),
        'sink_nodes': dist(lambda f: len(f['Dsinks'])),
        'closed_twin_pairs': dist(lambda f: f['n_closed_twins']),
        'bidirected_in_twin_pairs': dist(lambda f: f['n_bi_in_twins']),
        'aut': dist(lambda f: f['aut']),
        'interval_starts_at_2': sum(1 for i in idxs if rows[i]['lo_kind'] == 'beta_min'),
        'upper_end_spectral_Delta': sum(1 for i in idxs if rows[i]['hi_vanishing'] == ['Delta']),
        'contains_n7_cex': sum(1 for i in idxs if motif[i]['n7_sub']),
        'contains_n6_cex': sum(1 for i in idxs if motif[i]['n6_sub']),
        'catalogue_instance': sum(1 for i in idxs if i in cat_idx),
        'positive_for_all_large_beta': sum(1 for i in idxs if rows[i]['w_inf_positive']),
    }


def main():
    rows = [json.loads(l) for l in open(os.path.join(HERE, 'features.jsonl'))]
    motif = [json.loads(l) for l in open(os.path.join(HERE, 'motifs_per_hit.jsonl'))]
    fam = json.load(open(os.path.join(HERE, 'families.json')))
    Q = {r['d6']: r for r in json.load(open(lib.REPO + '/notes/verify_n8/v_structure/q_upper_ends_refined.json'))['all']}
    lifts = set(fam['lifts']['indices']); cat = set(fam['catalogue_indices'])
    P = [json.loads(l) for l in open(os.path.join(HERE, 'positivity.jsonl'))]
    for f, p in zip(rows, P):
        assert f['d6'] == p['d6']
        f['w_inf_positive'] = p['unbounded']      # exact positivity set unbounded above (positivity.py)
        f['q_hi'] = Q[f['d6']]['q_hi_refined_interval']
        f['narrow37'] = Q[f['d6']]['in_37_no_legal_stored']
    out = {}
    allidx = list(range(len(rows)))
    out['all'] = summary(rows, motif, lifts, cat, allidx)
    two9 = [i for i in allidx if rows[i]['arcs'] == 29]
    out['arcs29'] = []
    for i in two9:
        f = rows[i]
        out['arcs29'].append(dict(d6=f['d6'], A_rows=[''.join(map(str, r)) for r in lib.from_d6(f['d6'])],
                                  bi=f['bi'], single=f['single'], nonadj=f['nonadj'], omega_bi=f['omega_bi'],
                                  outdeg=f['outdeg'], indeg=f['indeg'], sinks=f['Dsinks'], aut=f['aut'],
                                  bi_in_twins=f['n_bi_in_twins'], closed_twins=f['n_closed_twins'],
                                  n7_sub=motif[i]['n7_sub'], n6_sub=motif[i]['n6_sub'],
                                  lo=[f['lo_kind'], f['lo_vanishing']], hi_vanishing=f['hi_vanishing'],
                                  hi_beta_interval=f['hi_beta_interval'], q_hi=f['q_hi'],
                                  w_inf_positive=f['w_inf_positive'], lift=i in lifts, catalogue=i in cat))
    big = [i for i in allidx if F(rows[i]['q_hi'][0]) > F(3, 2)]
    assert all(not (F(rows[i]['q_hi'][0]) <= F(3, 2) <= F(rows[i]['q_hi'][1])) for i in allidx)
    out['q_gt_3_2'] = summary(rows, motif, lifts, cat, big)
    out['q_gt_3_2_list'] = [dict(d6=rows[i]['d6'], q_hi=rows[i]['q_hi'], sinks=len(rows[i]['Dsinks']),
                                 hi=rows[i]['hi_vanishing'], lo=rows[i]['lo_kind'],
                                 w_inf_positive=rows[i]['w_inf_positive'],
                                 Dclass_rep=None) for i in big]
    # D-class relations inside the 16
    keys = {}
    for i in big:
        keys.setdefault(lib.D_class_key(lib.from_d6(rows[i]['d6'])), []).append(rows[i]['d6'])
    out['q_gt_3_2_D_classes'] = [v for v in keys.values()]
    nar = [i for i in allidx if rows[i]['narrow37']]
    out['narrow37'] = summary(rows, motif, lifts, cat, nar)
    out['narrow37_q_hi_max'] = max((rows[i]['q_hi'][1] for i in nar), key=F)
    out['narrow37_list'] = [rows[i]['d6'] for i in nar]
    out['lifts61'] = summary(rows, motif, lifts, cat, sorted(lifts))
    new = [i for i in allidx if i not in lifts]
    out['new9072'] = summary(rows, motif, lifts, cat, new)
    def qmed(idxs):
        qs = sorted(F(rows[i]['q_hi'][0]) for i in idxs)
        return dict(max_lower_bound=str(qs[-1]), median_lower_bound=str(qs[len(qs) // 2]))
    out['lifts61']['q_hi'] = qmed(sorted(lifts)); out['new9072']['q_hi'] = qmed(new)
    json.dump(out, open(os.path.join(HERE, 'special_sets.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1)[:12000])


if __name__ == '__main__':
    main()
