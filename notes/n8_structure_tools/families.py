"""Which n = 8 counterexamples are instances of the proven constructions?
 (1) the 61 closed-twin blow-ups of n = 7 counterexamples found in the list (notes/verify_n8_2/completeness_twins);
 (2) Theorem F's n = 8 member F_2;
 (3) the n = 8 members of the proved clique-class catalogue (notes/research_wf3/theory/families_catalogue_n6/_n7.json):
     n = 6 base with k = 2, n = 7 base with k = 1 (G[S,all;k]: class bidirected, sends to all, i -> class iff S_i = 1);
 (4) G[S,T;2] over an n = 6 counterexample (a closed-twin pair whose deletion leaves an n = 6 counterexample);
 (5) sink structure: core = G minus all nodes that send to every other node; is the core / G - v a smaller
     counterexample?
All identifications by nauty labelg canonical labels (exact).  Output: families.json"""
import json, os
from collections import Counter
import lib, small_lists as sl

HERE = os.path.dirname(os.path.abspath(__file__))
TH = lib.REPO + '/notes/research_wf3/theory/'


def class_ext(G, S, k, T=None):
    m = len(G); n = m + k
    T = T or [1] * m
    A = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            A[i][j] = G[i][j]
    for v in range(m, n):
        for i in range(m):
            A[i][v] = S[i]; A[v][i] = T[i]
        for w in range(m, n):
            if w != v: A[v][w] = 1
    return A


def main():
    hits = [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
    feats = [json.loads(l) for l in open(os.path.join(HERE, 'features.jsonl'))]
    hc = lib.labelg_canon(hits)
    assert len(set(hc)) == 9133
    idx = {c: i for i, c in enumerate(hc)}
    n7 = sl.n7_list(); n7c = lib.labelg_canon(n7); n7set = set(n7c)
    n6 = sl.n6_list(); n6c = lib.labelg_canon([d for _, d in n6]); n6set = set(n6c)
    out = {}
    # (1) lifts
    B = [json.loads(l) for l in open(lib.REPO + '/notes/verify_n8_2/completeness_twins/part2_blowups.jsonl')]
    lift_d6 = sorted({r['blowup_d6'] for r in B if r['in_hits_own']})
    lift_c = lib.labelg_canon(lift_d6)
    lifts = sorted({idx[c] for c in lift_c if c in idx})
    out['lifts'] = dict(count=len(lifts), indices=lifts,
                        parents=len({r['parent_d6'] for r in B if r['in_hits_own']}))
    # own recomputation of lifts: closed-twin blow-ups of all 123 at every node
    own = set()
    for d in n7:
        G = lib.from_d6(d)
        for v in range(7):
            A = class_ext([row[:] for row in G], [G[i][v] for i in range(7)], 1, [G[v][i] for i in range(7)])
            A[v][7] = A[7][v] = 1
            own.add(lib.to_d6(A))
    own = sorted(own); ownc = lib.labelg_canon(own)
    own_in = sorted({idx[c] for c in ownc if c in idx})
    out['lifts']['own_recount_distinct_blowups'] = len(set(ownc))
    out['lifts']['own_recount_in_list'] = len(own_in)
    out['lifts']['own_equals_campaign'] = own_in == lifts
    # (2) Theorem F, F_2
    F2 = lib.labelg_canon(['&GZy~VwuKW~^w'])[0]
    out['F2_in_list'] = F2 in idx
    out['F2_index'] = idx.get(F2)
    # (3) catalogue
    cat = []
    for fn, k in (('families_catalogue_n6.json', 2), ('families_catalogue_n7.json', 1)):
        for fam in json.load(open(TH + fn)):
            base = fam['base'].split(':')[-1]
            A = class_ext(lib.from_d6(base), fam['S'], k)
            cat.append((fn, base, fam['S'], k, lib.to_d6(A), fam['proved_for_all_k>=1']))
    catc = lib.labelg_canon([c[4] for c in cat])
    inl = [c in idx for c in catc]
    out['catalogue'] = dict(instances=len(cat), proved_for_k_used=sum(1 for c in cat if c[5]),
                            in_list=sum(inl), distinct_graphs=len(set(catc)),
                            distinct_in_list=len({c for c, f in zip(catc, inl) if f}),
                            not_in_list=[c[:4] for c, f in zip(cat, inl) if not f])
    cat_idx = sorted({idx[c] for c in catc if c in idx})
    # (4) closed-twin pair over n6, and closed-twin node over n7
    g62, tw7 = [], []
    for i, d in enumerate(hits):
        A = lib.from_d6(d)
        pairs = lib.twin_pairs(A)['closed']
        for u, v in pairs:
            rest = [x for x in range(8) if x not in (u, v)]
            if lib.labelg_canon([lib.to_d6(lib.induced(A, rest))])[0] in n6set:
                g62.append(i); break
        for u, v in pairs:
            rest = [x for x in range(8) if x != u]
            if lib.labelg_canon([lib.to_d6(lib.induced(A, rest))])[0] in n7set:
                tw7.append(i); break
    out['closed_twin_pair_over_n6'] = sorted(set(g62))
    out['closed_twin_over_n7'] = len(set(tw7))
    out['closed_twin_over_n7_equals_lifts'] = sorted(set(tw7)) == lifts
    # (5) sinks / cores
    core_d6 = [f['core_d6'] for f in feats]
    has_core = [i for i, f in enumerate(feats) if f['core_size'] < 8]
    cc = lib.labelg_canon([core_d6[i] for i in has_core])
    core_stat = Counter()
    core_is = {}
    for i, c in zip(has_core, cc):
        sz = feats[i]['core_size']
        tag = 'n7_cex' if c in n7set else ('n6_cex' if c in n6set else 'not_a_cex')
        core_stat[(sz, tag)] += 1
        core_is[i] = tag
    out['core_membership'] = {f'{k[0]}:{k[1]}': v for k, v in sorted(core_stat.items())}
    out['distinct_cores_by_size'] = {}
    for sz in (5, 6, 7):
        out['distinct_cores_by_size'][sz] = len({c for i, c in zip(has_core, cc) if feats[i]['core_size'] == sz})
    # G - v for every sink v
    gv = Counter(); any_gv = 0
    for i in has_core:
        A = lib.from_d6(hits[i])
        labs = lib.labelg_canon([lib.to_d6(lib.induced(A, [x for x in range(8) if x != v])) for v in feats[i]['Dsinks']])
        flags = [c in n7set for c in labs]
        gv[(len(flags), sum(flags))] += 1
        any_gv += any(flags)
    out['delete_one_sink'] = {f'{k[0]} sinks, {k[1]} deletions give an n7 cex': v for k, v in sorted(gv.items())}
    out['delete_some_sink_gives_n7_cex'] = any_gv
    # sink in-set deficiency (number of core nodes NOT sending to the sink) and sink class sizes
    defi = Counter()
    for f in feats:
        for S, vs in f['sink_classes']:
            defi[f['core_size'] - len(S)] += len(vs)
    out['sink_node_core_nonsenders'] = dict(sorted(defi.items()))
    out['catalogue_indices'] = cat_idx
    json.dump(out, open(os.path.join(HERE, 'families.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ('catalogue_indices',)}, indent=1)[:5000])


if __name__ == '__main__':
    main()
