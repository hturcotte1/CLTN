"""Induced-subgraph census of the 9,133 n = 8 counterexamples (nauty labelg canonical labels, exact).
For k = 3..7: which canonical k-node induced subgraphs occur in EVERY candidate (universal motifs), and
containment of the n = 6 / n = 7 counterexamples (and of K6, K7) as induced subgraphs.
Output: motifs.json, motifs_per_hit.jsonl (per hit: the 7- and 6-node induced counterexamples it contains)."""
import json, os
from itertools import combinations
from collections import Counter
import lib, small_lists as sl

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    hits = [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
    As = [lib.from_d6(d) for d in hits]
    n7c = lib.labelg_canon(sl.n7_list())
    fnflc = lib.labelg_canon([sl.FNFL])[0]
    n6 = sl.n6_list()
    n6c = lib.labelg_canon([d for _, d in n6])
    K = {k: lib.labelg_canon([lib.to_d6([[int(i != j) for j in range(k)] for i in range(k)])])[0] for k in (3, 4, 5, 6, 7)}
    n7set, n6set = set(n7c), set(n6c)
    assert len(n7set) == 123 and len(n6set) == 3
    out = {'universal': {}, 'n_distinct': {}, 'top_frequency': {}}
    per_hit = [dict(d6=d, n7_sub=[], n6_sub=[], fnfl_sub=[], K7=False, K6=0) for d in hits]
    for k in (3, 4, 5, 6, 7):
        subsets = list(combinations(range(8), k))
        d6s = [lib.to_d6(lib.induced(A, s)) for A in As for s in subsets]
        canon = lib.labelg_canon(d6s)
        m = len(subsets)
        freq = Counter()
        universal = None
        for h in range(len(hits)):
            cs = canon[h * m:(h + 1) * m]
            st = set(cs)
            freq.update(st)
            universal = st if universal is None else universal & st
            for s, c in zip(subsets, cs):
                if k == 7 and c in n7set: per_hit[h]['n7_sub'].append([list(s), sl.n7_list()[n7c.index(c)]])
                if k == 7 and c == fnflc: per_hit[h]['fnfl_sub'].append(list(s))
                if k == 6 and c in n6set: per_hit[h]['n6_sub'].append([list(s), n6[n6c.index(c)][0]])
                if k == 7 and c == K[7]: per_hit[h]['K7'] = True
                if k == 6 and c == K[6]: per_hit[h]['K6'] += 1
        out['universal'][k] = sorted(universal)
        out['n_distinct'][k] = len(freq)
        out['top_frequency'][k] = freq.most_common(12)
        print(k, 'distinct', len(freq), 'universal', sorted(universal), flush=True)
    c7 = sum(1 for p in per_hit if p['n7_sub'])
    c6 = sum(1 for p in per_hit if p['n6_sub'])
    out['contains_n7_counterexample'] = c7
    out['contains_FNFL'] = sum(1 for p in per_hit if p['fnfl_sub'])
    out['contains_n6_counterexample'] = c6
    out['contains_n6_or_n7'] = sum(1 for p in per_hit if p['n6_sub'] or p['n7_sub'])
    out['contains_K7'] = sum(1 for p in per_hit if p['K7'])
    out['contains_K6'] = sum(1 for p in per_hit if p['K6'])
    out['n7_counterexamples_used_as_induced_subgraph'] = len({x[1] for p in per_hit for x in p['n7_sub']})
    out['n6_counterexamples_used'] = sorted({x[1] for p in per_hit for x in p['n6_sub']})
    out['K_labels'] = K
    json.dump(out, open(os.path.join(HERE, 'motifs.json'), 'w'), indent=1)
    with open(os.path.join(HERE, 'motifs_per_hit.jsonl'), 'w') as f:
        for p in per_hit:
            f.write(json.dumps(p) + '\n')
    print({k: v for k, v in out.items() if k not in ('top_frequency', 'universal')})


if __name__ == '__main__':
    main()
