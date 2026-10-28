"""Vertex-deletion heredity: for every n = 8 candidate G and every vertex v, classify G - v against the
complete n = 7 limit map (n7/n7_limit_results.jsonl = all 1,889 strict n = 7 survivors with verdicts; anything
else failed domination-freeness or the strict spectral condition alpha < -1/2).  Also: for embeddings of an
n = 7 counterexample, does the extra node send to all seven (T = all)?  Output: heredity.json"""
import json, os
from collections import Counter
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    hits = [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
    res = {}
    for l in open(lib.REPO + '/n7/n7_limit_results.jsonl'):
        r = json.loads(l)
        v = r['verdict']
        if v == 'stable_somewhere':
            tag = 'K7' if r.get('complete') else 'cex'
        else:
            tag = 'never_positivity' if 'a_i' in r.get('reason', '') else 'never_full'
        res[r['d6']] = tag
    keys = list(res); kc = lib.labelg_canon(keys); cmap = {c: res[d] for c, d in zip(kc, keys)}
    subs = []
    for d in hits:
        A = lib.from_d6(d)
        for v in range(8):
            subs.append(lib.to_d6(lib.induced(A, [x for x in range(8) if x != v])))
    sc = lib.labelg_canon(subs)
    per = Counter(); best = Counter(); Tall = Counter()
    for h, d in enumerate(hits):
        A = lib.from_d6(d)
        tags = []
        for v in range(8):
            t = cmap.get(sc[8 * h + v], 'not_strict')
            tags.append(t)
            if t == 'cex':
                Tall['extra node sends to all 7' if all(A[v][j] for j in range(8) if j != v) else 'extra node misses an out-arc'] += 1
        per[sum(1 for t in tags if t != 'not_strict')] += 1
        order = ['cex', 'never_full', 'never_positivity', 'not_strict']
        best[min(tags, key=order.index)] += 1
    out = dict(number_of_deletions_that_are_strict_n7_survivors=dict(sorted(per.items())),
               best_deletion_class=dict(best), n7_cex_embeddings_extra_node=dict(Tall))
    json.dump(out, open(os.path.join(HERE, 'heredity.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
