"""Extra candidate universal statements, tested exactly on the n = 8 list and on the n = 6, 7 lists.
Output: more_universals.json"""
import json, os
from collections import Counter
import lib, small_lists as sl

HERE = os.path.dirname(os.path.abspath(__file__))


def props(A):
    n = len(A)
    one = [[1 if (A[i][j] and not A[j][i]) else 0 for j in range(n)] for i in range(n)]
    sinks = [v for v in range(n) if all(A[v][j] for j in range(n) if j != v)]
    non = [v for v in range(n) if v not in sinks]
    onC3 = set()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if len({i, j, k}) == 3 and one[i][j] and one[j][k] and one[k][i]:
                    onC3 |= {i, j, k}
    one_core = lib.induced(one, non)
    return {
        'every_nonsink_node_on_induced_C3': set(non) <= onC3,
        'oneway_digraph_on_nonsinks_strongly_connected': lib.strongly_connected(one_core),
        'every_node_has_oneway_out_or_is_sink_or_bidir_universal': all(
            any(one[v][j] for j in range(n)) or v in sinks for v in range(n)),
        'every_node_has_a_oneway_in_arc_or_is_sink': all(any(one[j][v] for j in range(n)) or v in sinks for v in range(n)),
        'max_outdeg>=5': max(lib.outdeg(A, i) for i in range(n)) >= n - 3,
        'min_outdeg>=2': min(lib.outdeg(A, i) for i in range(n)) >= 2,
        'min_indeg>=2': min(lib.indeg(A, i) for i in range(n)) >= 2,
        'A_has_hamiltonian_cycle': ham(A),
    }


def ham(A):
    n = len(A)
    from functools import lru_cache
    @lru_cache(None)
    def f(mask, v):
        if mask == (1 << n) - 1:
            return bool(A[v][0])
        return any(f(mask | 1 << w, w) for w in range(n) if not mask >> w & 1 and A[v][w])
    return f(1, 0)


def main():
    out = {}
    lists = {'n6': [d for _, d in sl.n6_list()], 'n7': sl.n7_list() + [sl.FNFL],
             'n8': [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]}
    for nm, L in lists.items():
        c = Counter()
        for d in L:
            for k, v in props(lib.from_d6(d)).items():
                c[k] += bool(v)
        out[nm] = dict(total=len(L), holds=dict(c))
    json.dump(out, open(os.path.join(HERE, 'more_universals.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
