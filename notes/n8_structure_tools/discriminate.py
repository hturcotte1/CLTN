"""Do the universal properties of the 9,133 also hold for the 338,505 strict spectral survivors that are
NOT limit counterexamples ('never')?  This measures whether a property is a consequence of the prefilter
(domination-free + alpha(A|perp) < -1/2) or genuinely discriminates the counterexamples.
Output: discriminate.json.  Exact integer code only."""
import json, os, sys
from multiprocessing import Pool
from collections import Counter
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def props(line):
    r = json.loads(line)
    if r.get('complete'):
        return None
    A = lib.from_d6(r['d6'])
    n = len(A)
    D = lib.nonedge(A)
    one = [[1 if (A[i][j] and not A[j][i]) else 0 for j in range(n)] for i in range(n)]
    ow3 = any(one[i][j] and one[j][k] and one[k][i]
              for i in range(n) for j in range(n) for k in range(n) if i < j and i < k)
    b, s, z = lib.pair_types(A)
    Dsrc = [v for v in range(n) if lib.indeg(D, v) == 0]
    src_iso = all(lib.outdeg(D, v) == 0 for v in Dsrc)
    omega, _ = lib.bidirected_clique_number(A)
    nsinks = sum(1 for v in range(n) if lib.outdeg(D, v) == 0)
    return (r['verdict'], r.get('reason', ''), dict(
        ow3=ow3, Astrong=lib.strongly_connected(A), oneway_cycle=lib.has_directed_cycle(one),
        src_iso=src_iso, omega_le5=omega <= 5, omega=omega, nonadj=z, nsinks=nsinks,
        nonadj_in_twins=len(lib.twin_pairs(A)['nonadj_in'])))


if __name__ == '__main__':
    lines = open(lib.REPO + '/n8/n8_limit_results.jsonl').read().splitlines()
    with Pool(2) as p:
        res = [x for x in p.imap(props, lines, chunksize=500) if x is not None]
    out = {}
    for verdict, reason, d in res:
        key = 'hit' if verdict == 'stable_somewhere' else ('never_positivity' if 'a_i' in reason else 'never_full')
        o = out.setdefault(key, {'total': 0})
        o['total'] += 1
        for k, v in d.items():
            o.setdefault(k, Counter())[str(v)] += 1
    out = {k: {kk: (dict(vv) if isinstance(vv, Counter) else vv) for kk, vv in v.items()} for k, v in out.items()}
    json.dump(out, open(os.path.join(HERE, 'discriminate.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))
