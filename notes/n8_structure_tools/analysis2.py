"""Further exact checks on the 9,133 (output analysis2.json):
 A. D-equivalence classes (Lemma J.2 key: charpoly + sorted per-vertex in-walk counts);
 B. psi_D(z) = chi_D - chi_{D+J} (= 8 chi_{D|perp}): multiplicity of the root z = 0 vs number of sink nodes;
 C. graphs with closed twins that are NOT lifts: status of G - u among the n = 7 strict survivors;
 D. automorphisms;
 E. Lemma X2' (sink classes) consistency at the certified interior beta of every candidate:
    w_G > 0, w_core > 0 and sigma_S = sum_{i in S} w_core,i > 1 for every sink class."""
import json, os
from fractions import Fraction as F
from collections import Counter
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def polysub(a, b):
    n = max(len(a), len(b)); a = [0] * (n - len(a)) + a; b = [0] * (n - len(b)) + b
    r = [x - y for x, y in zip(a, b)]
    while r and r[0] == 0: r.pop(0)
    return r


def main():
    feats = [json.loads(l) for l in open(os.path.join(HERE, 'features.jsonl'))]
    fam = json.load(open(os.path.join(HERE, 'families.json')))
    lifts = set(fam['lifts']['indices'])
    W = json.load(open(lib.REPO + '/certificates/n8_candidates_legal_witnesses.json'))['graphs']
    out = {}
    # A
    keys = Counter()
    for f in feats:
        A = lib.from_d6(f['d6'])
        keys[lib.D_class_key(A)] += 1
    sizes = Counter(keys.values())
    out['D_classes'] = dict(count=len(keys), size_distribution=dict(sorted(sizes.items())))
    # B
    zm = Counter()
    for f in feats:
        A = lib.from_d6(f['d6']); D = lib.nonedge(A)
        DJ = [[D[i][j] + 1 for j in range(8)] for i in range(8)]
        psi = polysub(lib.charpoly_int(D), lib.charpoly_int(DJ))
        assert len(psi) == 8 and psi[0] == 8
        mult = 0
        while psi and psi[-1] == 0:
            psi = psi[:-1]; mult += 1
        zm[(len(f['Dsinks']), mult)] += 1
    out['psiD_zero_root_multiplicity_by_sinks'] = {f'sinks={k[0]} mult={k[1]}': v for k, v in sorted(zm.items())}
    # C
    n7res = {}
    for l in open(lib.REPO + '/n7/n7_limit_results.jsonl'):
        r = json.loads(l); n7res[r['d6']] = (r['verdict'], r.get('reason', ''))
    n7keys = list(n7res); n7c = lib.labelg_canon(n7keys); n7map = {c: n7res[d] for c, d in zip(n7c, n7keys)}
    ct = Counter(); ct_list = []
    for i, f in enumerate(feats):
        if f['n_closed_twins'] and i not in lifts:
            A = lib.from_d6(f['d6'])
            u = lib.twin_pairs(A)['closed'][0][0]
            c = lib.labelg_canon([lib.to_d6(lib.induced(A, [x for x in range(8) if x != u]))])[0]
            st = n7map.get(c, ('not a strict n7 survivor', ''))
            key = st[0] + ('' if not st[1] else ' / ' + st[1][:40])
            ct[key] += 1
            ct_list.append(f['d6'])
    out['closed_twin_non_lifts'] = dict(count=len(ct_list), base_status=dict(ct))
    # D
    out['aut'] = dict(Counter(f['aut'] for f in feats))
    out['aut_gt1_without_closed_twins'] = [f['d6'] for f in feats if f['aut'] > 1 and not f['n_closed_twins']]
    # E
    ok = Counter()
    for f, w in zip(feats, W):
        assert w['d6'] == f['d6']
        beta = F(w['witness_in_first_stable_cell']['beta'])
        A = lib.from_d6(f['d6'])
        wG = lib.limit_fixed_point(A, beta)
        ok['w_G>0'] += all(x > 0 for x in wG)
        if f['Dsinks']:
            core = [v for v in range(8) if v not in f['Dsinks']]
            wc = lib.limit_fixed_point(lib.induced(A, core), beta)
            pos = all(x > 0 for x in wc)
            sig = all(sum(wc[core.index(i)] for i in S) > 1 for S, vs in f['sink_classes'])
            ok['with_sinks'] += 1
            ok['w_core>0'] += pos
            ok['all sigma_S>1'] += sig
    out['X2prime_check_at_interior_beta'] = dict(ok)
    json.dump(out, open(os.path.join(HERE, 'analysis2.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
