"""Independent recomputation of claim 16: which of the 123 n=7 candidates contain an n=6 hit (or K6) as an induced
subgraph.  Own digraph6 decoder; own brute-force canonical form (min over all 720 relabellings); cross-checked
against graphs.canonical_form."""
import sys, json, itertools
from collections import Counter
sys.path.insert(0, '/home/user/CLTN')
import graphs
from reprove import d6_decode
PERMS = list(itertools.permutations(range(6)))
def canon_bf(A):
    n = len(A); assert n == 6
    return min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in PERMS)
hits = json.load(open('/home/user/CLTN/certificates/n6_summary.json'))['hits']
names = {1000255: 'G1', 1387528: 'H3', 1424629: 'G2'}
ref_bf = {}; ref_g = {}
for h in hits:
    A = [list(r) for r in h['A']]
    ref_bf[canon_bf(A)] = names[h['idx']]; ref_g[graphs.canonical_form(A)] = names[h['idx']]
K6 = [[0 if i == j else 1 for j in range(6)] for i in range(6)]
ref_bf[canon_bf(K6)] = 'K6'; ref_g[graphs.canonical_form(K6)] = 'K6'
assert len(ref_bf) == 4 and len(ref_g) == 4
# Geneson's graph in the paper labelling must be G1 (or one of the hits)
gen = d6_decode(open('/home/user/CLTN/n7/geneson.d6').read().strip())
print('geneson.d6 is', ref_bf.get(canon_bf(gen)))
cands = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
assert len(cands) == 123
assert len(set(c['d6'] for c in cands)) == 123
cnt = Counter(); rows = []; nemb = 0; emb_T_all = 0; disagreements = 0
for c in cands:
    A = d6_decode(c['d6']); assert A == [list(r) for r in c['A']]
    subs = []
    for v in range(7):
        keep = [i for i in range(7) if i != v]
        Sb = [[A[i][j] for j in keep] for i in keep]
        a = ref_bf.get(canon_bf(Sb)); b = ref_g.get(graphs.canonical_form(Sb))
        if a != b: disagreements += 1
        if a:
            subs.append((v, a)); nemb += 1
            if all(A[v][i] for i in keep): emb_T_all += 1
    key = tuple(sorted(set(s for _, s in subs)))
    cnt[key] += 1
    if subs: rows.append((c['d6'], subs))
print('canonical-form disagreements (own vs graphs.canonical_form):', disagreements)
print('containment classes:', dict(cnt))
print('graphs containing some n=6 hit:', sum(v for kk, v in cnt.items() if kk and 'K6' not in kk))
print('graphs containing K6:', sum(v for kk, v in cnt.items() if 'K6' in kk))
print('embeddings:', nemb, ' of which added node sends to all six:', emb_T_all)
for r in rows: print(r)
json.dump({'classes': {','.join(kk) or 'none': v for kk, v in cnt.items()}, 'rows': rows, 'embeddings': nemb,
           'emb_T_all': emb_T_all, 'disagreements': disagreements}, open('sub6_check.json', 'w'), indent=1)
