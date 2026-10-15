"""Brute-force isomorphism classes of the 124 candidate graphs (own code, 7! perms)."""
import sys, json, itertools, types
src = open(__file__.replace('iso_check.py', 'verify.py')).read().replace('\nmain()\n', '\n')
mod = types.ModuleType('v'); exec(compile(src, 'verify.py', 'exec'), mod.__dict__)
ROOT = sys.argv[1]
fin = [json.loads(l) for l in open(f'{ROOT}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
def canon(A):
    n = len(A); best = None
    for p in itertools.permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key > best: best = key
    return best
cl = {}
for r in fin:
    cl.setdefault(canon(mod.decode_d6(r['d6'])), []).append(r['d6'])
dups = [v for v in cl.values() if len(v) > 1]
print(json.dumps({'n_graphs': len(fin), 'n_iso_classes': len(cl), 'duplicate_classes': dups}))
