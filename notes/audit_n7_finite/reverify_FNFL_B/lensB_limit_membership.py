"""Lens B, part 2: fastregion limit verdict; non-membership among the 123 limit candidates."""
import sys, json, subprocess, itertools
sys.path.insert(0, '/home/user/CLTN')
import fastregion
D = sys.argv[1]
rows = ["0011110","0011100","1101111","1110111","0011011","0100101","1100000"]
A = [[int(c) for c in r] for r in rows]
r = fastregion.region(A, list(range(7)), fastregion.INF)
print("fastregion INF verdict:", r['verdict'])
print("  stable_beta:", r.get('stable_beta'))
print("  cells:", json.dumps(r.get('cells'), default=str)[:800])
r2 = fastregion.region(A, list(range(7)), fastregion.INF, crosscheck=True, replay=True)
print("fastregion INF verdict (crosscheck, replay):", r2['verdict'])

def d6_decode(s):
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    return [[bits[i*n+j] for j in range(n)] for i in range(n)]
def d6_encode(M):
    n = len(M)
    bits = [M[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(63+n) + ''.join(chr(63 + int(''.join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))

cands = [json.loads(x) for x in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl') if x.strip()]
print("number of limit candidates:", len(cands), "distinct d6:", len({c['d6'] for c in cands}))
mine = d6_encode(A)
inp = '\n'.join([mine] + [c['d6'] for c in cands]) + '\n'
out = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input=inp, capture_output=True, text=True)
labs = out.stdout.split()
assert len(labs) == 1 + len(cands), (len(labs), out.stderr)
print("labelg canonical of graph:", labs[0], " cert labelg_canonical:", json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))['graph']['labelg_canonical'])
print("labelg: graph canonical form among 123 candidate canonical forms:", labs[0] in set(labs[1:]), "; distinct candidate canonical forms:", len(set(labs[1:])))
# labelg idempotence sanity: canonicalising a canonical form gives itself; a relabelled copy gives the same form
import random
random.seed(1)
perm = list(range(7)); random.shuffle(perm)
Ap = [[A[perm[i]][perm[j]] for j in range(7)] for i in range(7)]
o2 = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input=labs[0] + '\n' + d6_encode(Ap) + '\n', capture_output=True, text=True).stdout.split()
print("labelg sanity (idempotent, relabel-invariant):", o2 == [labs[0], labs[0]])

# Independent brute-force isomorphism test (no nauty): invariants then all 5040 permutations.
def inv(M):
    n = len(M)
    return (sum(map(sum, M)), sorted((sum(M[i]), sum(M[j][i] for j in range(n))) for i in range(n)))
iA = inv(A)
hits = []
for c in cands:
    B = d6_decode(c['d6'])
    if inv(B) != iA:
        continue
    for p in itertools.permutations(range(7)):
        if all(A[p[i]][p[j]] == B[i][j] for i in range(7) for j in range(7)):
            hits.append(c['d6']); break
print("candidates passing invariant filter / isomorphic by brute force:",
      sum(1 for c in cands if inv(d6_decode(c['d6'])) == iA), "/", hits)
print("edge counts of candidates:", sorted({c['edges'] for c in cands}))
