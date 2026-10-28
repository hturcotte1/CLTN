import sys, json, random, time, subprocess
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from fractions import Fraction as F
random.seed(1)
hits = [json.loads(l) for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
n7 = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
# roundtrip
for r in hits[:2000] + n7:
    A = d6_decode(r['d6']); assert d6_encode(A) == r['d6']; assert arcs(A) == r['edges'], (r['d6'], arcs(A), r['edges'])
print('roundtrip ok; n7 A field check:', all(d6_decode(r['d6']) == r['A'] for r in n7 if 'A' in r))
# canon invariance + cost
t = time.time(); costs = []
for r in hits[:300]:
    A = d6_decode(r['d6']); c = canon(A); costs.append(canon_cost(A))
    p = list(range(8)); random.shuffle(p); assert canon(permute(A, p)) == c
print('canon 300 graphs + 300 permuted:', round(time.time() - t, 1), 's; max cost', max(costs), 'mean', sum(costs) / len(costs))
# canon classes vs labelg classes on a sample with random relabellings
lines = []; truth = []
for r in hits[:200]:
    A = d6_decode(r['d6'])
    for k in range(3):
        p = list(range(8)); random.shuffle(p); B = permute(A, p); lines.append(d6_encode(B)); truth.append(canon(A))
out = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
own = [canon(d6_decode(l)) for l in lines]
# partitions agree?
from collections import defaultdict
pa = defaultdict(set); pb = defaultdict(set)
for i, (a, b) in enumerate(zip(own, out)): pa[a].add(i); pb[b].add(i)
print('partition agreement with labelg:', sorted(map(sorted, pa.values())) == sorted(map(sorted, pb.values())), len(pa), len(pb))
# decision engine vs stored witnesses
t = time.time(); nchk = 0; bad = 0
for r in hits[:40]:
    A = d6_decode(r['d6'])
    for w in r['witnesses']:
        eps = F(w['eps']); delta = F(w['delta']); d = decide(A, eps, delta); nchk += 1
        st = d['stability'] == 'stable'
        if st != w['stable'] or d['positive'] != w['positive']: bad += 1; print('MISMATCH', r['d6'], w, d['stability'], d['positive'])
print('witness checks', nchk, 'mismatches', bad, round(time.time() - t, 1), 's')
# Geneson
G = json.load(open('/home/user/CLTN/certificates/stage0_geneson_2607.21396.json')) if False else None
# never control: a never record at its interior beta
res = open('/home/user/CLTN/n8/n8_limit_results.jsonl')
cnt = 0
for l in res:
    r = json.loads(l)
    if r['verdict'] == 'never':
        A = d6_decode(r['d6']); d = decide(A, F(1, 10001), F(1001, 10001000)); print('never ctrl', r['d6'], d['stability'], d['positive']); cnt += 1
        if cnt >= 3: break
# K8
K = [[int(i != j) for j in range(8)] for i in range(8)]; d = decide(K, F(1, 3), F(2)); print('K8', d['stability'], d['positive'])
