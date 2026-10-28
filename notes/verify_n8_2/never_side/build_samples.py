"""Deterministic samples.  Rule for every pool: order the pool (stated per pool), N = pool size, take the elements
at indices floor(i*N/300), i = 0..299 (all 300 distinct when N >= 300)."""
import json, sys, os
sys.path.insert(0, '.')
from fractions import Fraction as F
R = '/home/user/CLTN/n8/n8_limit_results.jsonl'
def pick(pool, k=300):
    N = len(pool)
    if N <= k: return list(pool)
    return [pool[(i * N) // k] for i in range(k)]
log = {}
# ---- part 1: never verdicts, pools ordered by Python string order of d6
np_pool, nf_pool = [], []
for line in open(R):
    r = json.loads(line)
    if r['verdict'] != 'never': continue
    (nf_pool if 'cert' in r else np_pool).append(r)
np_pool.sort(key=lambda r: r['d6']); nf_pool.sort(key=lambda r: r['d6'])
log['never_positivity_pool'] = len(np_pool); log['never_full_pool'] = len(nf_pool)
with open('sample_never_positivity.jsonl', 'w') as fh:
    for r in pick(np_pool): fh.write(json.dumps({'d6': r['d6'], 'cls': 'never_positivity', 'reason': r.get('reason'), 'extra_betas': []}) + '\n')
with open('sample_never_full.jsonl', 'w') as fh:
    for r in pick(nf_pool):
        ex = set()
        for c in r['cert']['critical']:
            lo, hi = F(c['interval'][0]), F(c['interval'][1])
            for b in (lo, hi, (lo + hi) / 2): ex.add((b.numerator, b.denominator))
        fh.write(json.dumps({'d6': r['d6'], 'cls': 'never_full', 'n_critical': len(r['cert']['critical']),
                             'cells': [(c['positive'], c['stable']) for c in r['cert']['cells']],
                             'extra_betas': sorted(ex, key=lambda p: F(*p))}) + '\n')
# ---- part 2: filter8 regeneration pools in file order (orbitenum enumeration order), K8-e part 0/300 then GN part 0/1000
printed = {}
B_pool, S_pool = [], []
for f in ['k8e_p0_300.filter8.txt', 'gn_p0_1000.filter8.txt']:
    for line in open(f):
        d6, e, fl = line.split(); printed[d6] = fl
        (B_pool if fl == 'B' else S_pool if fl == 'S' else []).append(d6)
log['B_pool'] = len(B_pool); log['S_pool'] = len(S_pool); log['O_lines'] = sum(v == 'O' for v in printed.values())
dom_files = ['k8e_p0_300.domfree.txt', 'gn_p0_1000.domfree.txt']
N = 0
for f in dom_files:
    for line in open(f):
        if line.split()[0] not in printed: N += 1
log['nonprinted_domfree_pool'] = N
want = {(i * N) // 300 for i in range(300)}; NPDF = []; k = 0
for f in dom_files:
    for line in open(f):
        d6 = line.split()[0]
        if d6 in printed: continue
        if k in want: NPDF.append(d6)
        k += 1
assert len(NPDF) == 300
with open('sample_B.jsonl', 'w') as fh:
    for d in pick(B_pool): fh.write(json.dumps({'d6': d, 'cls': 'B', 'extra_betas': []}) + '\n')
with open('sample_S.jsonl', 'w') as fh:
    for d in pick(S_pool): fh.write(json.dumps({'d6': d, 'cls': 'S', 'extra_betas': []}) + '\n')
with open('sample_NPDF.jsonl', 'w') as fh:
    for d in NPDF: fh.write(json.dumps({'d6': d, 'cls': 'nonprinted_domfree', 'extra_betas': []}) + '\n')
# raw stride lines (every 50000th digraph of the orbitenum stream): classify by filter8 output membership
stride = [l.split()[0] for f in ['k8e_p0_300.stride50000.txt', 'gn_p0_1000.stride50000.txt'] for l in open(f)]
inD = set(l.strip() for l in open('stride_in_domfree.txt'))
NDF = [d for d in stride if d not in inD]
log['stride_lines'] = len(stride); log['stride_in_domfree'] = len(inD); log['stride_nondomfree'] = len(NDF)
with open('sample_NDF.jsonl', 'w') as fh:
    for d in NDF[:300]: fh.write(json.dumps({'d6': d, 'cls': 'nonprinted_nondomfree', 'filter8_flag': printed.get(d)}) + '\n')
with open('sample_stride_domfree.jsonl', 'w') as fh:
    for d in stride:
        if d in inD: fh.write(json.dumps({'d6': d, 'cls': 'stride_domfree', 'filter8_flag': printed.get(d)}) + '\n')
json.dump(log, open('samples_log.json', 'w'), indent=1); print(log)
