#!/usr/bin/env python3
"""Post-run ledger check for the K8 bidirected-subgraph batches: for every batch job (split decomp:...) in
n8/jobs/graph_counts.json, (1) the per-B line counts written by addarcs into job<j>.genlog must equal the
count-only dry run's counts (n8/k8decomp/k8_perB_counts.tsv) for exactly the B of that batch, (2) their sum
must equal the filter8 'total' in job<j>.stats and the job's expected value, and (3) over all batches every
one of the 12,345 graphs B must occur exactly once and the grand sum must be the exact K8 count.
Exit 0 iff all hold.  Integers only."""
import json, re, os, sys
gc = json.load(open('n8/jobs/graph_counts.json')); exact = gc['graph_counts']['G~~~~{']
dry = {}
for line in open('n8/k8decomp/k8_perB_counts.tsv'):
    B, H, e, c, s = line.split(); dry[B] = int(c)
seen = {}; ok = True; grand = 0; nb = 0
for j in gc['jobs']:
    if not str(j['split']).startswith('decomp:'): continue
    nb += 1; jid = j['job_id']; batch = j['split'][len('decomp:'):]
    Bs = [l.split()[0] for l in open(batch)]
    gen = dict(re.findall(r'^addarcs (\S+) lines (\d+)$', open(f'n8/jobs/job{jid}.genlog').read(), re.M))
    if set(gen) != set(Bs): ok = False; print(f'job {jid}: genlog B set != batch B set'); continue
    tot = 0
    for B in Bs:
        n = int(gen[B]); tot += n
        if dry.get(B) != n: ok = False; print(f'job {jid}: B {B} addarcs lines {n} != dry run {dry.get(B)}')
        seen[B] = seen.get(B, 0) + 1
    m = re.search(r'total (\d+)', open(f'n8/jobs/job{jid}.stats').read()); st = int(m.group(1)) if m else None
    if st != tot or tot != j['expected']: ok = False; print(f'job {jid}: stats total {st}, genlog sum {tot}, expected {j["expected"]}')
    grand += tot
if set(seen) != set(dry) or any(v != 1 for v in seen.values()): ok = False; print('B coverage is not exactly once each')
if grand != exact: ok = False; print(f'grand sum {grand} != exact {exact}')
print(f'k8decomp ledger: {nb} batches, {len(seen)} graphs B, grand sum {grand}, exact {exact}, {"OK" if ok else "PROBLEMS"}')
sys.exit(0 if ok else 1)
