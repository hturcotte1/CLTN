#!/usr/bin/env python3
"""Enumerate the bidirected subgraphs B for the K_n decomposition and pack them into batches.
   python3 gen_batches.py n emin:emax K outdir [counts.tsv]
geng -q n emin:emax lists every graph B on n vertices with emin..emax edges up to isomorphism (nauty, exactly
one representative per class).  Each line of outdir/batch<k>.txt is "B H e": B in graph6, H = complement of
B in graph6 (same labelled vertices), e = |E(B)|.  Without counts.tsv the batches are balanced greedily on
the estimate 2^(C-e) (C = n(n-1)/2).  With counts.tsv (lines "B H e count secs" from count_all_B.py, the
count-only dry run) they are balanced on the measured cost secs + count/5e5 and the manifest records each
batch's EXACT expected digraph count (sum of the dry-run counts of its B), which the driver checks.
manifest.json records the plan.  Integers only (the cost weights are only used for the packing order)."""
import sys, os, json, subprocess, heapq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g6tools import complement_g6, edges_g6
n = int(sys.argv[1]); erange = sys.argv[2]; K = int(sys.argv[3]); out = sys.argv[4]
counts_file = sys.argv[5] if len(sys.argv) > 5 else None
C = n * (n - 1) // 2
geng = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'n7', 'geng')
Bs = subprocess.run([geng, '-q', str(n), erange], capture_output=True, text=True, check=True).stdout.split()
assert len(set(Bs)) == len(Bs)
counts = {}
if counts_file:
    for line in open(counts_file):
        B, H, e, c, secs = line.split(); counts[B] = (int(c), float(secs))
        assert H == complement_g6(B) and int(e) == len(edges_g6(B))
    assert set(counts) == set(Bs), 'counts file does not cover exactly the geng list'
    items = sorted(((counts[B][1] + counts[B][0] / 5e5, B) for B in Bs), reverse=True)
else:
    items = sorted(((2 ** (C - len(edges_g6(B))), B) for B in Bs), reverse=True)
heap = [(0, k) for k in range(K)]; heapq.heapify(heap); batches = [[] for _ in range(K)]; load = [0] * K
for w, B in items:
    l, k = heapq.heappop(heap); batches[k].append(B); load[k] = l + w; heapq.heappush(heap, (l + w, k))
os.makedirs(out, exist_ok=True)
for k in range(K):
    with open(os.path.join(out, f'batch{k}.txt'), 'w') as fh:
        for B in batches[k]: fh.write(f'{B} {complement_g6(B)} {len(edges_g6(B))}\n')
man = {'n': n, 'edges': erange, 'K': K, 'n_B': len(Bs), 'geng': ' '.join([geng, '-q', str(n), erange]),
       'batch_sizes': [len(b) for b in batches], 'batch_weights': [round(l, 3) for l in load], 'counts_file': counts_file}
if counts:
    man['batch_expected'] = [sum(counts[B][0] for B in b) for b in batches]
    man['total_expected'] = sum(man['batch_expected'])
    man['batch_est_seconds'] = [round(sum(counts[B][1] for B in b), 1) for b in batches]
json.dump(man, open(os.path.join(out, 'manifest.json'), 'w'), indent=1)
print(f'{len(Bs)} graphs B, {K} batches, sizes {min(len(b) for b in batches)}..{max(len(b) for b in batches)}, weights {min(load):.1f}..{max(load):.1f}' + (f", exact expected total {man['total_expected']}, expected per batch {min(man['batch_expected'])}..{max(man['batch_expected'])}" if counts else ''))
