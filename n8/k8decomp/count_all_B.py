#!/usr/bin/env python3
"""Dry run of the K_n decomposition with counts only:  for every nonempty B (geng n 1:C), directg -u -o on the
complement H, parse the number of orientations, record "B H e count seconds" in <out>; the sum must equal
the exact Burnside count of the digraphs with underlying graph K_n and at least one bidirected pair.
   python3 count_all_B.py n expected_total out.tsv"""
import sys, os, re, subprocess, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g6tools import complement_g6, edges_g6
n = int(sys.argv[1]); expected = int(sys.argv[2]); out = sys.argv[3]; C = n * (n - 1) // 2
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
Bs = subprocess.run([os.path.join(root, 'n7', 'geng'), '-q', str(n), f'1:{C}'], capture_output=True, text=True, check=True).stdout.split()
tot = 0; t0 = time.time()
with open(out, 'w') as fh:
    for i, B in enumerate(Bs):
        H = complement_g6(B); s = time.time()
        r = subprocess.run([os.path.join(root, 'n7', 'directg'), '-u', '-o'], input=H + '\n', capture_output=True, text=True)
        m = re.search(r'(\d+) digraphs generated', r.stderr)
        if r.returncode != 0 or not m: print('FAILED', B, r.stderr); sys.exit(1)
        c = int(m.group(1)); tot += c
        fh.write(f'{B} {H} {len(edges_g6(B))} {c} {time.time()-s:.3f}\n'); fh.flush()
        if (i + 1) % 500 == 0: print(f'{i+1}/{len(Bs)} B done, running total {tot}, {time.time()-t0:.0f}s', flush=True)
print(f'DONE {len(Bs)} B, total {tot}, expected {expected}, {"MATCH" if tot == expected else "MISMATCH"}, {time.time()-t0:.0f}s')
