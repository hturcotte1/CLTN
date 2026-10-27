"""Unit test: addarcs (C) against the pure-Python reference on random digraphs and random B, n = 2..8;
also round-trip tests of the graph6/digraph6 codecs against nauty's own output."""
import random, subprocess, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g6tools import *
here = os.path.dirname(os.path.abspath(__file__))
random.seed(20261027); bad = 0; cases = 0
for n in range(2, 9):
    for _ in range(60):
        D = [[random.randint(0, 1) if i != j else 0 for j in range(n)] for i in range(n)]
        Badj = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                Badj[i][j] = Badj[j][i] = random.randint(0, 1)
        d6 = encode_d6(n, D); b6 = encode_g6(n, Badj)
        exp = add_bidirected(d6, b6)
        got = subprocess.run([os.path.join(here, 'addarcs'), b6], input=d6 + '\n', capture_output=True, text=True)
        cases += 1
        if got.returncode != 0 or got.stdout.strip() != exp: bad += 1; print('MISMATCH', n, d6, b6, exp, got.stdout.strip(), got.stderr.strip())
# codec round trips against nauty: geng graph6 lines decode/encode identically; directg digraph6 lines too
for n in (4, 5, 6):
    g = subprocess.run([os.path.join(here, '../../n7/geng'), '-q', str(n)], capture_output=True, text=True).stdout.split()
    for s in g:
        cases += 1
        if encode_g6(*decode_g6(s)) != s: bad += 1; print('g6 roundtrip fail', s)
    d = subprocess.run([os.path.join(here, '../../n7/directg'), '-q'], input=g[-1] + '\n', capture_output=True, text=True).stdout.split()
    for s in d[:2000]:
        cases += 1
        if encode_d6(*decode_d6(s)) != s: bad += 1; print('d6 roundtrip fail', s)
print(f'addarcs/codec tests: {cases} cases, {bad} failures'); sys.exit(1 if bad else 0)
