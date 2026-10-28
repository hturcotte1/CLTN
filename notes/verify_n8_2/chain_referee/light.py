import json, sys, collections
from fractions import Fraction as Fr
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/chain_referee')
from myfilter2 import decode_d6, dom_free, matmul, charpoly, hurwitz_minors, det
from math import comb
W = '/home/user/CLTN/notes/verify_n8_2/chain_referee'
# (1) the 23 B(filter8) vs N(mine) lines: show the first zero minor index and the first negative minor index
f8 = {l.split()[0]: l.split()[2] for l in open(f'{W}/variants_filter8.txt')}
mine = {l.split()[0]: l.split()[2] for l in open(f'{W}/variants_mine.txt')}
expl = collections.Counter()
for d6 in mine:
    if f8.get(d6, 'N') == 'B' and mine[d6] == 'N':
        n, A = decode_d6(d6)
        P = [[(Fr(1) if i == j else Fr(0)) - Fr(1, n) for j in range(n)] for i in range(n)]
        c = charpoly(matmul(matmul(P, [[Fr(x) for x in r] for r in A]), P))[:-1]; d = n - 1
        h = [Fr(0)] * (d + 1)
        for k, pk in enumerate(c):
            e = d - k
            for j in range(e + 1): h[j] += pk * comb(e, j) * Fr(-1, 2) ** (e - j)
        hc = list(reversed(h)); m = hurwitz_minors(hc)
        z = next((i for i, x in enumerate(m) if x == 0), None); ng = next((i for i, x in enumerate(m) if x < 0), None)
        expl[(z, ng)] += 1
print('(1) B-vs-N: (index of first zero minor, index of first negative minor):', dict(expl))
# (2) K8 per-B ledger
tot = 0; nB = 0
for l in open('n8/k8decomp/k8_perB_counts.tsv'):
    p = l.split()
    if not p or not p[-1].isdigit(): continue
    tot += int(p[-1]); nB += 1
print('(2) k8_perB_counts.tsv rows', nB, 'sum', tot, 'K8 exact 575009339:', tot == 575009339)
print(open('n8/k8decomp/k8_perB_counts.tsv').readline().rstrip())
# (3) referees json
r = json.load(open('notes/k8_decomposition_referees.json'))
print('(3) referees:', [(k, (v.get('verdict') if isinstance(v, dict) else str(v)[:60])) for k, v in (r.items() if isinstance(r, dict) else enumerate(r))][:8])
# (4) summary & witnesses
s = json.load(open('certificates/n8_limit_summary.json'))
print('(4) summary:', {k: s[k] for k in s if k not in ('edge_counts',)})
w = json.load(open('certificates/n8_candidates_legal_witnesses.json'))
print('witness file type', type(w).__name__, 'len', len(w), 'keys' if isinstance(w, dict) else '', list(w.keys())[:5] if isinstance(w, dict) else '')
