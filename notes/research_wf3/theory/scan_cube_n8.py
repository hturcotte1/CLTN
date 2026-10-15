"""Exact triage of the n=8 one-node extensions flagged (floats) as having an eigenvalue of A|perp at a primitive
cube root of unity with alpha ~ -1/2.  Keep those with: l^2+l+1 exactly once in n chi_perp, every other root
Re < -1/2 (exact Hurwitz), the root is moving (not a root of chi), R* > 1/8 (exact), positivity at beta = 2 + 1e-4."""
import json, sys, subprocess
from fractions import Fraction as F
import tlib, boundary_tools as bt, rstar, certify_top
d = json.load(open('ext_n8.json'))
H7 = {r['d6']: r['A'] for r in __import__('load_hits').n7_hits()}
H7['&FNFL~[udo?'] = tlib.from_d6('&FNFL~[udo?')
graphs_ = [certify_top.ext(H7[d6], pat) for d6, pat in d['cube']]
d6s = [tlib.to_d6(B) for B in graphs_]
p = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input='\n'.join(d6s) + '\n', capture_output=True, text=True)
cans = p.stdout.split()
uniq = {}
for B, c in zip(graphs_, cans): uniq.setdefault(c, B)
print('flagged', len(graphs_), 'distinct', len(uniq), flush=True)
from collections import Counter
cnt = Counter(); keep = []
for c, B in uniq.items():
    mult, ok = bt.boundary_split(B)
    if mult != 1: cnt['mult%d' % mult] += 1; continue
    if not ok: cnt['others_not_strict'] += 1; continue
    pos = tlib.positive_at(B, F(2) + F(1, 10**4))
    if not pos: cnt['not_pos_near_2'] += 1; continue
    R, val = rstar.R_at(B)
    if R is None: cnt['R_undefined:' + str(val)] += 1; continue
    if R <= F(1, 8): cnt['R<=1/8'] += 1; continue
    cnt['GOOD'] += 1; keep.append((c, tlib.to_d6(B), str(R)))
    print('GOOD', c, tlib.to_d6(B), 'R* =', R, flush=True)
print(cnt)
json.dump(keep, open('cube_n8_good.json', 'w'), indent=1)
