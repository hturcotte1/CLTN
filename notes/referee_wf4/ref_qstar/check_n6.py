"""Independent recomputation of q_6* from the region2d n=6 HITS certificates.
Uses sympy joint root isolation (not qstar.Alg). Exact only."""
import json, sys
from fractions import Fraction as F
import sympy as sp
x = sp.Symbol('x')
REPO = '/home/user/CLTN'
recs = []
for s in range(3):
    for l in open(f'{REPO}/certificates/region2d_n6_shard{s}of3_HITS.jsonl'):
        if l.strip(): recs.append(json.loads(l))
print('records', [(r['idx'], r['complete']) for r in recs])
qpoly = sp.Poly([797, -1686, 2159, -3004, 2499, -2158, 1521], x)
bpoly = sp.Poly([797, -6468, 22544, -44440, 53280, -38016, 13824], x)
print('bpoly(q+1) == qpoly:', sp.Poly(bpoly.as_expr().subs(x, x + 1), x) == qpoly)
print('qpoly irreducible:', qpoly.is_irreducible, ' roots in [561/481,7/6]:', qpoly.count_roots(sp.Rational(561, 481), sp.Rational(7, 6)))
best = []
for r in recs:
    c = r['cert']
    polys = [sp.Poly(e['minpoly_hi_to_lo'], x) for e in c['beta_critical']]
    prod = sp.Poly(sp.lcm_list([p.as_expr() for p in polys]), x)   # squarefree joint
    prod = sp.Poly(sp.sqf_part(prod.as_expr()), x)
    ivs = sorted(prod.intervals(eps=sp.Rational(1, 10**30)), key=lambda iv: iv[0][0])
    roots = [iv[0] for iv in ivs if iv[0][1] > 2]
    for (lo, hi, b0, tlo, thi) in c['stable_cells']:
        assert thi is None, 't-bounded cell present'
        hiR = sp.Rational(hi)
        loR = sp.Rational(lo)
        # roots strictly inside (lo,hi)?
        inside = [iv for iv in roots if iv[0] > loR and iv[1] < hiR]
        straddle = [iv for iv in roots if iv[0] < hiR < iv[1] or iv[0] < loR < iv[1]]
        assert not inside and not straddle, (inside, straddle)
        above = [iv for iv in roots if iv[0] >= hiR]
        bh = above[0]
        # which polynomial owns it
        own = [p for p in polys if p.count_roots(bh[0], bh[1]) == 1]
        best.append((bh, own[0].all_coeffs(), r['idx'], (lo, hi)))
        print(r['idx'], (lo, hi), 'beta_hi in', [str(bh[0]), str(bh[1])], 'q_sup ~', sp.N((bh[0] + bh[1]) / 2 - 1, 15), 'deg', own[0].degree())
best.sort(key=lambda z: z[0][0])
top = best[-1]
assert all(b[0][1] < top[0][0] or b[0] == top[0] for b in best), 'max not separated'
qa, qb = top[0][0] - 1, top[0][1] - 1
print('max beta_hi-1 interval', qa, qb, 'width', sp.N(qb - qa))
print('qpoly has exactly one root in max interval:', qpoly.count_roots(qa, qb) == 1)
print('max interval inside claimed (561/481,7/6):', qa >= sp.Rational(561, 481) and qb <= sp.Rational(7, 6))
print('attained by', sorted({b[2] for b in best if b[0] == top[0]}))
# no undecided t-cells
for r in recs:
    und = [(row['beta_cell'], tc['t_lo'], tc['t_hi']) for row in r['cert']['cells'] for tc in row['t_cells'] if tc.get('stable_fp') is None]
    print(r['idx'], 'undecided cells:', und)
