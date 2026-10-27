#!/usr/bin/env python3
"""Independent Burnside count of the isomorphism classes of digraphs whose underlying graph is K_n and that
have at least one bidirected pair (3 states per pair: i->j, j->i, both; tournaments excluded).
A permutation pi of the vertices fixes a digraph iff on every cycle of pi acting on the pairs the states are
consistent: if pi^k returns the pair {u,v} to itself with u,v fixed, any of the 3 states can be carried
round the cycle (3 choices); if it returns it swapped, only 'both' survives (1 choice).  Tournaments: 2 / 0.
   python3 polya_count.py n"""
import sys, itertools
n = int(sys.argv[1]); pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]; idx = {p: k for k, p in enumerate(pairs)}
tot3 = 0; tot2 = 0
for perm in itertools.permutations(range(n)):
    seen = [False] * len(pairs); f3 = 1; f2 = 1
    for k, (u, v) in enumerate(pairs):
        if seen[k]: continue
        a, b = u, v; flip = False
        while True:
            seen[idx[(min(a, b), max(a, b))]] = True
            a, b = perm[a], perm[b]
            if {a, b} == {u, v}:
                flip = (a != u); break
        if flip: f3 *= 1; f2 *= 0
        else: f3 *= 3; f2 *= 2
    tot3 += f3; tot2 += f2
import math
g = math.factorial(n); assert tot3 % g == 0 and tot2 % g == 0
print(f'n={n}: all 3-state classes {tot3//g}, tournaments {tot2//g}, with >=1 bidirected pair {tot3//g - tot2//g}')
