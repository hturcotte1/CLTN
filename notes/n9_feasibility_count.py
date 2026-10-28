#!/usr/bin/env python3
"""Exact Polya/Burnside count of isomorphism classes of digraphs on n vertices by arc count.
S_n acts on the n(n-1)/2 unordered pairs; each pair has 4 states (none, u->v, v->u, both) with
arc weights 0, 1, 1, 2.  For a permutation with a pair-cycle of length L: if the cycle returns
the pair with its ends preserved the fixed states contribute 1 + 2x^L + x^{2L}; if the ends are
swapped, only 'none' and 'both' survive: 1 + x^{2L}.  Sum over cycle types weighted by class
size, divide by n!.  Integers only.  Usage: polya_digraphs.py n lo   (count classes with >= lo arcs)
"""
import sys, math
from itertools import combinations
from collections import Counter

def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield []; return
    for k in range(min(n, m), 0, -1):
        for rest in partitions(n - k, k): yield [k] + rest

def class_size(part, n):
    c = Counter(part); d = 1
    for k, mult in c.items(): d *= (k ** mult) * math.factorial(mult)
    return math.factorial(n) // d

def perm_from_partition(part):
    perm = []; start = 0
    for k in part:
        perm.extend([start + (i + 1) % k for i in range(k)]); start += k
    return perm

def polymul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): c[i + j] += x * y
    return c

def fixed_poly(perm, n):
    pairs = list(combinations(range(n), 2)); seen = set(); poly = [1]
    for (u, v) in pairs:
        if (u, v) in seen: continue
        a, b = u, v; L = 0
        while True:
            seen.add((min(a, b), max(a, b))); a, b = perm[a], perm[b]; L += 1
            if {a, b} == {u, v}: swapped = (a != u); break
        cyc = [0] * (2 * L + 1); cyc[0] = 1; cyc[2 * L] += 1
        if not swapped: cyc[L] += 2
        poly = polymul(poly, cyc)
    return poly

def count(n):
    tot = None
    for part in partitions(n):
        p = fixed_poly(perm_from_partition(part), n); cs = class_size(part, n)
        if tot is None: tot = [0] * len(p)
        for i, x in enumerate(p): tot[i] += cs * x
    f = math.factorial(n); assert all(x % f == 0 for x in tot)
    return [x // f for x in tot]   # index = arc count

if __name__ == '__main__':
    n = int(sys.argv[1]); lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    by_arcs = count(n)
    print('n', n, 'all', sum(by_arcs), 'arcs>=%d' % lo, sum(by_arcs[lo:]))
    import json; print(json.dumps({'n': n, 'all_classes': sum(by_arcs), 'min_arcs': lo, 'classes_with_at_least_min_arcs': sum(by_arcs[lo:]), 'by_arc_count': by_arcs}))
