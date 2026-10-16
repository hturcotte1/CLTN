#!/usr/bin/env python3
"""Exact, independent count of the digraphs the n = 8 prefilter must see.

The prefilter runs  geng 8 15:28 k/K | directg -e29:56 | filter8  and filter8 reports how
many digraphs it read ('total' in part<k>.stats).  This script recomputes that number without
directg, by Burnside's lemma: for each underlying graph G (from geng), the number of
orientations with >= LO arcs up to Aut(G) is  (1/|Aut G|) * sum_{g in Aut G} #{orientations with
>= LO arcs fixed by g}.  An orientation gives each edge one of three states (->, <-, <->); g
permutes the edges, and on a cycle of g (length L) the orientation is fixed iff it is constant
along the cycle, which is possible for the two single-arc states only if the cycle preserves
the edge direction (parity of reversals even).  So a cycle contributes the arc-count
generating polynomial  2 x^L + x^{2L}  (preserving) or  x^{2L}  (reversing).

Validation: with LO = 0 and all graphs the count must be A000273(n):
  n = 6 -> 1,540,944   n = 7 -> 882,033,440.
Usage:
  count_restricted.py --validate            (n = 6, 7 checks)
  count_restricted.py --n 8 --edges 15:28 --lo 29 --K 4000 --out n8/n8_restricted_counts.json
All arithmetic is integer; no model, no floats.
"""
import argparse, itertools, json, subprocess, sys, os
import numpy as np

def pair_index(n):
    idx = {}; k = 0
    for u in range(n):
        for v in range(u + 1, n):
            idx[(u, v)] = k; k += 1
    return idx

def graph6_edges(s, n):
    """edge-present bool vector over pairs (u<v) in lexicographic order, from graph6 text."""
    assert ord(s[0]) - 63 == n
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits.extend(((v >> (5 - b)) & 1) for b in range(6))
    idx = pair_index(n)
    e = np.zeros(len(idx), dtype=bool)
    p = 0
    for j in range(1, n):          # graph6 order: column j, rows i < j
        for i in range(j):
            if bits[p]: e[idx[(i, j)]] = True
            p += 1
    return e, idx

def perm_tables(n, idx):
    perms = list(itertools.permutations(range(n)))
    m = len(idx)
    EP = np.zeros((len(perms), m), dtype=np.int16)   # image edge index
    EF = np.zeros((len(perms), m), dtype=np.int8)    # 1 if direction reversed
    pairs = sorted(idx, key=idx.get)
    for a, p in enumerate(perms):
        for (u, v) in pairs:
            pu, pv = p[u], p[v]
            if pu < pv: EP[a, idx[(u, v)]] = idx[(pu, pv)]; EF[a, idx[(u, v)]] = 0
            else:       EP[a, idx[(u, v)]] = idx[(pv, pu)]; EF[a, idx[(u, v)]] = 1
    return EP, EF

def count_graph(e, EP, EF, lo):
    """number of orientations of G (edges e) with >= lo arcs, up to Aut(G); exact int."""
    fixed = np.nonzero((e[EP] == e).all(axis=1))[0]   # automorphisms
    present = np.nonzero(e)[0]
    total = 0
    for a in fixed:
        ep, ef = EP[a], EF[a]
        seen = set(); poly = {0: 1}
        for e0 in present:
            if e0 in seen: continue
            L = 0; par = 0; x = e0
            while x not in seen:
                seen.add(x); par ^= int(ef[x]); x = int(ep[x]); L += 1
            cyc = {2 * L: 1} if par else {L: 2, 2 * L: 1}
            new = {}
            for d1, c1 in poly.items():
                for d2, c2 in cyc.items():
                    new[d1 + d2] = new.get(d1 + d2, 0) + c1 * c2
            poly = new
        total += sum(c for d, c in poly.items() if d >= lo)
    assert total % len(fixed) == 0
    return total // len(fixed), len(fixed)

def geng(n, edges, resmod=None):
    cmd = ['n7/geng', '-q', str(n)] + ([edges] if edges else []) + ([resmod] if resmod else [])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return [l for l in out.split('\n') if l]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--n', type=int, default=8); ap.add_argument('--edges', default='15:28')
    ap.add_argument('--lo', type=int, default=29); ap.add_argument('--K', type=int, default=4000)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.validate:
        for n, want in ((6, 1540944), (7, 882033440)):
            idx = pair_index(n); EP, EF = perm_tables(n, idx)
            tot = sum(count_graph(graph6_edges(g, n)[0], EP, EF, 0)[0] for g in geng(n, None))
            print(f'n={n}: all digraphs = {tot}  (A000273 = {want})  {"OK" if tot == want else "MISMATCH"}')
            if tot != want: sys.exit(1)
        return
    n = a.n
    idx = pair_index(n); EP, EF = perm_tables(n, idx)
    per_shard = {}; grand = 0; ngraphs = 0
    for k in range(a.K):
        s = 0
        gs = geng(n, a.edges, f'{k}/{a.K}')
        for g in gs:
            s += count_graph(graph6_edges(g, n)[0], EP, EF, a.lo)[0]
        per_shard[k] = {'graphs': len(gs), 'digraphs': s}; grand += s; ngraphs += len(gs)
        if k % 400 == 0: print(f'shard {k}: {len(gs)} graphs, {s} digraphs; running total {grand}', flush=True)
    res = {'n': n, 'geng_edges': a.edges, 'min_arcs': a.lo, 'K': a.K, 'underlying_graphs': ngraphs,
           'digraphs_total': grand, 'per_shard': per_shard,
           'method': 'Burnside over Aut(G) per underlying graph; cycle polynomial 2x^L+x^{2L} / x^{2L}'}
    print(json.dumps({k: v for k, v in res.items() if k != 'per_shard'}, indent=1))
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)
if __name__ == '__main__':
    main()
