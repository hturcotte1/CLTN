"""Referee's own Burnside count: digraphs (up to isomorphism) whose underlying simple graph is G (every edge of G
carries >= 1 arc; non-edges carry none) with >= LO arcs in total.  Isomorphism classes of such digraphs <-> orbits of
Aut(G) on edge-state assignments (3 states per edge), because an isomorphism of two digraphs with the same underlying
graph G is an automorphism of G.  Aut(G) by brute force over all 8! permutations; cycle polynomial per edge cycle."""
import sys, json, itertools, math, os
from collections import Counter
def g6_edges(s):
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - b)) & 1 for b in range(6)]
    E = []; p = 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]: E.append((i, j))
            p += 1
    return n, E
def count(s, lo):
    n, E = g6_edges(s)
    Eset = set(E); idx = {e: k for k, e in enumerate(E)}
    auts = []
    for perm in itertools.permutations(range(n)):
        ok = True
        for (i, j) in E:
            a, b = perm[i], perm[j]
            if (a, b) not in Eset and (b, a) not in Eset: ok = False; break
        if ok: auts.append(perm)
    total = 0
    for perm in auts:
        seen = set(); poly = {0: 1}
        for e in E:
            if e in seen: continue
            L = 0; rev = 0; cur = e
            while cur not in seen:
                seen.add(cur); L += 1
                a, b = perm[cur[0]], perm[cur[1]]
                if a < b: cur = (a, b)
                else: cur = (b, a); rev ^= 1
            # a fixed assignment is constant on the cycle (in the cycle's transported frame); if the cycle returns
            # reversed (rev=1) only the bidirected state is fixed
            cyc = {2 * L: 1} if rev else {L: 2, 2 * L: 1}
            new = {}
            for d1, c1 in poly.items():
                for d2, c2 in cyc.items(): new[d1 + d2] = new.get(d1 + d2, 0) + c1 * c2
            poly = new
        total += sum(c for d, c in poly.items() if d >= lo)
    assert total % len(auts) == 0
    return total // len(auts), len(auts)
if __name__ == '__main__':
    W = os.path.dirname(os.path.abspath(__file__))
    full = [l.strip() for l in open(f'{W}/geng_8_15_28.g6') if l.strip()]
    shard = set()
    for k in range(9): shard |= {l.strip() for l in open(f'{W}/shard{k}.g6') if l.strip()}
    rem = [g for g in full if g not in shard]
    gc = json.load(open('n8/jobs/graph_counts.json'))['graph_counts']
    rc = json.load(open('n8/n8_restricted_counts.json'))['per_shard']
    out = {'sample_rule': 'remaining graphs (geng 8 15:28 order minus shards 0..8), every 66th starting at 0 -> 30 graphs; plus every graph of shards 0 and 4', 'sample': [], 'shards': {}}
    bad = 0
    for i in range(0, 30 * 66, 66):
        g = rem[i]; c, aut = count(g, 29)
        rec = {'g6': g, 'edges': len(g6_edges(g)[1]), 'aut': aut, 'mine': c, 'graph_counts': gc[g], 'agree': c == gc[g]}
        bad += not rec['agree']; out['sample'].append(rec); print(rec, flush=True)
    for k in (0, 4):
        gs = [l.strip() for l in open(f'{W}/shard{k}.g6') if l.strip()]
        s = sum(count(g, 29)[0] for g in gs)
        out['shards'][k] = {'graphs': len(gs), 'mine': s, 'recorded': rc[str(k)]['digraphs'], 'agree': s == rc[str(k)]['digraphs']}
        bad += not out['shards'][k]['agree']; print(k, out['shards'][k], flush=True)
    # validation of the method itself: LO = 0 over all graphs on 5 nodes must give A000273(5) = 9608 (n7/geng 5)
    import subprocess
    g5 = subprocess.run(['n7/geng', '-q', '5'], capture_output=True, text=True).stdout.split()
    tot5 = sum(count(g, 0)[0] for g in g5)
    out['validation_n5_all_digraphs'] = {'mine': tot5, 'A000273': 9608, 'agree': tot5 == 9608}
    print(out['validation_n5_all_digraphs'])
    out['disagreements'] = bad
    json.dump(out, open(f'{W}/burnside_result.json', 'w'), indent=1)
    print('DISAGREEMENTS', bad)
