import load_hits, graphs, tlib
from collections import Counter
H6 = load_hits.n6_hits()
canon6 = {graphs.canonical_form(A): name for name, A in [('G1', H6[1000255]), ('H3', H6[1387528]), ('G2', H6[1424629])]}
K6 = [[0 if i == j else 1 for j in range(6)] for i in range(6)]
canon6[graphs.canonical_form(K6)] = 'K6'
H7 = load_hits.n7_hits()
cnt = Counter(); rows = []
for r in H7:
    A = r['A']
    subs = []
    for v in range(7):
        keep = [i for i in range(7) if i != v]
        S = [[A[i][j] for j in keep] for i in keep]
        c = graphs.canonical_form(S)
        if c in canon6:
            subs.append((v, canon6[c]))
    cnt[tuple(sorted(set(s for _, s in subs)))] += 1
    rows.append((r['d6'], r['edges'], subs, round(tlib.heur_alpha(A), 4)))
print(cnt)
for x in rows: print(x)
