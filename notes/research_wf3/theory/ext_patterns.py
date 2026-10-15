# For each n=7 candidate containing an n=6 hit as induced subgraph G7 - v, express v's in/out pattern in the
# coordinates of the n=6 hit (all isomorphisms), and check whether v is a twin of some node.
import load_hits, graphs, tlib
from itertools import permutations
H6 = load_hits.n6_hits()
named = {'G1': H6[1000255], 'H3': H6[1387528], 'G2': H6[1424629]}
H7 = load_hits.n7_hits()
def isos(S, T):
    n = len(S)
    for p in permutations(range(n)):
        if all(S[p[i]][p[j]] == T[i][j] for i in range(n) for j in range(n)):
            yield p   # T node i = S node p[i]
for r in H7:
    A = r['A']
    for v in range(7):
        keep = [i for i in range(7) if i != v]
        S = [[A[i][j] for j in keep] for i in keep]
        for nm, T in named.items():
            ps = list(isos(S, T))
            if not ps: continue
            pats = set()
            for p in ps:
                outv = tuple(A[v][keep[p[i]]] for i in range(6))
                inv = tuple(A[keep[p[i]]][v] for i in range(6))
                pats.add((outv, inv))
            # twin check in the 7-graph
            tw = []
            for u in range(7):
                if u == v: continue
                same_in = all(A[i][u] == A[i][v] for i in range(7) if i not in (u, v))
                same_out = all(A[u][i] == A[v][i] for i in range(7) if i not in (u, v))
                if same_in or same_out:
                    tw.append((u, 'in' if same_in else '', 'out' if same_out else '', 'bi' if A[u][v] and A[v][u] else ('u->v' if A[u][v] else ('v->u' if A[v][u] else 'none'))))
            print(r['d6'], r['edges'], nm, 'v=%d' % v, sorted(pats)[0], 'twins:', tw, 'alpha~%.4f' % tlib.heur_alpha(A))
