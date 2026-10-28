import json
recs = [json.loads(l) for l in open('certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
from collections import Counter
c = Counter(); bn = []
for r in recs:
    A = r['A']; n = len(A); arcs = sum(map(sum, A))
    B = sum(1 for i in range(n) for j in range(i+1, n) if A[i][j] and A[j][i])
    Nn = sum(1 for i in range(n) for j in range(i+1, n) if not A[i][j] and not A[j][i])
    O = n*(n-1)//2 - B - Nn
    assert arcs == 2*B + O and r['edges'] == arcs
    c[arcs] += 1; bn.append((B, Nn))
print('n=7 candidates', len(recs), 'arc counts', sorted(c.items()))
print('threshold m(m-1)/2+1 =', 7*6//2 + 1, '; min arcs observed', min(c), '; all B>Nn:', all(B > Nn for B, Nn in bn), '; min B-Nn', min(B-Nn for B,Nn in bn))
print('n=8 threshold', 8*7//2 + 1, '; underlying edges >= ', (8*7//2)//2 + 1)
