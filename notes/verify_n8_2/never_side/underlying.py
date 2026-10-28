"""Map every S digraph of n8_limit_results.jsonl to the graph6 of its underlying simple graph (own encoder)."""
import json, sys
def d6_to_A(s):
    assert s[0] == '&'; n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; assert 0 <= v < 64
        bits += [(v >> (5 - k)) & 1 for k in range(6)]
    assert len(bits) >= n * n
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]
def g6(adj):
    n = len(adj); bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(adj[i][j])
    while len(bits) % 6: bits.append(0)
    out = chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = v * 2 + b
        out += chr(63 + v)
    return out
for line in open(sys.argv[1]):
    r = json.loads(line); A = d6_to_A(r['d6']); n = len(A)
    U = [[1 if (A[i][j] or A[j][i]) else 0 for j in range(n)] for i in range(n)]
    print(g6(U), r['d6'], r['verdict'], 'cert' in r)
