"""graph6 / digraph6 helpers for the K8 bidirected-subgraph enumeration (pure Python, integers only).
Conventions (nauty formats.txt): graph6 = N(n) then the upper triangle column by column, bit (i,j) for
j = 1..n-1, i = 0..j-1; digraph6 = '&' N(n) then the n*n adjacency bits row-major, bit (i,j) = arc i->j;
six bits per character, value + 63, zero-padded."""
def _bits_to_chars(bits):
    bits = bits + [0] * ((6 - len(bits) % 6) % 6)
    return ''.join(chr(63 + int(''.join(map(str, bits[k:k + 6])), 2)) for k in range(0, len(bits), 6))
def _chars_to_bits(s, nbits):
    out = []
    for ch in s:
        v = ord(ch) - 63
        assert 0 <= v < 64, 'bad graph6 character'
        out += [(v >> (5 - b)) & 1 for b in range(6)]
    assert len(out) >= nbits, 'string too short'
    return out[:nbits]
def decode_g6(s):
    n = ord(s[0]) - 63; assert 1 <= n < 63
    bits = _chars_to_bits(s[1:], n * (n - 1) // 2)
    adj = [[0] * n for _ in range(n)]; k = 0
    for j in range(1, n):
        for i in range(j):
            adj[i][j] = adj[j][i] = bits[k]; k += 1
    return n, adj
def encode_g6(n, adj):
    bits = [adj[i][j] for j in range(1, n) for i in range(j)]
    return chr(63 + n) + _bits_to_chars(bits)
def decode_d6(s):
    assert s[0] == '&'; n = ord(s[1]) - 63; assert 1 <= n < 63
    bits = _chars_to_bits(s[2:], n * n)
    return n, [bits[i * n:(i + 1) * n] for i in range(n)]
def encode_d6(n, adj):
    return '&' + chr(63 + n) + _bits_to_chars([adj[i][j] for i in range(n) for j in range(n)])
def complement_g6(s):
    n, adj = decode_g6(s)
    return encode_g6(n, [[0 if i == j else 1 - adj[i][j] for j in range(n)] for i in range(n)])
def edges_g6(s):
    n, adj = decode_g6(s)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i][j]]
def complete_g6(n):
    return encode_g6(n, [[0 if i == j else 1 for j in range(n)] for i in range(n)])
def add_bidirected(d6, g6):
    n, D = decode_d6(d6); m, B = decode_g6(g6); assert n == m
    return encode_d6(n, [[1 if (D[i][j] or B[i][j]) else 0 for j in range(n)] for i in range(n)])
