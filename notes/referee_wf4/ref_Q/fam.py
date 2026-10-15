"""Independent construction of the Q, Q2, Y families (referee ref_Q). Pure integer code; no floats."""
import json
THEORY = '/home/user/CLTN/notes/research_wf3/theory/'

def d6_decode(s, header_len=None):
    assert s[0] == '&'
    if header_len is None:
        if s[1] == '~':
            n = ((ord(s[2]) - 63) << 12) | ((ord(s[3]) - 63) << 6) | (ord(s[4]) - 63); body = s[5:]
        else:
            n = ord(s[1]) - 63; body = s[2:]
    else:
        n = ord(s[1]) - 63; body = s[header_len:]
    bits = []
    for ch in body:
        v = ord(ch) - 63
        assert 0 <= v < 64
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    assert len(bits) >= n * n and not any(bits[n * n:]), 'padding'
    assert len(bits) - n * n < 6
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def d6_encode(A):
    n = len(A); assert n <= 62
    bits = [A[i][j] for i in range(n) for j in range(n)]
    bits += [0] * ((-len(bits)) % 6)
    return '&' + chr(n + 63) + ''.join(chr(63 + int(''.join(map(str, bits[i:i + 6])), 2)) for i in range(0, len(bits), 6))

def rows(rs): return [[int(c) for c in r] for r in rs]

def with_classes(B, classes):
    """classes: list of (S, size). class node v: v -> every base node and every other class node
    (classes are bidirected cliques, mutually bidirected); i -> v iff S[i] == 1."""
    m = len(B); n = m + sum(sz for _, sz in classes)
    A = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): A[i][j] = B[i][j]
    cls = []
    for S, sz in classes:
        cls += [S] * sz
    for a, S in enumerate(cls):
        v = m + a
        for j in range(n):
            if j != v: A[v][j] = 1          # T = all, and class-to-class all edges
        for i in range(m):
            A[i][v] = S[i]
    return A

def is_clique(A):
    n = len(A); return all(A[i][j] and A[j][i] for i in range(n) for j in range(n) if i != j)

BQ = rows(['0111110', '1010011', '1100011', '1110100', '1110011', '0001101', '1111010'])
SQ = (1, 1, 1, 1, 0, 1, 1)
BQ2 = rows(['0111110', '1011011', '1101101', '1100111', '1101001', '0010101', '1100100'])
SQ2a = (1, 1, 0, 1, 1, 1, 1); SQ2b = (1, 1, 1, 0, 1, 1, 1)
BY = rows(['01111101', '10110110', '11011011', '11001111', '11010011', '00101011', '11001001', '01011100'])
SY1 = (0, 1, 1, 1, 1, 1, 1, 1); SY2 = (1, 1, 0, 1, 1, 1, 1, 1)

def Q(k): return with_classes(BQ, [(SQ, k)])
def Q2(k): return with_classes(BQ2, [(SQ2a, k), (SQ2b, 3 * k)])
def Y(k): return with_classes(BY, [(SY1, 4 * k), (SY2, k)])

if __name__ == '__main__':
    # base digraph6 strings named in the notes vs the rows used
    for name, B in [('&F^S{^ReL|?', BQ), ('&F^U|m^QTq?', BQ2), ('&G^ZZZr|KjqTo', BY)]:
        print(name, 'decode==rows:', d6_decode(name) == B, 'encode:', d6_encode(B) == name)
    # stored certificate rows / classes
    cq = json.load(open(THEORY + 'family_Q_certificate.json'))
    c2 = json.load(open(THEORY + 'family_Q2_certificate.json'))
    cy = json.load(open(THEORY + 'family_Y_certificate.json'))
    print('Q rows', rows(cq['B_rows']) == BQ, tuple(cq['S']) == SQ)
    print('Q2 rows', rows(c2['B_rows']) == BQ2, c2['classes'])
    print('Y rows', rows(cy['C8_rows']) == BY, cy['classes'])
    # witnesses
    W = [('Q', 25, Q(25), cq['Q5_witness']['d6'], None),
         ('Q2', 3, Q2(3), c2['W5_witnesses'][0]['d6'], None),
         ('Q2', 10, Q2(10), c2['W5_witnesses'][1]['d6'], None),
         ('Y', 13, Y(13), cy['Y5_witness']['d6'], 2)]
    for fam, k, A, d6, hl in W:
        D = d6_decode(d6, hl)
        print(fam, k, 'n=', len(A), 'decoded n=', len(D), 'equal:', D == A, 'non-clique:', not is_clique(A))
