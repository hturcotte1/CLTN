def my_d6(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    assert s[1] != '~' and 0 < n
    if n >= 63: print('  NOTE: non-standard digraph6 header chr(%d) for n=%d (standard form needs ~ + 18 bits)' % (ord(s[1]), n))
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    assert len(bits) >= n * n and not any(bits[n * n:])
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

