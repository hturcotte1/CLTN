import itertools, tlib, classext, boundary_tools as bt
from collections import Counter
from fractions import Fraction as F
A = tlib.from_d6('&FNFL~[udo?')
for C in [(0, 1, 1, 1, 0, 0, 0), (1, 0, 1, 1, 1, 0, 0)]:
    T = tuple(1 - c for c in C)
    cnt = Counter()
    for S in itertools.product([0,1], repeat=7):
        B = classext.blow(A, S, T, 1)
        mult, ok = bt.boundary_split(B)
        pos = tlib.positive_at(B, F(2) + F(1, 10**4))
        cnt[(mult, ok, pos)] += 1
    print(C, cnt)
