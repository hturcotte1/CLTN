import tlib, weighted, boundary_tools as bt, sys
from fractions import Fraction as F
A = tlib.from_d6('&FNFL~[udo?')
keep = [0,1,3,4,5,6]
H = [[A[i][j] for j in keep] for i in keep]
for mu in [(1,2,1,1,2,1)] + [(a,1,2,b,1,1) for a in range(1,5) for b in range(1,5)]:
    B = weighted.blowup(H, mu)
    r = bt.analyse(B)
    print(mu, len(B), r, flush=True)
