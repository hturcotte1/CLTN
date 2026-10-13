"""Undecided (stable is None) rate for n=6 full support at several legal points, same 3000 random graphs.
Also: how many undecided cases already FAIL positivity (so 'not a stable FP support' is decided by (1))."""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
import ctln
rng = random.Random(7)
n = 6
graphs_ = [[[0 if i == j else rng.randint(0, 1) for j in range(n)] for i in range(n)] for _ in range(3000)]
pts = [(F(1,2), F(3,2)), (F(1,2), F(3,2) + F(1,997)), (F(500,1001), F(1501,1000)), (F(1,10), F(19,90)*F(10,1)), (F(1,100), F(29,2500))]
for eps, delta in pts:
    und = und_pos_false = 0
    for A in graphs_:
        r = ctln.analyze_support(A, list(range(n)), eps, delta)
        if r["stable"] is None:
            und += 1
            if r["x"] is not None and not r["positive"]:
                und_pos_false += 1
    print(f"eps={eps} delta={delta} q={ctln.q_of(eps,delta)} beta={1+delta/eps} t={(1-eps)/eps}: undecided {und}/3000, of which positivity already fails: {und_pos_false}", flush=True)
