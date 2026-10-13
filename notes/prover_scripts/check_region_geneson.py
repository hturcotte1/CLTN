import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import fastregion as fr
A = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
sig = list(range(6))
for t in (fr.INF, F(100), F(60), F(40)):
    c = fr.region(A, sig, t, crosscheck=(t is not fr.INF))
    print("t =", t, "verdict:", c["verdict"], "Delta_degree:", c.get("Delta_degree"))
    for cell in c.get("cells", []):
        print("   cell", cell["lo"], cell["hi"], "sample_beta", cell["sample_beta"], "pos", cell["positive"], "off", cell["off_ok"], "stable", cell["stable"])
    print("   stable q:", fr.stable_q_intervals(c))
    for cr in c.get("critical", []):
        print("   crit", cr["interval"], "deg", len(cr["minpoly_hi_to_lo"]) - 1, cr["vanishing"], cr["minpoly_hi_to_lo"] if len(cr["minpoly_hi_to_lo"]) < 9 else "")
