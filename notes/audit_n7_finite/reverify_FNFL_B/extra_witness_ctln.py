"""Extra (outside Lens B): ctln.analyze_support at the three witness points."""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
rows = ["0011110","0011100","1101111","1110111","0011011","0100101","1100000"]
A = [[int(c) for c in r] for r in rows]
for e, d in [("1/49","33/1568"),("1/36","17/588"),("1/101","33/3232")]:
    e, d = F(e), F(d)
    r = ctln.analyze_support(A, list(range(7)), e, d)
    t = (1-e)/e; b = 1 + d/e; q = d*(1-e)/e
    print(e, d, "t", t, "beta", b, "q", q, "legal", q > 1, "stable", r.get('stable'), "|", r.get('stability_detail'),
          "x>0", all(F(xi) > 0 for xi in r["x"]) if r.get('x') else None)
