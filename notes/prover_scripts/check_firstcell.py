"""Probe sympy's isolator at the inf boundary, and fastregion._cells."""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ
import fastregion as fr

def iso(polys_lo_hi, bmin):
    return dup_isolate_real_roots_list([list(reversed(p)) for p in polys_lo_hi], ZZ, inf=bmin, strict=True, basis=True)

bmin = F(7, 3)
cases = {
  "root exactly at b_min (3b-7)(b-5)": [[-35, 22, -3]],            # (3b-7)(b-5) = 3b^2-22b+35 -> lo->hi [35,-22,3]
  "root 1e-9 above b_min": [[-(7*10**9+1), 3*10**9]],               # 3e9 b - (7e9+1)
  "root 1e-9 below b_min": [[-(7*10**9-1), 3*10**9]],
  "irrational root just above (b^2 - 5.44444445)": [[-490000001, 0, 90000000]],  # b = sqrt(490000001/9e7) ~ 2.33333334
  "irrational root just below": [[-489999999, 0, 90000000]],
  "two roots straddling region": [[-35, 22, -3], [-(7*10**9+1), 3*10**9]],
}
cases["root exactly at b_min (3b-7)(b-5)"] = [[35, -22, 3]]
for name, ps in cases.items():
    out = iso(ps, bmin)
    print(name)
    for (a, b), idx, f in out:
        print("   interval", (str(a), str(b)), "a>=bmin:", a >= bmin, "root-in-open:", a != b, "factor", f)
    out2 = [iv for iv in out if not (iv[0][0] == iv[0][1] == bmin)]
    cells = fr._cells(out2, bmin)
    print("   cells:", [(str(lo), None if hi is None else str(hi), str(s)) for _, _, lo, hi, s in cells])
