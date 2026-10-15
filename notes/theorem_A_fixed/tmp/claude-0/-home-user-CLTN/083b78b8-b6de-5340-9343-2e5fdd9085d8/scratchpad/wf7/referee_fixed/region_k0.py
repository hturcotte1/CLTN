"""k = 0 (R itself) at every region_Aplus point via ctln; plus one control above 43/20."""
import json
from fractions import Fraction as Fr
from reflib import *
import ctln
pts = json.load(open('region_Aplus.json'))
A = d6_decode(R_D6); out = []
for p in pts:
    e, d = Fr(p['eps']), Fr(p['delta'])
    r = ctln.analyze_support(A, list(range(7)), e, d)
    ok = ctln.legal(e, d) and r['stable_fixed_point'] and r['stable'] is True
    out.append(dict(beta=p['beta'], t=p['t'], k0_sfp=ok, all_k_ge_1=p['all_k_SFP']))
    print(p['beta'], p['t'], ok)
# control: beta = 11/5 near legal boundary t = 5 + 1/100: is some small k unstable? (claim A+ should not extend)
beta, t = Fr(11, 5), Fr(501, 100); e = 1/(1+t); d = (beta-1)*e
ctrl = []
for k in range(0, 6):
    Ak = build_Rk(k); r = ctln.analyze_support(Ak, list(range(len(Ak))), e, d)
    ctrl.append(dict(k=k, sfp=r['stable_fixed_point'], detail=r['stability_detail']))
    print('control beta=11/5 t=5.01 k', k, r['stable_fixed_point'], r['stability_detail'])
json.dump(dict(points=out, control=ctrl), open('region_k0.json', 'w'), indent=1)
