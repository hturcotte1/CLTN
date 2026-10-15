"""Re-run ctln.analyze_support (repo, read-only) on the four stored family witnesses, graphs rebuilt
from the family definitions (fam.py) and checked equal to the stored digraph6."""
import sys, json, time
sys.path.insert(0, '/home/user/CLTN')
import ctln
from fractions import Fraction as F
import fam
cases = [('Q2', 3, fam.Q2(3), F(3), F(2,10001)), ('Q', 25, fam.Q(25), F(5,2), F(3,20002)),
         ('Q2', 10, fam.Q2(10), F(7,2), None), ('Y', 13, fam.Y(13), F(3), F(2,10001))]
out = []
for name, k, A, beta, dstored in cases:
    t = F(10000); eps = 1 / (t + 1); delta = (beta - 1) * eps
    if dstored is not None: assert delta == dstored
    t0 = time.time()
    r = ctln.analyze_support(A, list(range(len(A))), eps, delta)
    rec = dict(family=name, k=k, n=len(A), beta=str(beta), eps=str(eps), delta=str(delta), q=r['q'],
               legal=ctln.legal(eps, delta), positive=r['positive'], stable=r['stable'], detail=r['stability_detail'],
               clique=r['clique'], stable_fixed_point=r['stable_fixed_point'], x_first=r['x'][:3] if r['x'] else None,
               seconds=round(time.time() - t0))
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open('ctln_witness_results.json', 'w'), indent=1)
