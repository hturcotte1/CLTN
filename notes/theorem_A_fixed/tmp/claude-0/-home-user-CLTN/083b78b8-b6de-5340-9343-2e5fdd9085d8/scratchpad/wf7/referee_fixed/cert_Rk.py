"""Exact ctln certificates for R_k at a single fixed (eps, delta). Usage: python3 cert_Rk.py out.json k1 k2 ..."""
import sys, json, time, hashlib
from fractions import Fraction as Fr
from reflib import *
import ctln
eps, delta = Fr(1,21), Fr(11,210)
out = sys.argv[1]; ks = [int(x) for x in sys.argv[2:]]
recs = []
for k in ks:
    A = build_Rk(k); n = len(A)
    t0 = time.time()
    rec = ctln.analyze_support(A, list(range(n)), eps, delta)
    ok = (ctln.legal(eps, delta) and rec['positive'] and rec['stable'] is True
          and rec['stable_fixed_point'] and not rec['clique'] and rec['off_support_ok'])
    r = dict(k=k, n=n, eps=str(eps), delta=str(delta), legal=ctln.legal(eps, delta),
             positive=rec['positive'], stable=rec['stable'], detail=rec['stability_detail'],
             clique=rec['clique'], sfp=rec['stable_fixed_point'], CERTIFIED=ok,
             A=[''.join(map(str, row)) for row in A], x=rec['x'], charpoly=rec['charpoly'],
             seconds=round(time.time()-t0, 1))
    recs.append(r)
    print(k, n, ok, rec['stability_detail'], r['seconds'], flush=True)
    json.dump(recs, open(out, 'w'), indent=1)
