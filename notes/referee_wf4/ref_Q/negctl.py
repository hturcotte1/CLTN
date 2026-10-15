"""Negative controls: one below each witnessed threshold, same beta; both the independent checks and ctln."""
import sys, json
sys.path.insert(0, '/home/user/CLTN')
import ctln, inst, sym
from fractions import Fraction as F
Phis = {nm: sym.family_Phi(sym.FAMS[nm]['B'], sym.FAMS[nm]['classes'])[0] for nm in ('Q', 'Q2', 'Y')}
for nm, kk, beta, ppk in [('Q', 24, F(5, 2), 1), ('Q2', 2, F(3), 4), ('Q2', 9, F(7, 2), 4), ('Y', 12, F(3), 5)]:
    A = sym.FAMS[nm]['build'](kk)
    rec = inst.run(nm, kk, A, beta, Phis[nm], ppk * kk)
    eps = F(1, 10001); r = ctln.analyze_support(A, list(range(len(A))), eps, (beta - 1) * eps)
    rec['ctln_positive'] = r['positive']; rec['ctln_stable'] = r['stable']; rec['ctln_detail'] = r['stability_detail']
    print(json.dumps(rec), flush=True)
