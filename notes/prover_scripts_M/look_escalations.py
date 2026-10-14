import sys, numpy as np
sys.path.insert(0, '/home/user/CLTN')
from curve1d import curve_test, coeff_polys
from fastregion import hurwitz_det_poly
from run_scan_small import iso
from scan_graph_level import nu_inf
As = list(iso(5))
for A in As:
    A = A.astype(int).tolist()
    r = curve_test(A)
    if r['verdict'] != 'escalate': continue
    P = coeff_polys(A)
    zero_c = [k for k, p in enumerate(P) if not p]
    zero_D = [k for k in range(1, 6) if not hurwitz_det_poly(P, k)]
    print(A, 'zero coeff idx', zero_c, 'identically-zero Hurwitz minors', zero_D, 'nu_inf', round(float(nu_inf(np.array([A], float))[0]), 4))
