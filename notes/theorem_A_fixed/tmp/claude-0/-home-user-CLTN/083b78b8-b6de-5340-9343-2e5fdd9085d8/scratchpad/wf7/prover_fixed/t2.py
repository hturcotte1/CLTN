import fixlib as fl, time
from fractions import Fraction as F
A = fl.from_d6('&F^U|m^QTq?'); S=(1,1,0,1,1,1,1)
for kk in (1, 20, 40):
    t0=time.time()
    rec = fl.certify_member(A,[S],[kk],F(21,10),20)
    print(kk, rec['n'], rec['legal'], rec['positive'], rec['stable'], rec['stability_detail'], rec['clique'], round(time.time()-t0,1), flush=True)
