import fixlib as fl, sympy as sp, time
from fractions import Fraction as F
G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
S = (1,1,1,1,1,0)
t0=time.time()
print('identity F', [fl.check_identity(G1,[S],[kk]) for kk in (1,2,3)], time.time()-t0)
B = fl.from_d6('&F^S{^ReL|?'); SQ=(1,1,1,1,0,1,1)
print('identity Q', [fl.check_identity(B,[SQ],[kk]) for kk in (1,2,3)], time.time()-t0)
# F at beta=21/10, t=20
for kv in (37,38,39,40,41,45,60):
    print(kv, fl.stable_at_k(G1,[S],[1],F(21,10),20,kv))
k0, mins = fl.k_threshold(G1,[S],[1],F(21,10),20)
print('F k0', k0)
for mn in mins: print(sp.Poly(mn, fl.k).LC() if sp.Poly(mn,fl.k).degree()>0 else mn)
