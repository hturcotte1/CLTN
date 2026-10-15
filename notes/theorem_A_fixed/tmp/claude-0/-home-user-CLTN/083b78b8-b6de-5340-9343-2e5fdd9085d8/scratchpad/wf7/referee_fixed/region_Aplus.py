"""Test Theorem A+ (R_k SFP for all k >= 0 at every legal point with beta <= 43/20) at many exact
rational legal points, each with the referee's own all-k quotient proof (proof_Rk.py logic) and
ctln at k = 0. Also probe just above 43/20 to see the claim is not vacuous."""
import subprocess, json, sys
from fractions import Fraction as Fr
import random
random.seed(7)
pts = []
# corners: beta = 43/20 near legal boundary t -> 20/3+, and eps tiny
for beta in [Fr(43,20), Fr(21,10), Fr(41,20), Fr(201,100)]:
    tb = 1/(beta-2)
    for t in [tb + Fr(1,1000), tb + Fr(1,10), tb*2, tb*10, Fr(10000)]:
        pts.append((beta, t))
for _ in range(8):
    beta = Fr(2) + Fr(random.randint(1, 150), 1000)
    tb = 1/(beta-2); t = tb + Fr(random.randint(1, 2000), 37)
    pts.append((beta, t))
res = []
for beta, t in pts:
    eps = 1/(1+t); delta = (beta-1)*eps
    assert delta*(1-eps)/eps > 1
    out = f'pt_{beta.numerator}_{beta.denominator}_{t.numerator}_{t.denominator}.json'.replace('/', '_')
    p = subprocess.run(['python3', 'proof_Rk.py', str(eps), str(delta), 'pts/'+out], capture_output=True, text=True)
    v = 'VERDICT True' in p.stdout
    res.append(dict(beta=str(beta), t=str(t), eps=str(eps), delta=str(delta), all_k_SFP=v))
    print(beta, float(t), v, flush=True)
json.dump(res, open('region_Aplus.json', 'w'), indent=1)
