"""Theorem SJ at further legal points (beta <= 21/10) via proof_SJ.py, a control at beta = 43/20, and ctln instances
of R^{vk} at (1/21, 11/210) for k = 2, 3, 5."""
import subprocess, json
from fractions import Fraction as Fr
from reflib import *
import ctln
res = {'points': [], 'instances': []}
for beta, t in [(Fr(21,10), Fr(10001,1000)), (Fr(21,10), Fr(1000)), (Fr(41,20), Fr(201,10)), (Fr(43,20), Fr(7)), (Fr(43,20), Fr(100))]:
    e = 1/(1+t); d = (beta-1)*e
    p = subprocess.run(['python3', 'proof_SJ.py', str(e), str(d), f'pts/sj_{beta.numerator}_{t.numerator}.json'], capture_output=True, text=True)
    v = 'VERDICT True' in p.stdout
    res['points'].append(dict(beta=str(beta), t=str(t), all_k_SFP=v, tail=p.stdout.strip().splitlines()[-2:]))
    print(beta, t, v, p.stdout.strip().splitlines()[1:2], flush=True)
A0 = d6_decode(R_D6)
for k in [2, 3, 5]:
    A = selfjoin(A0, k); n = len(A)
    r = ctln.analyze_support(A, list(range(n)), Fr(1,21), Fr(11,210))
    res['instances'].append(dict(k=k, n=n, sfp=r['stable_fixed_point'], stable=r['stable'], positive=r['positive'], clique=r['clique'], detail=r['stability_detail']))
    print('SJ ctln', k, n, r['stable_fixed_point'], r['stability_detail'], r['clique'], flush=True)
json.dump(res, open('sj_points.json', 'w'), indent=1)
