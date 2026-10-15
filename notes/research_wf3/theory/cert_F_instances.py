"""Exact SFP witnesses for the family F_k (family_F.py) via ctln.analyze_support (the repo verifier).
Writes certificates/F_family_witnesses.json (in this scratch dir)."""
import json, sys, time
from fractions import Fraction as F
import tlib, ctln
from family_F import F_k
out = []
for k, beta in [(1, F(2)+F(1,10)), (2, F(58,25)), (3, F(58,25)), (4, F(58,25)), (2, F(21,10)), (3, F(21,10)), (4, F(21,10)), (6, F(21,10)), (10, F(21,10))]:
    A = F_k(k); n = len(A)
    t0 = time.time()
    lim = tlib.limit_ok_at(A, beta)
    wit = tlib.sfp_witness(A, beta, ts=(20, 50, 100, 200, 500, 1000, 3000, 10**4, 10**5)) if lim else None
    rec = {'family': 'F_k (Geneson G1 + bidirected clique class of size k; class <- G1 nodes 0..4, class -> all G1)',
           'k': k, 'n': n, 'd6': tlib.to_d6(A), 'A_rows': [''.join(map(str, r)) for r in A], 'beta': str(beta),
           'limit_ok_at_beta_exact': lim}
    if wit:
        t, eps, delta, r = wit
        rec.update({'t': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)),
                    'legal': ctln.legal(eps, delta), 'positive': r['positive'], 'stable': r['stable'],
                    'stability_detail': r['stability_detail'], 'stable_fixed_point': r['stable_fixed_point'],
                    'clique': r['clique'], 'target_free': r['target_free'], 'x': r['x'], 'charpoly': r['charpoly']})
    out.append(rec)
    print(k, n, str(beta), lim, None if not wit else (str(wit[0]), str(ctln.q_of(wit[1], wit[2])), wit[3]['stability_detail']), '%.1fs' % (time.time() - t0), flush=True)
json.dump(out, open('F_family_witnesses.json', 'w'), indent=1)
