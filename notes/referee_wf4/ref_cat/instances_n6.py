"""Certify explicit instances k = 1, 2, 50 of the 10 random catalogue entries with ctln.analyze_support
at explicit legal rational parameters eps = 1/(1+t), delta = (beta-1)/(1+t).  beta chosen by the exact
limit test (all roots of Phi_k Re > -1/beta) and exact positivity (w_G(beta) > 0, sigma_S(beta) > 1)."""
import sys, json, random, time
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import ctln
import reprove as R
cat6 = json.load(open(R.CAT + 'families_catalogue_n6.json'))
cat7 = json.load(open(R.CAT + 'families_catalogue_n7.json'))
allc = [('n6', i, e) for i, e in enumerate(cat6)] + [('n7', i, e) for i, e in enumerate(cat7)]
rng = random.Random(20261015)
pick = [allc[0], allc[3], allc[6]]  # one n6 family per hit (G1, H3, G2)
res = []
for tag, i, e in pick:
    A = R.d6_decode(e['base'].split(':')[-1]); S = list(e['S']); m = len(A)
    Phi = R.Phi_W(A, S)
    for kk in (1, 2, 50):
        t0 = time.time()
        rec = {'cat': tag, 'i': i, 'base': e['base'], 'S': S, 'k': kk, 'n': m + kk}
        beta = None
        for d in (10, 20, 40, 100, 200, 1000, 10000):
            b = 2 + F(1, d)
            w = R.w_G(A, b)
            if w is None or not all(v > 0 for v in w) or not sum(w[j] for j in range(m) if S[j]) > 1: continue
            if R.limit_stable(Phi, m, kk, b) is True: beta = b; break
        rec['beta'] = str(beta)
        if beta is None:
            rec['result'] = 'no beta found'; res.append(rec); print(rec, flush=True); continue
        B = R.blow(A, S, kk); n = len(B)
        tmin = int(1 / (beta - 2)) + 1
        found = None
        for mult in (2, 10, 50, 250, 1000):
            t = F(tmin * mult)
            eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
            assert ctln.legal(eps, delta) and delta * (1 - eps) / eps > 1 and t * (beta - 2) > 1
            r = ctln.analyze_support(B, list(range(n)), eps, delta)
            if r['stable_fixed_point'] and r['stable'] is True and r['positive'] and not r['clique'] and not r['target_free']:
                found = (t, eps, delta, r); break
            rec.setdefault('tried', []).append([str(t), r['positive'], r['stable'], r['stability_detail']])
        if found:
            t, eps, delta, r = found
            rec.update(result='SFP', t=str(t), eps=str(eps), delta=str(delta), q=r['q'], stable=r['stable'],
                       detail=r['stability_detail'], positive=r['positive'], clique=r['clique'],
                       x_min=str(min(F(v) for v in r['x'])))
        else:
            rec['result'] = 'no SFP found at tried t'
        rec['secs'] = round(time.time() - t0, 1)
        print({kk2: rec[kk2] for kk2 in rec if kk2 not in ('x_min',)}, flush=True)
        res.append(rec)
json.dump(res, open('instances_n6_results.json', 'w'), indent=1)
print('SFP count', sum(1 for r in res if r['result'] == 'SFP'), 'of', len(res))
