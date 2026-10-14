"""Exact certification of FLAG_POS / flag_lit records from scan2.py.
For each flagged (graph, beta): choose a rational beta (limit_denominator) at which the float pattern
persists, a rational t1 in a stable interval and a rational t2 > t1 in an unstable interval, then decide
both points EXACTLY with ctln.analyze_support (eps = 1/(t+1), delta = (beta-1) eps), and decide
alpha(A|perp) vs -(1-1/beta) exactly (alpha_exact.py, sympy root isolation).
A None verdict is never counted.  Output: one JSON line per record.
usage: certify.py <scan.jsonl> <out.jsonl> [kinds=FLAG_POS] [max_records]"""
import sys, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, '/home/user/CLTN')
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf/lemmaM_refute')
import ctln
from scan import from_d6
from scan2 import analyze_batch
from alpha_exact import alpha_vs

def exact_point(A, beta, t):
    eps = 1 / (t + 1); delta = (beta - 1) * eps
    r = ctln.analyze_support(A, list(range(len(A))), eps, delta)
    return dict(beta=str(beta), t=str(t), eps=str(eps), delta=str(delta), q=str(ctln.q_of(eps, delta)),
                legal=bool(ctln.q_of(eps, delta) > 1), stable=r['stable'], detail=r['stability_detail'],
                positive=r['positive'], stable_fixed_point=r['stable_fixed_point'], x=r['x'])

def try_beta(A, beta_f, need_pos):
    An = np.array([A], float)
    for den in (8, 16, 25, 50, 100, 200, 500, 1000, 4000, 20000, 10 ** 5, 10 ** 6):
        b = F(beta_f).limit_denominator(den)
        if b <= 2: continue
        R = analyze_batch(An, np.array([[float(b)]]))
        v = R['valid'][0, 0]
        samp = R['samp'][0, 0][v]; st = R['stab'][0, 0][v]; un = R['unst'][0, 0][v]
        pos = bool(R['pos'][0, 0, 0])
        if need_pos and not pos: continue
        # stable finite sample followed by unstable sample
        for i in range(len(samp)):
            if st[i] and np.isfinite(samp[i]):
                for j in range(i + 1, len(samp)):
                    if un[j]:
                        yield b, samp[i], (samp[j] if np.isfinite(samp[j]) else None), R['bounds'][0, 0]
                        break
                break

def rat_t(x, lo_bound):
    for den in (1, 2, 4, 10, 100, 1000, 10 ** 4, 10 ** 6):
        t = F(x).limit_denominator(den)
        if abs(float(t) - x) <= 1e-3 * x and t > lo_bound: return t
    return F(x)

def certify(A, beta_f, need_pos=True):
    for b, s1, s2, bounds in try_beta(A, beta_f, need_pos):
        tmin = 1 / (b - 2)
        t1 = rat_t(s1, tmin)
        if s2 is None:
            fin = [x for x in bounds if np.isfinite(x)]
            s2 = max(fin) * 4 + 100
        t2 = rat_t(s2, t1)
        p1 = exact_point(A, b, t1); p2 = exact_point(A, b, t2)
        if p1['stable'] is True and p2['stable'] is False and p1['legal'] and p2['legal'] and t2 > t1:
            P, al = alpha_vs(A, str(-(1 - 1 / b)))
            sides = [o['side'] for o in al]
            alpha_side = '>' if '>' in sides else ('<' if all(s == '<' for s in sides) else '?')
            return dict(ok=True, beta=str(b), p_stable=p1, p_unstable=p2,
                        alpha_vs_threshold=alpha_side, threshold=str(-(1 - 1 / b)),
                        chi_perp_times_m=str(P.as_expr()),
                        alpha_roots=[o for o in al if o['side'] != '<'])
        last = (p1, p2)
    return dict(ok=False)

if __name__ == '__main__':
    inp, outp = sys.argv[1], sys.argv[2]
    kinds = sys.argv[3].split(',') if len(sys.argv) > 3 else ['FLAG_POS']
    mx = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
    out = open(outp, 'w'); n = 0
    for l in open(inp):
        r = json.loads(l)
        if r.get('kind') not in kinds: continue
        n += 1
        if n > mx: break
        A = from_d6(r['d6'])
        need_pos = r['kind'] == 'FLAG_POS'
        hits = [h for h in r['hits'] if (h['a'] or h['b']) and (h['pos'] or not need_pos)]
        # try the hit with the biggest stable margin first
        hits.sort(key=lambda h: -max([m for m, s in zip(h['minre'], h['stab']) if s] or [0]))
        res = dict(ok=False)
        for h in hits[:4]:
            res = certify(A, h['beta'], need_pos)
            if res['ok']: break
        rec = dict(d6=r['d6'], src=r['src'], kind=r['kind'], adjacency_rows=[''.join(map(str, row)) for row in A],
                   alpha_float=r['alpha'], **res)
        out.write(json.dumps(rec) + '\n'); out.flush()
        print(r['d6'], res.get('ok'), res.get('beta'), res.get('p_stable', {}).get('t'), res.get('p_unstable', {}).get('t'),
              res.get('alpha_vs_threshold'), flush=True)
