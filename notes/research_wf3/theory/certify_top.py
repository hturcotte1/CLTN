import json, sys, numpy as np
from fractions import Fraction as F
import tlib, ctln, graphs
def ext(A, pat):
    m = len(A); n = m + 1
    B = [[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    for i in range(m): B[i][m] = pat[i]; B[m][i] = pat[m + i]
    return B
def heur_range(B):
    a = tlib.heur_alpha(B); bc = 1/(1+a)
    ok = [b for b in np.linspace(2.0001, bc - 1e-9, 4000) if tlib.heur_w(B, b).min() > 0]
    return a, bc, (min(ok), max(ok)) if ok else None
def best_rational_beta(B, hi, lo=F(2)):
    """largest rational beta with denominators up to 1000 below hi for which limit_ok_at holds exactly."""
    for den in (1000, 400, 200, 100):
        b = F(int(hi * den), den)
        for _ in range(60):
            if b <= lo: break
            if tlib.limit_ok_at(B, b): return b
            b -= F(1, den)
    return None
if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    seen = set(); out = []
    H7 = {r['d6']: r['A'] for r in __import__('load_hits').n7_hits()}
    H7['&FNFL~[udo?'] = tlib.from_d6('&FNFL~[udo?')
    for q, d6, pat, a in d['best'][:int(sys.argv[2])]:
        B = ext(H7[d6], pat)
        c = graphs.canonical_form(B)
        if c in seen: continue
        seen.add(c)
        a, bc, rng = heur_range(B)
        if not rng: continue
        b = best_rational_beta(B, min(bc, rng[1]))
        wit = tlib.sfp_witness(B, b, ts=(100, 1000, 10**4, 10**5, 10**6)) if b else None
        rec = {'d6': tlib.to_d6(B), 'base': d6, 'pattern': pat, 'heur_alpha': a, 'heur_beta_c': bc, 'heur_pos_range': rng,
               'beta_exact': str(b), 'limit_ok_exact': b is not None}
        if wit:
            t, eps, delta, r = wit
            rec.update({'t': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)), 'q_float': float(ctln.q_of(eps, delta)),
                        'stable': r['stable'], 'positive': r['positive'], 'detail': r['stability_detail'], 'clique': r['clique']})
        out.append(rec)
        print(rec['d6'], 'beta', b, 'q' , rec.get('q_float'), rec.get('detail'), flush=True)
    json.dump(out, open(sys.argv[3], 'w'), indent=1)
