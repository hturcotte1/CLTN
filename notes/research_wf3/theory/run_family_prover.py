"""Run family_prover over bases: the three n=6 hits (all 64 in-sets S) and the 123 n=7 candidates + FNFL (all 128 S,
float prefilter on k in {1,3,10,40} first).  Writes families_catalogue.json."""
import json, itertools, sys, numpy as np
import tlib, classext, family_prover as fp, load_hits
from heur_blowclass import alpha_np
def quick(A, S):
    m = len(A); T = (1,) * m
    Q, P = classext.Phi(A, S, T)
    q = np.array([float(c) for c in Q.all_coeffs()]); p = np.array([float(c) for c in P.all_coeffs()])
    L = max(len(q), len(p)); q = np.concatenate([np.zeros(L - len(q)), q]); p = np.concatenate([np.zeros(L - len(p)), p])
    for kk in (1, 3, 10, 40, 200):
        rts = np.roots(q + kk * p)
        if len(rts) and rts.real.min() <= -0.5: return False
    w = tlib.heur_w(A, 2.0)
    return w.min() > 0 and sum(w[i] for i in range(m) if S[i]) > 1
H6 = load_hits.n6_hits()
bases = [('G1:' + tlib.to_d6(H6[1000255]), H6[1000255]), ('H3:' + tlib.to_d6(H6[1387528]), H6[1387528]), ('G2:' + tlib.to_d6(H6[1424629]), H6[1424629])]
if len(sys.argv) > 1 and sys.argv[1] == 'n7':
    bases = [(r['d6'], r['A']) for r in load_hits.n7_hits()] + [('&FNFL~[udo?', tlib.from_d6('&FNFL~[udo?'))]
out = []
for name, A in bases:
    m = len(A)
    for S in itertools.product([0, 1], repeat=m):
        if not quick(A, S): continue
        r = fp.prove_family(A, S, k0=1)
        r0 = fp.prove_family(A, S, k0=0) if r['ok'] else {'ok': False}
        rec = {'base': name, 'S': S, 'proved_for_all_k>=1': r['ok'], 'proved_for_all_k>=0': r0['ok'], 'detail': r}
        out.append(rec)
        print(name, S, r['ok'], r0['ok'], r.get('chi_L', r.get('why')), flush=True)
json.dump(out, open('families_catalogue_%s.json' % (sys.argv[1] if len(sys.argv) > 1 else 'n6'), 'w'), indent=1)
print('proved families:', sum(1 for o in out if o['proved_for_all_k>=1']))
