"""HEURISTIC (floats only steer): min Re of the limit polynomial chi_L for every catalogue family; keep > -1/2."""
import json, sympy as sp, numpy as np
TH='/home/user/CLTN/notes/research_wf3/theory/'
z=sp.symbols('z')
out=[]
for fn in ['families_catalogue_n6.json','families_catalogue_n7.json']:
    for rec in json.load(open(TH+fn)):
        chiL = sp.sympify(rec['detail']['chi_L'])
        P = sp.Poly(chiL, z)
        rts = np.roots([float(c) for c in P.all_coeffs()])
        nz = [r for r in rts if abs(r) > 1e-9]
        mr = min(r.real for r in nz) if nz else 0.0
        base = rec['base'].split(':')[-1]
        out.append({'base': base, 'S': rec['S'], 'chi_L': rec['detail']['chi_L'], 'minRe': mr, 'k>=0': rec['proved_for_all_k>=0']})
good=[o for o in out if o['minRe'] > -0.5 + 1e-9]
print(len(out), 'families;', len(good), 'with min Re chi_L > -1/2 (float)')
good.sort(key=lambda o: o['minRe'], reverse=True)
for o in good[:15]: print(o['base'], o['S'], round(o['minRe'],4), o['chi_L'])
json.dump(good, open('heur_catalogue_good.json','w'), indent=0)
