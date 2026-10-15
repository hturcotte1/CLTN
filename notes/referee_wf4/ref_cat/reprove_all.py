"""Run reprove.reprove on EVERY catalogue entry (n6 + n7), k0 = 1, and k0 = 0 where the catalogue claims k >= 0."""
import json, time, sys
import reprove as R
cat6 = json.load(open(R.CAT + 'families_catalogue_n6.json'))
cat7 = json.load(open(R.CAT + 'families_catalogue_n7.json'))
out = []; bad = []
for tag, cat in (('n6', cat6), ('n7', cat7)):
    for i, e in enumerate(cat):
        r, A, S, Phi = R.reprove(e)
        r0 = R.reprove(e, k0=0)[0] if e['proved_for_all_k>=0'] else None
        rec = {'cat': tag, 'i': i, 'base': e['base'], 'S': e['S'], 'proved_k>=1': r['PROVED'],
               'phi_match': r['Phi_matches_catalogue'], 'minors_match': r['minors_match_catalogue'],
               'proved_k>=0': (r0['PROVED'] if r0 else None), 'claims_k>=0': e['proved_for_all_k>=0'],
               'chi_L_limit': None}
        if not (r['PROVED'] and r['Phi_matches_catalogue'] and r['minors_match_catalogue'] and (r0 is None or r0['PROVED'])):
            bad.append(rec)
        out.append(rec)
        if i % 50 == 0: print(tag, i, len(bad), flush=True)
json.dump({'records': out, 'bad': bad}, open('reprove_all_results.json', 'w'), indent=0)
print('total', len(out), 'bad', len(bad))
