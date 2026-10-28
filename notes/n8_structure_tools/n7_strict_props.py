"""The same graph properties on ALL 1,889 strict n = 7 survivors (n7/n7_limit_results.jsonl), by verdict.
Output: n7_strict_props.json"""
import json, os
from collections import Counter
import lib
from family_test import checks
from more_universals import ham
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for l in open(lib.REPO + '/n7/n7_limit_results.jsonl'):
    r = json.loads(l)
    if r.get('complete'):
        continue
    k = 'hit' if r['verdict'] == 'stable_somewhere' else ('never_positivity' if 'a_i' in r.get('reason', '') else 'never_full')
    A = lib.from_d6(r['d6'])
    c = checks(A); c['hamiltonian'] = ham(A)
    o = out.setdefault(k, Counter()); o['total'] += 1
    for kk in ('omega_le_n_minus_3', 'strongly_connected', 'induced_C3', 'recv_all_implies_send_all', 'hamiltonian'):
        o[kk] += bool(c[kk])
json.dump({k: dict(v) for k, v in out.items()}, open(os.path.join(HERE, 'n7_strict_props.json'), 'w'), indent=1)
print(out)
