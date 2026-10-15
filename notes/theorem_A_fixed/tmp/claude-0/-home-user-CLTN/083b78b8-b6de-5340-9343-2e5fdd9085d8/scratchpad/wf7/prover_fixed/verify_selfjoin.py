"""
verify_selfjoin.py -- round-trip check of cert_selfjoin_R_eps1_21.json: rebuild each R^{vk} from the core, compare with the stored
digraph6, re-run ctln.analyze_support at the stored exact (eps, delta), compare positive / stable / clique / x / sha256 of the exact
characteristic polynomial.  Members with n > 40 are skipped unless --full.  Prints ALL CHECKS PASSED or raises.
"""
import json, sys, hashlib
from fractions import Fraction as F
import fixlib as fl
from selfjoin import selfjoin

full = '--full' in sys.argv
d = json.load(open('cert_selfjoin_R_eps1_21.json'))
A = fl.from_d6(d['base']); eps, delta = F(d['eps']), F(d['delta'])
assert fl.ctln.legal(eps, delta)
n_ok = 0
for m in d['members']:
    G = selfjoin(A, m['k'])
    assert fl.from_d6(m['d6']) == G and len(G) == m['n']
    if m['n'] > 40 and not full:
        print('k', m['k'], 'n', m['n'], 'skipped (use --full)'); continue
    r = fl.ctln.analyze_support(G, list(range(len(G))), eps, delta)
    h = hashlib.sha256(json.dumps(r['charpoly']).encode()).hexdigest()
    assert r['positive'] and r['stable'] is True and not r['clique'] and r['stable_fixed_point']
    assert h == m['charpoly_sha256'] and r['x'] == m['x']
    n_ok += 1; print('k', m['k'], 'n', m['n'], 'OK', flush=True)
print('ALL CHECKS PASSED (%d self-join records re-verified)' % n_ok)
