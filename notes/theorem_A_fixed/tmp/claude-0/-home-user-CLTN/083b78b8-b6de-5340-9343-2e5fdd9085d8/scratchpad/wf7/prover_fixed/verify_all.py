"""
verify_all.py -- round-trip check of every instance certificate in PF/ (CLAUDE.md: no claim without a re-runnable certificate).
For every member record in cert_*.json: rebuild the explicit graph from (core, classes, sizes), check it decodes from the stored
d6, re-run ctln.analyze_support at the stored exact (eps, delta), and compare legal / positive / stable / clique / x and the
sha256 of the exact characteristic polynomial.  Members with n > 60 are skipped unless --full is given (n = 107 takes ~13 min).
Also: Geneson's G1 at (1/21, 11/210) (the n = 6 entry of the size-spectrum corollary), and the join instances of
conjP_refutation.json.  Prints ALL CHECKS PASSED or raises.
"""
import json, glob, sys, hashlib
from fractions import Fraction as F
import fixlib as fl

full = '--full' in sys.argv
n_checked = 0
for fn in sorted(glob.glob('cert_*.json')):
    d = json.load(open(fn))
    A = fl.from_d6(d['base']); classes = [tuple(S) for S in d['classes']]
    eps, delta = F(d['eps']), F(d['delta'])
    for m in d['members']:
        G = fl.build(A, classes, m['sizes'])
        assert fl.from_d6(m['d6']) == G, (fn, m['k'])
        if m['n'] > 60 and not full:
            print(fn, 'k', m['k'], 'n', m['n'], 'skipped (use --full)'); continue
        r = fl.ctln.analyze_support(G, list(range(len(G))), eps, delta)
        h = hashlib.sha256(json.dumps(r['charpoly']).encode()).hexdigest()
        assert fl.ctln.legal(eps, delta) and r['positive'] and r['stable'] is True and not r['clique'] and r['stable_fixed_point']
        assert h == m['charpoly_sha256'] and r['x'] == m['x'], (fn, m['k'])
        n_checked += 1
        print(fn, 'k', m['k'], 'n', m['n'], 'OK', flush=True)
# Geneson's G1 at (1/21, 11/210)
G1 = [[0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0]]
r = fl.ctln.analyze_support(G1, list(range(6)), F(1, 21), F(11, 210))
assert r['positive'] and r['stable'] is True and not r['clique'] and fl.ctln.legal(F(1, 21), F(11, 210))
print('G1 (n = 6) at (1/21, 11/210): SFP, not a clique  OK'); n_checked += 1
# join instances
d = json.load(open('conjP_refutation.json'))
R = fl.from_d6('&F^U|m^QTq?')
for inst in d['P3_instances']:
    G = fl.build(R, [(1,) * 7], [inst['k']])
    assert fl.from_d6(inst['d6']) == G
    r = fl.ctln.analyze_support(G, list(range(len(G))), F(1, 21), F(11, 210))
    assert r['positive'] and r['stable'] is True and not r['clique']
    print('join R v K_%d (n = %d) OK' % (inst['k'], len(G))); n_checked += 1
print('ALL CHECKS PASSED (%d instance records re-verified)' % n_checked)
