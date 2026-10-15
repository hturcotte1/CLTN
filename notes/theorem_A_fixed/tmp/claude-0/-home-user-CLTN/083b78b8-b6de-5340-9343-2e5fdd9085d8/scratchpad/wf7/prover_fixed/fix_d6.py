"""
fix_d6.py -- the first version of fixlib.to_d6 wrote a one-byte digraph6 header even for n >= 63 (wrong; nauty uses '~' + 3 bytes).
Only the stored 'd6' strings of members with n >= 63 were affected; every exact computation used the explicit matrix built by
fixlib.build(core, classes, sizes).  This script recomputes the d6 field of every member record in the cert_*.json files from
(core, classes, sizes) with the corrected encoder, checks the round trip from_d6(to_d6(G)) == G, and rewrites the files.
"""
import json, glob
import fixlib as fl

for fn in sorted(glob.glob('cert_*.json')):
    d = json.load(open(fn))
    A = fl.from_d6(d['base']); classes = [tuple(S) for S in d['classes']]
    changed = 0
    for m in d['members']:
        G = fl.build(A, classes, m['sizes'])
        s = fl.to_d6(G)
        assert fl.from_d6(s) == G and len(G) == m['n']
        if m['d6'] != s:
            m['d6'] = s; changed += 1
    json.dump(d, open(fn, 'w'), indent=1)
    print(fn, 'members', len(d['members']), 'd6 fields corrected', changed)
