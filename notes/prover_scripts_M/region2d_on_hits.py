"""EXACT: run the repo's region2d (2-D CAD restricted to the positivity set) on every graph in the given jsonl files
whose curve verdict is 'stable_somewhere' with alpha >= -1/2. An SFP verdict here would be a finite-eps-only
counterexample candidate (not limit-stable anywhere): STOP and hand-verify."""
import sys, json
sys.path.insert(0, '/home/user/CLTN')
import region2d
from certify_M_violation import from_d6
out = []
for f in sys.argv[1:]:
    for l in open(f):
        r = json.loads(l)
        if not r['verdict'].startswith('stable_somewhere') or 'alpha<-1/2' in r['verdict']: continue
        A = from_d6(r['d6'])
        cert = region2d.region2d(A, list(range(len(A))), crosscheck=True)
        rec = {"d6": r['d6'], "curve": r['verdict'], "region2d": cert['verdict'], "positivity_cells": cert.get('positivity_cells'),
               "stable_cells": cert.get('stable_cells')}
        out.append(rec)
        print(('!!! STOP: SFP with alpha >= -1/2  ' if cert['verdict'] != 'never' else '') + json.dumps(rec, default=str)[:400], flush=True)
json.dump(out, open('region2d_on_hits.json', 'w'), indent=1, default=str)
from collections import Counter
print('summary', Counter(r['region2d'] for r in out))
