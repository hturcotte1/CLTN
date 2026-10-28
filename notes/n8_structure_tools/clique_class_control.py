"""Positive control for clique_class_test.strict: every strict survivor of the production n = 7 and n = 8 maps
(n7/n7_limit_results.jsonl all 1,889; n8/n8_limit_results.jsonl every 25th line + all 9,133 hits) must pass,
and their maximum bidirected clique (excluding K_n) is reported.  Output: clique_class_control.json"""
import json, os
from collections import Counter
import lib
from clique_class_test import strict
from family_test import max_clique
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for n, path, step in ((7, '/n7/n7_limit_results.jsonl', 1), (8, '/n8/n8_limit_results.jsonl', 25)):
    c = Counter(); om = Counter()
    for k, l in enumerate(open(lib.REPO + path)):
        r = json.loads(l)
        if k % step and r['verdict'] != 'stable_somewhere':
            continue
        A = lib.from_d6(r['d6'])
        c['tested'] += 1; c['strict_pass'] += strict(A)
        if not r.get('complete'):
            om[max_clique(A)] += 1
    out[f'n{n}'] = dict(c, omega_noncomplete=dict(sorted(om.items())))
json.dump(out, open(os.path.join(HERE, 'clique_class_control.json'), 'w'), indent=1)
print(out)
