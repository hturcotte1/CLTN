"""Sanity (exact): every n=7 graph whose limit verdict is 'stable_somewhere' is legal-stable somewhere, so by Lemma R
the exact curve test must return 'stable_somewhere' for all of them (124 = 123 non-clique + K7)."""
import json, time
from curve1d import curve_test
from certify_M_violation import from_d6
tally = {}; t0 = time.time(); bad = []
for line in open('/home/user/CLTN/n7/n7_limit_results.jsonl'):
    r = json.loads(line)
    if r['verdict'] != 'stable_somewhere': continue
    v = curve_test(from_d6(r['d6']))['verdict']
    tally[v] = tally.get(v, 0) + 1
    if v != 'stable_somewhere': bad.append(r['d6'])
print(tally, 'bad', bad, 'sec', time.time() - t0)
