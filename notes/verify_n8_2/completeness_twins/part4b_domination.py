"""Part 4 follow-up: the 232 switches absent from the results file -- own domination check (filter8's F1 definition) and
source reasons; and own domination check that every present switch is domination-free."""
import sys, json
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from part2b_spectral import domination_free, classify_limit
from collections import Counter
R = '/home/user/CLTN/'; D = R + 'notes/verify_n8_2/completeness_twins/'
res = {}
for l in open(R + 'n8/n8_limit_results.jsonl'):
    r = json.loads(l); res[r['d6']] = r
rows = [json.loads(l) for l in open(D + 'part4_never_closure_rows.jsonl')]
c = Counter()
for x in rows:
    B = d6_decode(x['switch_d6']); df = domination_free(B); cl = classify_limit(B)
    src = res[x['source_d6']]; reason = src.get('reason', 'never_full')
    c[(x['in_results'], df, cl, reason)] += 1
print({str(k): v for k, v in c.items()})
json.dump({'(in_results, own_domination_free, own_limit_class, source_reason) -> count': {str(k): v for k, v in c.items()}}, open(D + 'part4b_summary.json', 'w'), indent=1)
