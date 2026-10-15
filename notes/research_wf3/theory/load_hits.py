import json, sys
sys.path.insert(0, '/home/user/CLTN')
import graphs, tlib
def n6_hits():
    d = json.load(open('/home/user/CLTN/certificates/n6_summary.json'))
    return {h['idx']: [list(r) for r in h['A']] for h in d['hits']}
def n7_hits():
    out = []
    for line in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl'):
        r = json.loads(line)
        out.append(r)
    return out
