"""The n = 6 and n = 7 counterexample lists used for comparison (all from repository certificate files).
n6: the three signed-off limit counterexamples (certificates/region_n6_limit_shard0of1_HITS.jsonl keys,
    decoded with enum_cache.key_from_int + graphs.from_key);
n7: the 123 limit counterexamples (certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl) and the
    finite-eps-only &FNFL~[udo? (certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json)."""
import sys, json
sys.path.insert(0, '/home/user/CLTN')
import lib

def n6_list():
    import enum_cache, graphs
    out = []
    for l in open(lib.REPO + '/certificates/region_n6_limit_shard0of1_HITS.jsonl'):
        r = json.loads(l)
        if r['complete']:
            continue
        key = enum_cache.key_from_int(r['key'], 6)
        A = [list(row) for row in graphs.from_key(key, 6)]
        out.append((r['idx'], lib.to_d6(A)))
    return out

def n7_list():
    return [json.loads(l)['d6'] for l in open(lib.REPO + '/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]

FNFL = '&FNFL~[udo?'
