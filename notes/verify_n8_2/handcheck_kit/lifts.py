import sys, json, subprocess
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/handcheck_kit')
from engine import d6_decode, d6_encode
LABELG = '/home/user/CLTN/n7/labelg'
HITS = '/home/user/CLTN/n8/n8_limit_HITS.jsonl'
N7 = '/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl'

def labelg(lines):
    p = subprocess.run([LABELG, '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = p.stdout.split()
    assert len(out) == len(lines), (len(out), len(lines), p.stderr[:300])
    return out

def blowup(A, v):
    n = len(A)
    B = [row[:] + [row[v]] for row in A]          # new column: arcs into v' copy arcs into v
    B.append(A[v][:] + [0])                        # new row: arcs out of v' copy arcs out of v
    B[v][n] = 1; B[n][v] = 1                       # v <-> v'
    B[n][n] = 0
    return B

hits = [json.loads(l) for l in open(HITS)]
hit_d6 = [h['d6'] for h in hits]
hit_can = labelg(hit_d6)
assert len(set(hit_can)) == len(hit_can), 'hits not pairwise non-isomorphic?!'
can2idx = {c: i for i, c in enumerate(hit_can)}

# --- Geneson double blow-ups
Ag = d6_decode('&EZjtwWE')
found = []
dbl = {}
for v in range(6):
    A7 = blowup(Ag, v)
    for w in range(7):
        A8 = blowup(A7, w)
        dbl[(v, w)] = d6_encode(A8)
keys = list(dbl)
cans = labelg([dbl[k] for k in keys])
for k, c in zip(keys, cans):
    if c in can2idx:
        found.append((k, dbl[k], can2idx[c], hit_d6[can2idx[c]], hits[can2idx[c]]['edges']))
print('Geneson double blow-ups in hits:', found)
print('distinct double blow-ups:', len(set(cans)))
# also single blow-ups to 7 nodes vs n7 list
n7 = [json.loads(l) for l in open(N7)]
n7_can = labelg([r['d6'] for r in n7])
assert len(set(n7_can)) == 123
sing = labelg([d6_encode(blowup(Ag, v)) for v in range(6)])
print('Geneson single blow-ups among the 123:', [(v, n7[n7_can.index(c)]['d6']) for v, c in enumerate(sing) if c in n7_can])

# --- 123 -> 8 lift relation
rows = []
in_hits = 0
for i, r in enumerate(n7):
    A7 = d6_decode(r['d6'])
    assert len(A7) == 7
    bl = [d6_encode(blowup(A7, v)) for v in range(7)]
    cs = labelg(bl)
    lift = []
    for v, c in enumerate(cs):
        if c in can2idx:
            j = can2idx[c]
            lift.append({'v': v, 'hit_index': j, 'hit_d6': hit_d6[j], 'hit_arcs': hits[j]['edges']})
    in_hits += len({l['hit_index'] for l in lift})
    rows.append({'n7_index': i, 'n7_d6': r['d6'], 'n7_arcs': r['edges'], 'lifts_in_hits': lift})
json.dump(rows, open('lift_123_to_8.json', 'w'), indent=1)
hit_idx_from_lift = sorted({l['hit_index'] for r in rows for l in r['lifts_in_hits']})
print('parents with >=1 lift in hits:', sum(1 for r in rows if r['lifts_in_hits']), 'of 123')
print('distinct hits that are closed-twin blow-ups of a 123-parent:', len(hit_idx_from_lift))
json.dump({'hit_indices': hit_idx_from_lift}, open('lift_hit_indices.json', 'w'))
