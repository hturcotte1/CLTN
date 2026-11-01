"""Recount the 'primitive' counterexamples (no smaller counterexample as an induced subgraph)
at n = 6, 7, 8 from the certificate files, exactly, and compute the three shares as Fractions.

Inputs (all under certificates/):
  n6_summary.json                               -> the 3 non-complete limit-stable six-node graphs (A matrices)
  COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl      -> the 123 seven-node limit counterexamples (d6)
  COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json  -> the finite-eps-only seven-node graph &FNFL~[udo? (d6)
  n8_counterexamples_classes.json               -> the 9,133 eight-node limit counterexamples (members_d6 per class)

'Primitive' as used in notes/n8_structure.md and notes/REPORT_n8.md: at n = 7, contains no six-node counterexample
as an induced subgraph; at n = 8, contains no six-node and no seven-node LIMIT counterexample as an induced
subgraph (FNFL is reported separately there).  Both variants (with / without FNFL) are computed here.
Isomorphism is decided by graphs.canonical_form (exact, exhaustive over the refined partition).
Run from the repository root:  nice -n 10 python3 notes/referee_external/fixes_2026-11-01/doc_search/primitive_recount.py
"""
import json, sys, os, itertools, time
from fractions import Fraction
ROOT = '/home/user/CLTN'
sys.path.insert(0, ROOT)
from graphs import canonical_form, to_key
OUT = os.path.join(ROOT, 'notes/referee_external/fixes_2026-11-01/doc_search/primitive_recount.json')

def from_d6(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    for i in range(n): assert A[i][i] == 0
    return A

def induced(A, S):
    return [[A[i][j] for j in S] for i in S]

def inv(A):
    n = len(A)
    return (n, sum(map(sum, A)), tuple(sorted((sum(A[i]), sum(A[j][i] for j in range(n))) for i in range(n))))

# ---- n = 6: the three graphs
n6 = json.load(open(os.path.join(ROOT, 'certificates/n6_summary.json')))
hits6 = [h['A'] for h in n6['hits'] if not h['features'].get('is_complete', False)]
hits6 = [A for A in hits6 if sum(map(sum, A)) < 30]   # drop K6 (30 arcs) if present
assert len(hits6) == 3 and n6['limit']['limit_stable_noncomplete'] == 3, len(hits6)
C6 = {canonical_form(A): inv(A) for A in hits6}
assert len(C6) == 3

# ---- n = 7: 123 limit + FNFL
hits7 = []
with open(os.path.join(ROOT, 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')) as f:
    for line in f:
        r = json.loads(line); A = from_d6(r['d6'])
        assert len(A) == 7 and sum(map(sum, A)) == r['edges'] and not r['complete']
        hits7.append((r['d6'], A))
assert len(hits7) == 123
fnfl = json.load(open(os.path.join(ROOT, 'certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json')))
fnfl_d6 = fnfl.get('d6') or '&FNFL~[udo?'
assert fnfl_d6 == '&FNFL~[udo?', fnfl_d6
A_fnfl = from_d6(fnfl_d6)
C7 = {canonical_form(A): inv(A) for _, A in hits7}
assert len(C7) == 123
C7F = dict(C7); C7F[canonical_form(A_fnfl)] = inv(A_fnfl)
assert len(C7F) == 124

def contains(A, targets, k):
    """does A have an induced k-subgraph isomorphic to a member of targets (dict canon -> invariant)?"""
    invs = set(targets.values())
    for S in itertools.combinations(range(len(A)), k):
        B = induced(A, S)
        if inv(B) in invs and canonical_form(B) in targets:
            return True
    return False

# n = 7 primitive: no six-node counterexample induced
prim7 = [d6 for d6, A in hits7 if not contains(A, C6, 6)]
cont7 = 123 - len(prim7)

# ---- n = 8: 9,133 graphs from the aggregation certificate
cls = json.load(open(os.path.join(ROOT, 'certificates/n8_counterexamples_classes.json')))
hits8 = []
for row in cls['rows']:
    for d6 in row['members_d6']:
        A = from_d6(d6); assert len(A) == 8 and sum(map(sum, A)) == row['arcs']
        hits8.append((d6, A))
assert len(hits8) == 9133 and len(set(d for d, _ in hits8)) == 9133
t0 = time.time()
c6 = c7 = c7f = 0; prim8 = 0; prim8_strict = 0
for i, (d6, A) in enumerate(hits8):
    has6 = contains(A, C6, 6)
    has7 = contains(A, C7, 7)
    has7f = has7 or contains(A, {canonical_form(A_fnfl): inv(A_fnfl)}, 7)
    c6 += has6; c7 += has7; c7f += has7f
    if not (has6 or has7): prim8 += 1
    if not (has6 or has7f): prim8_strict += 1
    if i % 1000 == 0: print(i, time.time() - t0, flush=True)

res = {
  'counts': {'n6_counterexamples': 3, 'n7_limit_counterexamples': 123, 'n7_all_finite_eps': 124, 'n8_limit_counterexamples': 9133},
  'n7': {'contain_an_n6_counterexample': cont7, 'primitive_no_n6_induced': len(prim7)},
  'n8': {'contain_an_n6_counterexample': c6, 'contain_an_n7_limit_counterexample': c7,
         'contain_an_n7_limit_or_FNFL': c7f,
         'primitive_no_n6_or_n7_limit_induced': prim8,
         'primitive_no_n6_or_n7_limit_or_FNFL_induced': prim8_strict},
  'shares_exact': {'n6': str(Fraction(3, 3)), 'n7': str(Fraction(len(prim7), 123)), 'n8': str(Fraction(prim8, 9133)),
                   'n8_excluding_FNFL_containers': str(Fraction(prim8_strict, 9133))},
  'shares_decimal_display_only': {'n6': float(Fraction(3, 3)), 'n7': float(Fraction(len(prim7), 123)), 'n8': float(Fraction(prim8, 9133))},
  'monotone': {'absolute_counts_increase': 3 < len(prim7) < prim8,
               'shares_decrease': Fraction(3, 3) > Fraction(len(prim7), 123) > Fraction(prim8, 9133)},
}
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
