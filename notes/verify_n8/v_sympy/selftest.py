"""Positive and negative controls for vsympy8.check. Usage: python3 -I selftest.py <out_dir>"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vsympy8 import *
out = {}
# decoder: n = 6 Geneson graph from the n = 7 verifier self-test ('&EZjtwWE'); stable non-clique at eps=1/100, delta=29/2500
G = decode_d6('&EZjtwWE')
out['geneson_rows'] = [''.join(map(str, r)) for r in G]
out['geneson_reencode'] = encode_d6(G) == '&EZjtwWE'
r = check(G, '1/100', '29/2500'); out['geneson_control'] = {k: r[k] for k in ('legal', 'positive', 'stable', 'routh_agrees', 'closed_form_ok', 'charpoly_crosscheck')}
out['geneson_complete'] = is_complete_bidirected(G)
# K8 at a legal point: stable and positive (the conjecture's "if" direction)
K8 = [[0 if i == j else 1 for j in range(8)] for i in range(8)]
r = check(K8, *point_from_beta('3', 1000)); out['K8_control'] = {k: r[k] for k in ('legal', 'positive', 'stable', 'q')}
# K8 minus one arc at large q (q = 10 > n-2 = 6: Geneson's theorem says not stable-positive)
K8m = [row[:] for row in K8]; K8m[0][1] = 0
r = check(K8m, *point_from_beta('12', 1000)); out['K8_minus_arc_q11_control'] = {k: r[k] for k in ('legal', 'positive', 'stable', 'q')}
# d6 encoding spec check: '&' + N(8) = 'G' then 11 chars; identity-free small example: single arc 0->1 => first bit row 0 col 1
A = [[0] * 8 for _ in range(8)]; A[0][1] = 1
out['single_arc_0_to_1_d6'] = encode_d6(A)
out['single_arc_roundtrip'] = decode_d6(encode_d6(A)) == A
# first record of the input must decode to a graph with the recorded arc count
first = json.loads(open('/home/user/CLTN/n8/n8_limit_HITS.jsonl').readline())
Af = decode_d6(first['d6']); out['first_record_arcs'] = (sum(map(sum, Af)), first['edges'])
# Sturm: x^2 - 2 has exactly one root in (1, 2) and none in (2, 3)
P, seq = sturm_seq([1, 0, -2]); out['sturm_x2_minus_2'] = [count_roots_open(P, seq, 1, 2), count_roots_open(P, seq, 2, 3), [str(v) for v in refine_root(P, seq, 1, 2, 20)]]
# Hurwitz: lam^2 + 1 (eigenvalues +-i) must be UNDECIDED; lam^2 - 3 lam + 2 unstable; (lam+1)^8 stable
from sympy import S
m = hurwitz_minors([S(1), S(0), S(1)]); out['hurwitz_lam2_plus_1_minors'] = [str(v) for v in m]
out['hurwitz_unstable_minors'] = [str(v) for v in hurwitz_minors([S(1), S(-3), S(2)])]
from sympy import Poly, symbols
l = symbols('l'); c = [S(v) for v in Poly((l + 1) ** 8, l).all_coeffs()]
out['hurwitz_lam_plus_1_pow8_all_positive'] = all(v > 0 for v in hurwitz_minors(c))
json.dump(out, open(f'{sys.argv[1]}/selftest.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
