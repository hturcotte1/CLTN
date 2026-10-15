"""
validate_small.py -- run the generic exact primitives of exact_core.py against KNOWN small-n numbers (no floats).
 V1  n=6, Geneson's graph G1:  q_6* in (561/481, 7/6)  (REPORT_n7 item 8; limit-stable, positive at q = 561/481; not
     limit-stable at q = 7/6).
 V2  n=7, F_1 = G1[S,all;1]:   beta_c(1) = 2.2747... (notes_theory Theorem F table): bracket [22746/10000, 22748/10000].
 V3  n=8, F_2: positivity ends at the plastic number + 1 (b^3 - 3b^2 + 2b - 1 = 0, b = 2.32471795...): positive at
     23247/10000, not at 23248/10000; still limit-stable at 23248/10000 (beta_c(2) = 2.3391).
 V4  n=7 core B = &F^S{^ReL|? with S = (1,1,1,1,0,1,1): chi_L = z^3 (z^4 + z^3 + 3z^2 + z + 1) (Theorem Q), positive at
     beta = 3 and 10.
 V5  n=7 core B' = &F^U|m^QTq? (Theorem Q2): 4 chi_{D - bbar 1^T} = z^3 (2z^2 + z + 2)^2, bbar = (e_2 + 3 e_3)/4.
Run: python3 validate_small.py  (writes validate_small.json)
"""
import json, sys
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
from exact_core import *
from tlib import from_d6
import sympy as sp

G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
S1 = (1, 1, 1, 1, 1, 0)

def charpoly_rational(M):
    """det(zI - M) for a Fraction matrix via sympy (exact)."""
    z = sp.symbols('z')
    Mm = sp.Matrix([[sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else x for x in row] for row in M])
    return sp.Poly(sp.factor((z * sp.eye(len(M)) - Mm).det(method='berkowitz')), z)

out = {}
# V1
b_lo, b_hi = 1 + F(561, 481), 1 + F(7, 6)
v1 = {'stable_at_q=561/481': str(limit_stable(G1, b_lo)[0]), 'positive_at_q=561/481': positivity(G1, b_lo)['ok'],
      'stable_at_q=7/6': str(limit_stable(G1, b_hi)[0])}
v1['pass'] = v1['stable_at_q=561/481'] == 'True' and v1['positive_at_q=561/481'] and v1['stable_at_q=7/6'] == 'False'
out['V1_n6_G1_qstar_bracket'] = v1
# V2
F1 = sink_class_extension(G1, S1, 1)
v2 = {'n': len(F1), 'stable_at_2.2746': str(limit_stable(F1, F(22746, 10000))[0]),
      'stable_at_2.2748': str(limit_stable(F1, F(22748, 10000))[0]), 'positive_at_2.2746': positivity(F1, F(22746, 10000))['ok']}
v2['pass'] = v2['stable_at_2.2746'] == 'True' and v2['stable_at_2.2748'] == 'False' and v2['positive_at_2.2746']
out['V2_n7_F1_beta_c'] = v2
# V3
F2 = sink_class_extension(G1, S1, 2)
v3 = {'n': len(F2), 'positive_at_2.3247': positivity(F2, F(23247, 10000))['ok'], 'positive_at_2.3248': positivity(F2, F(23248, 10000))['ok'],
      'stable_at_2.3248': str(limit_stable(F2, F(23248, 10000))[0]),
      'plastic_poly_sign_change': (lambda f: (f(F(23247, 10000)) < 0 < f(F(23248, 10000))))(lambda b: b**3 - 3*b*b + 2*b - 1)}
v3['pass'] = v3['positive_at_2.3247'] and not v3['positive_at_2.3248'] and v3['stable_at_2.3248'] == 'True' and v3['plastic_poly_sign_change']
out['V3_n8_F2_plastic_positivity'] = v3
# V4
B = from_d6('&F^S{^ReL|?'); S = (1, 1, 1, 1, 0, 1, 1)
D = D_of_A(B); m = len(D); b = [1 - s for s in S]
L = [[F(D[i][j] - b[i]) for j in range(m)] for i in range(m)]
z = sp.symbols('z')
chiL = charpoly_rational(L)
v4 = {'chi_L': str(chiL.as_expr()), 'pass_chi_L': sp.expand(chiL.as_expr() - z**3 * (z**4 + z**3 + 3*z**2 + z + 1)) == 0,
      'positive_beta3': positivity(B, 3)['ok'], 'positive_beta10': positivity(B, 10)['ok']}
v4['pass'] = v4['pass_chi_L'] and v4['positive_beta3'] and v4['positive_beta10']
out['V4_n7_core_Q'] = v4
# V5
B2 = from_d6('&F^U|m^QTq?'); D2 = D_of_A(B2); m = len(D2)
bbar = [F(0)] * m; bbar[2] += F(1, 4); bbar[3] += F(3, 4)
L2 = [[F(D2[i][j]) - bbar[i] for j in range(m)] for i in range(m)]
chiL2 = charpoly_rational(L2)
v5 = {'chi_L': str(chiL2.as_expr()), 'pass': sp.expand(4 * chiL2.as_expr() - z**3 * (2*z**2 + z + 2)**2) == 0}
out['V5_n7_core_Q2'] = v5
out['ALL_PASS'] = all(v['pass'] for v in out.values() if isinstance(v, dict))
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/validate_small.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
