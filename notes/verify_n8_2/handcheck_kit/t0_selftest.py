import sys, json
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/handcheck_kit')
from engine import *
from fractions import Fraction as Fr
# Geneson stage-0 check
G = json.load(open('/home/user/CLTN/certificates/stage0_geneson_2607.21396.json'))
A = G['A_ours']
print('Geneson d6 (my encoder):', d6_encode(A), ' decode roundtrip ok:', d6_decode(d6_encode(A)) == A)
r = analyze(A, Fr(1, 10), Fr(29, 250))
print('Geneson eps=1/10 delta=29/250: q=', r['q'], 'positive', r['positive'], r['stability'])
print(' x =', [str(v) for v in r['x']])
print(' stored x[0] =', G['records_by_eps']['1/10']['x'][:2])
# complete graph K4 bidirected must be stable at legal point; 3-cycle must be unstable/not positive
K4 = [[int(i != j) for j in range(4)] for i in range(4)]
print('K4:', analyze(K4, Fr(1, 4), Fr(1, 2))['stability'], analyze(K4, Fr(1, 4), Fr(1, 2))['positive'])
C3 = [[0,1,0],[0,0,1],[1,0,0]]
r = analyze(C3, Fr(1, 4), Fr(1, 2)); print('C3:', r['stability'], r['positive'])
# Hurwitz sanity on known polynomial s^2 + s + 1 (stable) and s^2 - s + 1 (unstable), s^3+s^2+s+1 (boundary: roots +-i)
print(classify(hurwitz_minors([Fr(1),Fr(1),Fr(1)])[1]), classify(hurwitz_minors([Fr(1),Fr(-1),Fr(1)])[1]), classify(hurwitz_minors([Fr(1),Fr(1),Fr(1),Fr(1)])[1]))
