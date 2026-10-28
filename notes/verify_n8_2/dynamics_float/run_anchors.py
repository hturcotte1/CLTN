import sys, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/dynamics_float')
from dyn import *
import run_cands
from run_cands import one

def encode_d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(n + 63) + ''.join(chr(63 + int(''.join(map(str, bits[k:k + 6])), 2)) for k in range(0, len(bits), 6))

gen = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
assert encode_d6(gen) == '&EZjtwWE', encode_d6(gen)
K8 = [[int(i != j) for j in range(8)] for i in range(8)]
C3 = [[int(j == (i + 1) % 3) for j in range(3)] for i in range(3)]
anchors = [
    dict(name='geneson_6node', idx=1, d6='&EZjtwWE', eps='1/100', delta='29/2500', t='99', expect='converge'),
    dict(name='K8_complete', idx=2, d6=encode_d6(K8), eps='1/4', delta='1/2', t='3', expect='converge'),
    dict(name='3cycle', idx=3, d6=encode_d6(C3), eps='1/4', delta='1/2', t='3', expect='no_converge'),
    dict(name='3cycle_small_eps', idx=4, d6=encode_d6(C3), eps='1/100', delta='29/2500', t='99', expect='no_converge'),
]
for a in anchors:
    n, A = parse_d6(a['d6']); assert parse_d6(encode_d6(A))[1] == A
res = []
for a in anchors:
    t0 = time.time(); o = one(a); o['name'] = a['name']; o['expect'] = a['expect']
    print(a['name'], a['d6'], 'legal', o['legal'], 'q', o['q_check'], 'x_pos', o['exact']['all_pos'], 'max_re %.3g' % o['float']['max_re'],
          'T_long', o['T_long'], 'horizon_ok', o['horizon_sufficient'], 'RED', o['red_reasons'],
          'short_conv', o['n_short_converged'], '/', len(o['short']), 'long_conv', o['n_long_converged'], '/', len(o['long']),
          'long_ends_full_fp', o['n_long_ends_full_fp'], 'long errs', ['%.2g' % r['err_end'] for r in o['long']],
          'end supports', [r['end']['size'] for r in o['long']], 'wall', round(time.time() - t0, 1), flush=True)
    res.append(o)
json.dump(res, open('anchors.json', 'w'))
