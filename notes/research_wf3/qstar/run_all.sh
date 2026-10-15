#!/bin/sh
# Re-run every qstar computation and check (exact; single process). Run from anywhere.
D=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/qstar
cd /home/user/CLTN || exit 1
set -e
python3 $D/qstar.py --n6 --sharpen --out $D/qstar_n6.json > /dev/null
python3 $D/qstar.py --fnfl --sharpen --out $D/qstar_fnfl.json > /dev/null
python3 $D/qstar.py --n7-limit --out $D/qstar_n7limit.json > /dev/null
python3 $D/verify_fnfl_sup.py $D/qstar_fnfl.json > $D/verify_fnfl_sup.out
tail -1 $D/verify_fnfl_sup.out
python3 -c "
import json
d=json.load(open('$D/qstar_n6.json'))['n6_full_support']['aggregate']
assert d['exact'] and d['q_star']['minpoly_hi_to_lo']==[797,-1686,2159,-3004,2499,-2158,1521], d
assert d['upper_attained_by']==['idx1000255','idx1424629']
print('n6 q* = sextic root', d['q_star']['isolating_interval'], 'approx', d['q_star']['approx_decimal_NOT_A_CERTIFICATE'])
f=json.load(open('$D/qstar_fnfl.json'))['n7_finite_FNFL_CANDIDATE']['cells'][0]['sharpened']
assert f['status']=='exact' and f['sup']['degree']==25
print('FNFL sup (CANDIDATE) deg 25 in', f['sup']['isolating_interval'], 'approx', f['sup']['approx_decimal_NOT_A_CERTIFICATE'])
l=json.load(open('$D/qstar_n7limit.json'))['n7_limit_CANDIDATES_lower_bound']
assert l['max_beta_c_minus_1']['minpoly_hi_to_lo']==[1,0,-1,-1]
print('n7 limit lower bound (CANDIDATES): plastic number, attained by', len(l['attained_by']), 'graphs')
print('ALL QSTAR CHECKS PASSED')
"
