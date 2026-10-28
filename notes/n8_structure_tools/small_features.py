"""Same exact features (features.gfeat) for the n = 6 and n = 7 counterexamples. Output: small_features.json"""
import json, os
import lib, features, small_lists as sl
HERE = os.path.dirname(os.path.abspath(__file__))
out = {'n6': [], 'n7': [], 'n7_FNFL': None}
for idx, d6 in sl.n6_list():
    r = features.gfeat(lib.from_d6(d6)); r['idx'] = idx; out['n6'].append(r)
for d6 in sl.n7_list():
    out['n7'].append(features.gfeat(lib.from_d6(d6)))
out['n7_FNFL'] = features.gfeat(lib.from_d6(sl.FNFL))
json.dump(out, open(os.path.join(HERE, 'small_features.json'), 'w'))
print('n6', len(out['n6']), 'n7', len(out['n7']))
