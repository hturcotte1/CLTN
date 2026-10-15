"""r5_N6.py -- referee: N = 6, W = 2 instance: C1 (limit polynomial strip test) and positivity only (k0 not computed)."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q')
from r3_instances import rounded, positivity_tree
from ref_lib import *
out = {}
for W in [1, 2]:
    c = rounded(6, W); p, G, m = reduced_from_c(c)
    out[W] = {'c': c, 'm': m, 'C1': strip_right_of(p, sp.Rational(1, 6)), 'pos': positivity_tree(c, 6)['ok']}
print(json.dumps(out, default=str))
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/referee_q/r5_N6.json', 'w'), default=str, indent=1)
