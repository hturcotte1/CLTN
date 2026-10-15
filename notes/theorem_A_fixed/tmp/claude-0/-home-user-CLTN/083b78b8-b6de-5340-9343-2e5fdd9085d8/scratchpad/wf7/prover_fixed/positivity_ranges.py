"""positivity_ranges.py -- EXACT: positivity of every family member (Lemma X2') on [2, B] for B covering each family's ceiling:
R, Q on [2, 3] (beta_inf = 2.844...), Q2 on [2, 4], Y on [2, 5] (beta_Y = 4.9607...).  Uses theorem_B.positivity_interval
(det, adj(Mb)1 entries and sigma_j - 1 numerators: no root in [2, B], sign of det).  Writes positivity_ranges.json."""
import json
import fixlib as fl
from theorem_B import FAMILIES, positivity_interval
out = {}
for nm, B in (('R', 3), ('Q', 3), ('Q2', 4), ('Y', 5)):
    base, classes, rhos = FAMILIES[nm]
    ok, rec = positivity_interval(fl.from_d6(base), classes, B)
    out[nm] = {'interval': [2, B], 'positive_on_interval': ok, 'detail': rec}
    print(nm, B, ok); assert ok
json.dump(out, open('positivity_ranges.json', 'w'), indent=1)
