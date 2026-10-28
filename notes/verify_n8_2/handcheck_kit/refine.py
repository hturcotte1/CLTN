"""Exact refinement of stored isolating intervals so the witness beta is separated from the endpoints."""
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/handcheck_kit')

def peval(coef_hi_to_lo, x):
    v = Fr(0)
    for c in coef_hi_to_lo:
        v = v * x + c
    return v

def refine(coef, lo, hi, beta, max_steps=200):
    """Bisect [lo,hi] (an isolating interval, sign change assumed) until beta is outside; return (lo,hi,steps,side)
    side = 'below' if beta < root interval, 'above' if beta > root interval."""
    slo, shi = peval(coef, lo), peval(coef, hi)
    assert slo != 0 and shi != 0 and (slo < 0) != (shi < 0), ('no sign change', lo, hi, slo, shi)
    steps = 0
    while not (beta < lo or beta > hi):
        if beta == lo or beta == hi:
            raise ValueError('beta coincides with endpoint')
        mid = (lo + hi) / 2
        sm = peval(coef, mid)
        if sm == 0:
            raise ValueError('beta root hit exactly')
        if (sm < 0) == (slo < 0):
            lo, slo = mid, sm
        else:
            hi, shi = mid, sm
        steps += 1
        if steps > max_steps:
            raise ValueError('too many steps')
    return lo, hi, steps, ('below' if beta < lo else 'above')

if __name__ == '__main__':
    hits = [json.loads(l) for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
    wit = {g['d6']: g for g in json.load(open('/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'))['graphs']}
    res = json.load(open('results.json'))
    out = {}
    for tag, r in res.items():
        h = hits[r['idx']]
        beta = Fr(r['beta'])
        iv = h['cert']['stable_intervals'][0]
        rep = {}
        for name in ('lo', 'hi'):
            e = iv[name]
            if 'beta' in e:
                b = Fr(e['beta'])
                rep[name] = dict(exact=str(b), ok=(beta > b if name == 'lo' else beta < b))
            else:
                lo, hi = Fr(e['beta_interval'][0]), Fr(e['beta_interval'][1])
                try:
                    lo2, hi2, steps, side = refine(e['beta_minpoly_hi_to_lo'], lo, hi, beta)
                    rep[name] = dict(refined=[str(lo2), str(hi2)], steps=steps, side=side,
                                     ok=(side == 'above' if name == 'lo' else side == 'below'))
                except (AssertionError, ValueError) as ex:
                    rep[name] = dict(error=str(ex))
        out[tag] = rep
        print(tag, r['d6'], 'beta', beta, json.dumps(rep))
    json.dump(out, open('interval_refine.json', 'w'), indent=1)
