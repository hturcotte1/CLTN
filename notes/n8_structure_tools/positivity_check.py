"""Consistency of the own exact positivity sets (positivity.jsonl) with the stored stable intervals, and the
mechanism classification of each stable interval (lower end: 2 or a positivity root; upper end: spectral
boundary beta_c with positivity continuing, or a positivity root).  Output: positivity_check.json"""
import json, os
from collections import Counter
import sympy as sp
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def val(s):
    return sp.oo if s == 'oo' else sp.sympify(s)


def main():
    P = [json.loads(l) for l in open(os.path.join(HERE, 'positivity.jsonl'))]
    H = [json.loads(l) for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
    W = json.load(open(lib.REPO + '/certificates/n8_candidates_legal_witnesses.json'))['graphs']
    c = Counter(); bad = []
    for p, h, w in zip(P, H, W):
        assert p['d6'] == h['d6'] == w['d6']
        beta = sp.Rational(w['witness_in_first_stable_cell']['beta'])
        comp = [k for k in p['P'] if val(k[0]) < beta < val(k[1])]
        if len(comp) != 1:
            bad.append(h['d6']); continue
        lo_, hi_ = val(comp[0][0]), val(comp[0][1])
        st = h['cert']['stable_intervals'][0]
        lo, hi = st['lo'], st['hi']
        # lower end
        if lo['kind'] == 'beta_min':
            ok_lo = lo_ == 2
            ltype = 'starts at 2'
        else:
            a, b_ = [sp.Rational(x) for x in lo['beta_interval']]
            ok_lo = a <= lo_ <= b_
            ltype = 'starts at a positivity root > 2'
        a, b_ = [sp.Rational(x) for x in hi['beta_interval']]
        if hi['vanishing'] == ['Delta']:
            ok_hi = hi_ > b_          # positivity continues beyond the spectral end
            utype = 'spectral end, positivity continues' + (' to oo' if hi_ == sp.oo else ' (bounded)')
        else:
            ok_hi = a <= hi_ <= b_
            utype = 'ends at a positivity root'
        ok_lo, ok_hi = bool(ok_lo), bool(ok_hi)
        c["ok"] += int(ok_lo and ok_hi)
        c[(ltype, utype)] += 1
        if not (ok_lo and ok_hi):
            bad.append(h['d6'])
    out = dict(consistent=c.pop('ok'), inconsistent=bad,
               mechanism={f'{k[0]} | {k[1]}': v for k, v in sorted(c.items())},
               positivity_components=dict(Counter(p['components'] for p in P)),
               positivity_unbounded=sum(1 for p in P if p['unbounded']))
    json.dump(out, open(os.path.join(HERE, 'positivity_check.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
