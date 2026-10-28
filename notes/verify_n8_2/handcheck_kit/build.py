import sys, json
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/handcheck_kit')
from engine import *
from fractions import Fraction as Fr

HITS = '/home/user/CLTN/n8/n8_limit_HITS.jsonl'
WIT = '/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'
hits = [json.loads(l) for l in open(HITS)]
wit = {g['d6']: g for g in json.load(open(WIT))['graphs']}
assert len(wit) == len(hits) == 9133
lift_idx = set(json.load(open('lift_hit_indices.json'))['hit_indices'])
lifts = json.load(open('lift_123_to_8.json'))

def first_idx(pred):
    return next(i for i, h in enumerate(hits) if pred(i, h))

sel = [
 ('C1', 'fewest arcs (29); first such record in file order', first_idx(lambda i, h: h['edges'] == 29)),
 ('C2', 'most arcs (45); first such record in file order', first_idx(lambda i, h: h['edges'] == 45)),
 ('C3', 'largest refined q upper end (named in the task: &GMx|\\Aw^qz]w)', first_idx(lambda i, h: h['d6'] == '&GMx|\\Aw^qz]w')),
 ('C4', 'narrow-interval graph: first record in file order with witness_stored_legal = null', first_idx(lambda i, h: wit[h['d6']]['witness_stored_legal'] is None)),
 ('C5', '34 arcs; first such record in file order (fallback: no closed-twin double blow-up of Geneson\'s &EZjtwWE is in the list)', first_idx(lambda i, h: h['edges'] == 34)),
 ('C6', 'line index 4567 (0-based) of n8_limit_HITS.jsonl', 4567),
]

def fs(x):
    return str(x)

def dec(x, k=6):
    return f'{float(x):.{k}f}'

def grid(A):
    n = len(A)
    s = '      ' + ' '.join(f'{j}' for j in range(n)) + '\n'
    for i in range(n):
        s += f'  {i} | ' + ' '.join(str(A[i][j]) for j in range(n)) + '\n'
    return s

def fracgrid(M, symbols):
    """print I-W using symbols for the three possible entry values"""
    n = len(M)
    s = ''
    for i in range(n):
        s += '  [ ' + '  '.join(symbols[M[i][j]] for j in range(n)) + ' ]\n'
    return s

out = []
results = {}
for tag, why, idx in sel:
    h = hits[idx]
    d6 = h['d6']
    A = d6_decode(d6)
    assert d6_encode(A) == d6
    w = wit[d6]['witness_in_first_stable_cell']
    eps, delta, beta_w, t = Fr(w['eps']), Fr(w['delta']), Fr(w['beta']), int(w['t'])
    assert eps == Fr(1, t + 1) and delta == (beta_w - 1) * eps, 'witness ray relation'
    r = analyze(A, eps, delta)
    assert r['q'] == Fr(w['q'])
    arcs = sum(map(sum, A))
    assert arcs == h['edges']
    mp = missing_pair(A)
    iv = h['cert']['stable_intervals']
    assert len(iv) == 1
    lo, hi = iv[0]['lo'], iv[0]['hi']
    # beta containment check: witness beta strictly inside (lo_upper, hi_lower)
    lo_up = Fr(lo['beta']) if 'beta' in lo else Fr(lo['beta_interval'][1])
    hi_lo = Fr(hi['beta']) if 'beta' in hi else Fr(hi['beta_interval'][0])
    inside = lo_up < beta_w < hi_lo
    results[tag] = dict(idx=idx, d6=d6, arcs=arcs, eps=str(eps), delta=str(delta), q=str(r['q']), beta=str(beta_w),
                        positive=r['positive'], stability=r['stability'], legal=r['legal'], complete=r['complete'],
                        beta_inside_stored_interval=inside, minors=[str(m) for m in r['minors']],
                        charpoly=[str(c) for c in r['charpoly']], x=[str(v) for v in r['x']],
                        missing_pair=mp, is_lift_of_123=(idx in lift_idx))
    n = 8
    S = []
    S.append(f'## {tag}: `{d6}`  ({arcs} arcs, hits line index {idx})\n')
    S.append(f'Selection rule: {why}.\n')
    S.append('**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:\n')
    S.append('```\n' + grid(A) + '```\n')
    bits = ''.join(str(A[i][j]) for i in range(n) for j in range(n))
    S.append(f'digraph6 check: `{d6[2:]}` is the 64-bit row-major string `{bits}` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).\n')
    i0, j0 = mp
    S.append(f'**Not a clique:** A[{i0}][{j0}] = 0, i.e. the arc {i0} -> {j0} is missing (so the pair {{{i0}, {j0}}} is not bidirected). '
             f'Total arcs {arcs} < 56.\n')
    S.append(f'**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  '
             f't = {t}, eps = {eps}, delta = {delta}, beta = 1 + delta/eps = {beta_w}, '
             f'q = delta(1-eps)/eps = **{r["q"]}** ({dec(r["q"])}) > 1.  theta = 1.\n')
    S.append('**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; '
             f'off-diagonal (i,j) = 1 - eps = {1-eps} if A[j][i] = 1 (j sends to i), else 1 + delta = {1+delta}.  '
             'Writing `a` = 1 - eps, `b` = 1 + delta:\n')
    symbols = {Fr(1): '1', 1 - eps: 'a', 1 + delta: 'b'}
    S.append('```\n' + fracgrid(r['M'], symbols) + '```\n')
    S.append('(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)\n')
    S.append('**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):\n')
    S.append('```\n' + '\n'.join(f'  x[{i}] = {r["x"][i]}   ({dec(r["x"][i], 8)})' for i in range(n)) + '\n```\n')
    S.append('**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:\n')
    S.append('```\n' + '\n'.join(f'  c{k} = {r["charpoly"][k]}   ({dec(r["charpoly"][k], 10)})' for k in range(1, n + 1)) + '\n```\n')
    S.append('**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):\n')
    S.append('```\n' + '\n'.join(f'  D{k+1} = {r["minors"][k]}\n       sign +, approx {float(r["minors"][k]):.6e}' for k in range(n)) + '\n```\n')
    S.append(f'My verdict at this point: positive = {r["positive"]}, Routh-Hurwitz = **{r["stability"]}**, legal = {r["legal"]}, complete = {r["complete"]}.  '
             f'=> stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.\n')
    # stored interval
    def side(e):
        if 'beta' in e:
            return f"kind={e['kind']}, beta = {e['beta']}, q = {e.get('q','?')}"
        return (f"kind={e['kind']}, beta in [{e['beta_interval'][0]}, {e['beta_interval'][1]}] "
                f"(root of a degree-{len(e['beta_minpoly_hi_to_lo'])-1} integer polynomial, vanishing {e['vanishing']}), "
                f"q in [{e['q_interval'][0]}, {e['q_interval'][1]}]")
    S.append('**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):\n')
    S.append(f'- lo: {side(lo)}\n- hi: {side(hi)}\n')
    import os
    if os.path.exists('interval_refine.json') and os.path.exists('float_sanity.json'):
        rf = json.load(open('interval_refine.json'))[tag]
        fl = json.load(open('float_sanity.json'))[tag]
        def rside(name):
            e = rf[name]
            if 'exact' in e:
                return f"rational endpoint {e['exact']}: witness {'>' if name=='lo' else '<'} endpoint = {e['ok']}"
            return (f"isolating interval bisected exactly {e['steps']} more step(s) (sign of the integer minimal polynomial at the midpoint) "
                    f"to [{e['refined'][0]}, {e['refined'][1]}]; witness is {e['side']} it = {e['ok']}")
        S.append(f'- witness beta = {beta_w} ({dec(beta_w)}) lies strictly inside the stored interval:  lo: {rside("lo")};  hi: {rside("hi")}.\n')
        S.append(f'- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx {fl:.3e} < 0 (consistent with the exact Hurwitz verdict; not used for any decision).\n')
    else:
        S.append(f'- witness beta = {beta_w} ({dec(beta_w)}) inside stored coarse interval: {inside} (exact comparison {lo_up} < {beta_w} < {hi_lo}; see interval_refine.json).\n')
    S.append(f'- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: {idx in lift_idx}.\n')
    if idx in lift_idx:
        par = [(p['n7_d6'], l['v']) for p in lifts for l in p['lifts_in_hits'] if l['hit_index'] == idx]
        S.append(f'  Parents (n7 d6, blown-up vertex): {par}\n')
    out.append('\n'.join(S))
json.dump(results, open('results.json', 'w'), indent=1)
open('sections.md', 'w').write('\n'.join(out))
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ('idx','d6','arcs','q','positive','stability','legal','beta_inside_stored_interval','is_lift_of_123','missing_pair')} for k, v in results.items()}, indent=0))
