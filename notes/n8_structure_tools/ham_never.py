"""Hamiltonian directed cycle: tested exactly on all 347,638 non-complete strict survivors of n = 8 (hits and never).
Output: ham_never.json"""
import json, os
from multiprocessing import Pool
from collections import Counter
import lib
from more_universals import ham

HERE = os.path.dirname(os.path.abspath(__file__))

def f(line):
    r = json.loads(line)
    if r.get('complete'):
        return None
    k = 'hit' if r['verdict'] == 'stable_somewhere' else ('never_positivity' if 'a_i' in r.get('reason', '') else 'never_full')
    return k, ham(lib.from_d6(r['d6']))

if __name__ == '__main__':
    lines = open(lib.REPO + '/n8/n8_limit_results.jsonl').read().splitlines()
    with Pool(2) as p:
        c = Counter(x for x in p.imap(f, lines, chunksize=500) if x)
    out = {f'{k[0]} ham={k[1]}': v for k, v in sorted(c.items())}
    json.dump(out, open(os.path.join(HERE, 'ham_never.json'), 'w'), indent=1)
    print(out)
