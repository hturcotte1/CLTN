"""Small exact consistency checks (output consistency_extra.json):
 * the finite-eps-only n = 8 candidates (Theorem E member E_1; the 15 'cube' candidates of notes_theory 4.4)
   must NOT be in the limit list (their alpha = -1/2 exactly);
 * the two extremal graphs: relation (bidirected in-twin switch? D-equivalent?)."""
import json, os
import lib
HERE = os.path.dirname(os.path.abspath(__file__))
TH = lib.REPO + '/notes/research_wf3/theory/'
hits = [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
hc = set(lib.labelg_canon(hits))
E = [r['d6'] for r in json.load(open(TH + 'family_E_witnesses.json')) if len(lib.from_d6(r['d6'])) == 8]
C = [r['d6'] for r in json.load(open(TH + 'cube_n8_witnesses.json'))]
out = dict(E_n8=len(E), E_n8_in_limit_list=sum(c in hc for c in lib.labelg_canon(E)),
           cube_n8=len(C), cube_in_limit_list=sum(c in hc for c in lib.labelg_canon(C)) if C else 0)
a, b = '&GMx|\\Aw^qz]w', '&GVrm]A{^pZ]w'
A, B = lib.from_d6(a), lib.from_d6(b)
out['extremal_same_D_class'] = lib.D_class_key(A) == lib.D_class_key(B)
out['extremal_isomorphic'] = len(set(lib.labelg_canon([a, b]))) == 1
# search an in-twin switch: twins u,v (bidirected, same in-sets); switch the out-arc to some w between u and v
found = []
for (u, v) in lib.twin_pairs(A)['bi_in']:
    for w in range(8):
        if w in (u, v) or A[u][w] == A[v][w]:
            continue
        S = [row[:] for row in A]; S[u][w], S[v][w] = A[v][w], A[u][w]
        if lib.labelg_canon([lib.to_d6(S)])[0] == lib.labelg_canon([b])[0]:
            found.append(dict(twins=[u, v], w=w))
out['extremal_in_twin_switches_A_to_B'] = found
out['extremal_A_bidirected_in_twins'] = lib.twin_pairs(A)['bi_in']
json.dump(out, open(os.path.join(HERE, 'consistency_extra.json'), 'w'), indent=1)
print(out)
