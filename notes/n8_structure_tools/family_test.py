"""Test the n = 6, 7, 8 universal statements on the PROVED families and certified instances at n >= 9:
 omega_bi <= n - 3 (bidirected clique bound), strong connectivity, induced oriented 3-cycle, and
 'receives from all => sends to all'.
 (a) explicit certified instances: notes/research_wf3/theory (F witnesses, high-q n = 8..11, Q/Q2/Y instances
     n = 19..73), Theorem U explicit n = 28, 41 (notes/theorem_U/...), R_k k = 0..40 and self-joins k = 1..4
     (Theorem A, notes/theorem_A_fixed/...);
 (b) symbolically for ALL k: every clique-class family G[S, all; k] (catalogue n6 + n7, Theorems F, Q, A) and the
     two-class families Q2, Y: omega = max over subsets J of classes of (sum_{j in J} k_j + omega(G[cap_J S_j])).
Output: family_test.json"""
import json, os, glob
from itertools import combinations
import lib
from families import class_ext

HERE = os.path.dirname(os.path.abspath(__file__))
TH = lib.REPO + '/notes/research_wf3/theory/'
UP = glob.glob(lib.REPO + '/notes/theorem_U/tmp/*/*/*/scratchpad/wf8/prover_q/')[0]
AP = glob.glob(lib.REPO + '/notes/theorem_A_fixed/tmp/*/*/*/scratchpad/wf7/prover_fixed/')[0]


def max_clique(A, nodes=None):
    n = len(A)
    nodes = list(range(n)) if nodes is None else list(nodes)
    N = {v: {w for w in nodes if w != v and A[v][w] and A[w][v]} for v in nodes}
    best = [0]
    def bk(R, P, X):
        if not P and not X:
            best[0] = max(best[0], len(R)); return
        if len(R) + len(P) <= best[0]:
            return
        u = max(P | X, key=lambda x: len(N[x] & P))
        for v in list(P - N[u]):
            bk(R | {v}, P & N[v], X & N[v]); P = P - {v}; X = X | {v}
    bk(set(), set(nodes), set())
    return best[0]


def has_C3(A):
    n = len(A)
    one = [[1 if (A[i][j] and not A[j][i]) else 0 for j in range(n)] for i in range(n)]
    outs = [[j for j in range(n) if one[i][j]] for i in range(n)]
    return any(one[k][i] for i in range(n) for j in outs[i] for k in outs[j])


def checks(A):
    n = len(A)
    recv_all = [v for v in range(n) if all(A[i][v] for i in range(n) if i != v)]
    return dict(n=n, omega=max_clique(A), omega_le_n_minus_3=max_clique(A) <= n - 3,
                strongly_connected=lib.strongly_connected(A), induced_C3=has_C3(A),
                recv_all_implies_send_all=all(all(A[v][j] for j in range(n) if j != v) for v in recv_all))


def main():
    inst = []
    cf = json.load(open(TH + 'family_F_certificate.json'))
    inst += [('F witness', w['d6']) for w in cf['F6_witnesses']]
    for f in ['top_n8.json', 'cert_q_910.json', 'cert_q_11.json']:
        inst += [('high-q ' + f, r['d6']) for r in json.load(open(TH + f)) if 'eps' in r]
    for f, key in [('family_Q_certificate.json', 'Q5_witness'), ('family_Q2_certificate.json', 'W5_witnesses'),
                   ('family_Y_certificate.json', 'Y5_witness')]:
        ws = json.load(open(TH + f))[key]
        ws = [ws] if isinstance(ws, dict) else ws
        inst += [(f, w['d6']) for w in ws]
    out = {'instances': []}
    for lab, d in inst:
        r = checks(lib.from_d6(d)); r.update(source=lab, d6=d); out['instances'].append(r)
    for r in json.load(open(UP + 'certificates_explicit.json')):
        rows = r['adjacency_rows']
        A = [[int(c) for c in row] for row in rows]
        c = checks(A); c.update(source='Theorem U explicit', certified=r['certified']); out['instances'].append(c)
    ta = json.load(open(AP + 'theorem_A_certificate.json'))
    R = [[int(c) for c in row] for row in ta['base_rows']]
    for k in range(0, 41):
        c = checks(class_ext([row[:] for row in R], ta['S'], k) if k else R); c.update(source='Theorem A R_k', k=k)
        out['instances'].append(c)
    for k in range(1, 5):
        n = 7 * k
        A = [[0] * n for _ in range(n)]
        for a in range(k):
            for b in range(k):
                for i in range(7):
                    for j in range(7):
                        u, v = 7 * a + i, 7 * b + j
                        if u != v:
                            A[u][v] = R[i][j] if a == b else 1
        c = checks(A); c.update(source='Theorem A self-join', k=k); out['instances'].append(c)
    # (b) symbolic families
    fams = []
    for fn in ('families_catalogue_n6.json', 'families_catalogue_n7.json'):
        for fam in json.load(open(TH + fn)):
            fams.append(('catalogue ' + fn, fam['base'].split(':')[-1], [(fam['S'], 1)], 0 if fam['proved_for_all_k>=0'] else 1))
    fams.append(('Theorem F', '&EWBFrtm', [([0, 1, 1, 1, 1, 1], 1)], 0))   # G1 in the catalogue labelling, S from catalogue
    fams.append(('Theorem Q', '&F^S{^ReL|?', [([1, 1, 1, 1, 0, 1, 1], 1)], 25))
    fams.append(('Theorem A', None, [(ta['S'], 1)], 0))
    fams.append(('Theorem Q2', '&F^U|m^QTq?', [([1, 1, 0, 1, 1, 1, 1], 1), ([1, 1, 1, 0, 1, 1, 1], 3)], 3))
    fams.append(('Theorem Y', '&G^ZZZr|KjqTo', [([0, 1, 1, 1, 1, 1, 1, 1], 4), ([1, 1, 0, 1, 1, 1, 1, 1], 1)], 13))
    sym = []
    for lab, d6, classes, k0 in fams:
        G = R if d6 is None else lib.from_d6(d6)
        m = len(G)
        # omega(k) = max_J (k * sum_{j in J} r_j + omega(G[cap S_j]));  n(k) = m + k * sum r_j
        # bound holds for all k >= k0 iff for every J:  k*sum_J r_j + w_J <= m + k*sum_all r_j - 3, all k >= max(k0,1)
        tot = sum(r for _, r in classes)
        ok = True; worst = None
        for size in range(len(classes) + 1):
            for J in combinations(range(len(classes)), size):
                nodes = [i for i in range(m) if all(classes[j][0][i] for j in J)]
                wJ = max_clique(G, nodes) if nodes else 0
                sJ = sum(classes[j][1] for j in J)
                kk = max(k0, 1)
                # linear in k with slope sJ - tot <= 0: the binding k is the smallest one
                slack = (m + kk * tot - 3) - (kk * sJ + wJ)
                if slack < 0: ok = False
                if worst is None or slack < worst: worst = slack
        sym.append(dict(family=lab, base=d6 or 'R (Theorem A)', m=m, k_from=max(k0, 1), bound_holds_for_all_k=ok, min_slack=worst))
    out['symbolic_families'] = dict(total=len(sym), holds=sum(1 for s in sym if s['bound_holds_for_all_k']),
                                    min_slack=min(s['min_slack'] for s in sym),
                                    tight=[s for s in sym if s['min_slack'] == 0][:20],
                                    n_tight=sum(1 for s in sym if s['min_slack'] == 0),
                                    failures=[s for s in sym if not s['bound_holds_for_all_k']])
    I = out['instances']
    out['instances_summary'] = {k: sum(1 for r in I if r[k]) for k in
                                ('omega_le_n_minus_3', 'strongly_connected', 'induced_C3', 'recv_all_implies_send_all')}
    out['instances_summary']['total'] = len(I)
    out['instances_summary']['n_range'] = [min(r['n'] for r in I), max(r['n'] for r in I)]
    out['instances_tight_omega'] = sum(1 for r in I if r['omega'] == r['n'] - 3)
    json.dump(out, open(os.path.join(HERE, 'family_test.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'instances'}, indent=1))
    print([(r['source'], r['n'], r['omega']) for r in I if not all(r[k] for k in ('omega_le_n_minus_3', 'strongly_connected', 'induced_C3', 'recv_all_implies_send_all'))])


if __name__ == '__main__':
    main()
