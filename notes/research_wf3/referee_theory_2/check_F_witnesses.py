"""Referee: recompute every Theorem F witness with ctln.analyze_support at the stated (eps, delta); add own witnesses.
Exact Fractions only."""
import sys, json
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as Fr
import ctln
G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
def Fk(kk):
    n = 6 + kk; A = [[0]*n for _ in range(n)]
    for i in range(6):
        for j in range(6): A[i][j] = G1[i][j]
    for v in range(6, n):
        for i in range(6):
            A[v][i] = 1; A[i][v] = 1 if i <= 4 else 0
        for u in range(6, n):
            if u != v: A[v][u] = 1
    return A
def d6dec(s):
    assert s[0] == '&'; n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - i)) & 1 for i in range(6)]
    return [[bits[i*n + j] for j in range(n)] for i in range(n)]
def check(A, eps, delta, label):
    n = len(A)
    rec = ctln.analyze_support(A, list(range(n)), eps, delta)
    lg = ctln.legal(eps, delta)
    ok = lg and rec['positive'] and rec['stable'] is True and rec['stability_detail'] == 'all Hurwitz minors positive' and not rec['clique'] and rec['stable_fixed_point']
    # identity cross-check of I - W closed form
    assert ctln.identity_check_IminusW(A, list(range(n)), eps, delta)
    print(label, 'n=%d eps=%s delta=%s q=%s legal=%s pos=%s stable=%s (%s) clique=%s -> %s' % (n, eps, delta, ctln.q_of(eps, delta), lg, rec['positive'], rec['stable'], rec['stability_detail'], rec['clique'], 'OK' if ok else 'FAIL'))
    return ok, rec
cert = json.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/family_F_certificate.json'))
allok = True
for wtn in cert['F6_witnesses']:
    A = d6dec(wtn['d6']); kk = wtn['k']
    assert A == Fk(kk), ('d6 does not decode to F_k', kk)
    eps = Fr(wtn['eps']); delta = Fr(wtn['delta'])
    t = (1 - eps) / eps; beta = 1 + delta / eps
    assert t == Fr(wtn['t']) and beta == Fr(wtn['beta'])
    ok, _ = check(A, eps, delta, 'cert k=%d' % kk); allok &= ok
fw = json.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/F_family_witnesses.json'))
for wtn in fw:
    A = d6dec(wtn['d6']); kk = wtn['k']
    assert A == Fk(kk) and [''.join(map(str, r)) for r in A] == wtn['A_rows']
    ok, rec = check(A, Fr(wtn['eps']), Fr(wtn['delta']), 'F_family_witnesses k=%d' % kk); allok &= ok
    assert rec['x'] == wtn['x'], 'x mismatch'
own = []
def wit(kk, beta, ts):
    A = Fk(kk)
    for t in ts:
        t = Fr(t); eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        if not ctln.legal(eps, delta): continue
        ok, rec = check(A, eps, delta, 'own k=%d beta=%s t=%s' % (kk, beta, t))
        if ok:
            own.append({'k': kk, 'n': 6 + kk, 'beta': str(beta), 't': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)), 'detail': rec['stability_detail']})
            return True
    return False
res = {}
res['k0_beta_2.1'] = wit(0, Fr(21, 10), [20, 100, 1000])
res['k7_beta_2.1'] = wit(7, Fr(21, 10), [20, 100, 1000])
res['k9_beta_2.1'] = wit(9, Fr(21, 10), [20, 100, 1000])
res['k12_beta_2.1'] = wit(12, Fr(21, 10), [20, 100, 1000])
res['k20_beta_2.1'] = wit(20, Fr(21, 10), [20, 100, 1000])
# near the plastic boundary: beta = 23247/10000 < b_p, for k = 2..8 (claimed q-sup = plastic)
for kk in [2, 5, 8]:
    res['k%d_beta_2.3247' % kk] = wit(kk, Fr(23247, 10000), [1000, 10**4, 10**5, 10**6])
# just above b_p positivity must fail (any t)
for kk in [2, 5]:
    A = Fk(kk)
    for t in [100, 10**4]:
        t = Fr(t); beta = Fr(2325, 1000); eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        rec = ctln.analyze_support(A, list(range(len(A))), eps, delta)
        print('above b_p k=%d t=%s positive=%s' % (kk, t, rec['positive'])); assert rec['positive'] is False
print(res)
json.dump({'own_witnesses': own, 'results': res, 'all_cert_witnesses_ok': allok}, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/referee_theory_2/witness_check.json', 'w'), indent=1)
print('ALL CERT WITNESSES OK' if allok else 'SOME CERT WITNESS FAILED')
