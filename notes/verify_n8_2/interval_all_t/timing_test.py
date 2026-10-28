import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/interval_all_t')
import engine as E
recs=[json.loads(l) for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
cw=json.load(open('/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'))['graphs']
for idx in [0, 103]:
    r=recs[idx]; g=cw[idx]
    n,A=E.parse_digraph6(r['d6']); assert n==8
    assert sum(map(sum,A))==r['edges']
    beta=F(g['witness_in_first_stable_cell']['beta'])
    t0=time.time(); print('Wcheck', E.build_K_vs_W_check(A, beta.numerator, beta.denominator), time.time()-t0)
    t0=time.time(); conds,extra=E.symbolic_conditions(A,beta); t1=time.time()
    print('symbolic', t1-t0, {k:len(v)-1 for k,v in conds.items()})
    t0=time.time(); res=E.sign_set(conds); t2=time.time()
    print('signset', t2-t0)
    print('nroots', len(res['roots']), 'nfacts', len(res['facts']), [f.degree() for f in res['facts']])
    for (lo,hi,fi),v in zip(res['roots'],res['root_conds']): print(' root', float(lo), float(hi), 'fact',fi, v)
    print('S', [( (float(l[0]) if l else 0), (float(r_[1]) if r_ else 'inf')) for l,r_ in res['S']])
    # pointwise cross-check at t=1000
    for t in [100,1000,10000]:
        eps=F(1,t+1); delta=(beta-1)*eps
        pw=E.pointwise(A,eps,delta)
        # compare Hurwitz minors sign with symbolic
        sy=[E.peval(extra['H'][k],F(t)) for k in range(8)]
        print(' t',t,pw['stable'],pw['positive'],pw['detail'], 'signs agree', all((m>0)==(s>0) and (m<0)==(s<0) for m,s in zip(pw['minors'],sy)),
              'x agree', all((xi>0)==(E.peval(extra['N'][i],F(t))*E.peval(extra['D'],F(t))>0) for i,xi in enumerate(pw['x'])))
    print('beta in interval', E.beta_in_interval(beta, r['cert']['stable_intervals'][0]))
