import sys, json, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/interval_all_t')
import engine as E
from run import analyse_beta, in_S
out={}
# K8: complete graph must be stable+positive for every t>0 (S_sp = (0,inf)), legal for t > t_L
A=[[1 if i!=j else 0 for j in range(8)] for i in range(8)]
for beta in [F(2001,1000), F(5,2), F(7,2)]:
    an,_,_=analyse_beta(A,beta,'K8')
    out['K8 beta=%s'%beta]={'stabpos_roots':an['stab_pos_only']['n_roots'],'stabpos_unbounded':an['stab_pos_only']['S_contains_all_large_t'],
                            'full_S':an['full']['S'],'t_legal':an['legal_threshold_t']}
# empty graph: must never be stable (S_sp empty)
A0=[[0]*8 for _ in range(8)]
an,_,_=analyse_beta(A0,F(5,2),'empty'); out['empty beta=5/2']={'S_sp':an['stab_pos_only']['S'],'unbounded':an['stab_pos_only']['S_contains_all_large_t']}
# negative controls: 'never' graphs from n8_limit_results.jsonl (reservoir sample, fixed seed)
random.seed(12345)
never_full=[]; never_pos=[]
with open('/home/user/CLTN/n8/n8_limit_results.jsonl') as f:
    for k,l in enumerate(f):
        r=json.loads(l)
        if r.get('verdict')!='never': continue
        tgt = never_full if 'cert' in r else never_pos
        if len(tgt)<40: tgt.append(r)
        else:
            j=random.randrange(k+1)
            if j<40: tgt[j]=r
        if k>2000000: break
res={'never_full':[], 'never_positivity':[]}
for name,lst in [('never_full',never_full),('never_positivity',never_pos)]:
    for r in lst:
        n,A=E.parse_digraph6(r['d6'])
        row={'d6':r['d6']}
        for beta in [F(2001,1000),F(21,10),F(5,2),F(4)]:
            an,_,_=analyse_beta(A,beta,'ctrl')
            row[str(beta)]={'unbounded':an['stab_pos_only']['S_contains_all_large_t'],'S_sp_components':len(an['stab_pos_only']['S'])}
        res[name].append(row)
out['never_controls']={k:{'n':len(v),'any_unbounded':sum(1 for row in v for b,x in row.items() if b!='d6' and x['unbounded']),
                          'any_S_nonempty':sum(1 for row in v for b,x in row.items() if b!='d6' and x['S_sp_components']>0)} for k,v in res.items()}
# n=7 limit counterexamples: witnesses beta, check unbounded at a stored witness beta that is stable
n7=[json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
random.seed(7); samp=random.sample(n7,12)
n7res=[]
for r in samp:
    n,A=E.parse_digraph6(r['d6']); assert n==7
    cells=[c for c in r['cert']['cells'] if c['stable'] and c['positive']]
    bad=[c for c in r['cert']['cells'] if not (c['stable'] and c['positive'])]
    beta=F(cells[0]['sample_beta']); an,_,_=analyse_beta(A,beta,'n7')
    row={'d6':r['d6'],'beta_stable_cell':str(beta),'unbounded_at_stable_cell_beta':an['stab_pos_only']['S_contains_all_large_t'],'t_star':an['full']['t_star']}
    if bad:
        b2=F(bad[-1]['sample_beta']); an2,_,_=analyse_beta(A,b2,'n7b')
        row.update({'beta_unstable_cell':str(b2),'unbounded_at_unstable_cell_beta':an2['stab_pos_only']['S_contains_all_large_t']})
    n7res.append(row)
out['n7_sample']=n7res
json.dump(out,open('controls.json','w'),indent=1,default=str)
print(json.dumps(out,indent=1,default=str)[:6000])
