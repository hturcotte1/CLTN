import json, sys
sys.path.insert(0,'/home/user/CLTN/notes/verify_n8_2/dynamics_float')
hits=[json.loads(l) for l in open('/home/user/CLTN/n8/n8_limit_HITS.jsonl')]
wit=json.load(open('/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'))['graphs']
assert [h['d6'] for h in hits]==[w['d6'] for w in wit]
sel={}
for i in range(0,len(hits),23):
    sel[i]=('every23',)
for i,w in enumerate(wit):
    if w['witness_stored_legal'] is None:
        sel[i]=sel.get(i,())+('null_stored_legal',)
cands=[]
for i in sorted(sel):
    w=wit[i]['witness_in_first_stable_cell']
    cands.append(dict(idx=i,d6=hits[i]['d6'],why=list(sel[i]),eps=w['eps'],delta=w['delta'],beta=w['beta'],t=w['t'],q=w['q'],complete=hits[i]['complete']))
json.dump(cands,open('candidates_sample.json','w'),indent=0)
print('candidates',len(cands),'every23',sum(1 for c in cands if 'every23' in c['why']),'null',sum(1 for c in cands if 'null_stored_legal' in c['why']))
# never sample: every 1692nd never record (0-based index among never records), k=0..199
nev=[]; k=0
for line_no,l in enumerate(open('/home/user/CLTN/n8/n8_limit_results.jsonl')):
    r=json.loads(l)
    if r['verdict']!='never': continue
    if k%1692==0 and len(nev)<200:
        nev.append(dict(never_idx=k,line=line_no,d6=r['d6'],reason=r.get('reason'),has_cert='cert' in r))
    k+=1
json.dump(nev,open('never_sample.json','w'),indent=0)
from collections import Counter
print('never',len(nev),Counter(n['reason'] for n in nev))
