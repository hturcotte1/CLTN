# Task 3 on G1 vs G2p (same labelling): every support, several legal (eps, delta). Exact (ctln.analyze_support).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
from itertools import combinations
n=6
G2p=[list(r) for r in graphs.permute(G2,(4,3,5,2,1,0))]
params=[(F(1,2),F(3)),(F(1,3),F(1)),(F(1,4),F(1,2)),(F(1,10),F(1,8)),(F(1,101),F(27,2626)),(F(1,5),F(1)),(F(1,20),F(1,10))]
allsame=True
for eps,delta in params:
    assert ctln.legal(eps,delta), (eps,delta)
    fp1=[];fp2=[];st1=[];st2=[]; diffs=0; bnd=0
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r1=ctln.analyze_support(G1,list(s),eps,delta); r2=ctln.analyze_support(G2p,list(s),eps,delta)
            if r1.get('BOUNDARY') or r2.get('BOUNDARY') or r1['stable'] is None or r2['stable'] is None: bnd+=1
            if r1['is_fixed_point']: fp1.append(s)
            if r2['is_fixed_point']: fp2.append(s)
            if r1['stable_fixed_point']: st1.append(s)
            if r2['stable_fixed_point']: st2.append(s)
            same=(r1['is_fixed_point']==r2['is_fixed_point'] and r1['stable']==r2['stable'] and r1['stable_fixed_point']==r2['stable_fixed_point']
                  and (not r1['is_fixed_point'] or (r1['x']==r2['x'] and r1['charpoly']==r2['charpoly'])))
            if not same: diffs+=1
    allsame &= (diffs==0)
    sc1=ctln.scan_supports(G1,eps,delta); sc2=ctln.scan_supports(G2p,eps,delta)
    print(f'eps={eps} delta={delta} q={ctln.q_of(eps,delta)}: FP(G1)={fp1}\n   FP(G2p)={fp2}\n   stable G1={st1} G2p={st2}  differing supports={diffs} undecided={bnd}'
          f'  scan_supports stable: {[r["sigma"] for r in sc1["stable_fixed_points"]]} vs {[r["sigma"] for r in sc2["stable_fixed_points"]]}')
print('ALL SUPPORTS IDENTICAL (FP, stability, x, charpoly) AT ALL PARAMS:', allsame)
