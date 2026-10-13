exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
from itertools import combinations
n=6
G2p=[list(r) for r in graphs.permute(G2,(4,3,5,2,1,0))]
params=[(F(1,2),F(3)),(F(1,3),F(1)),(F(1,4),F(1,2)),(F(1,10),F(1,8)),(F(1,101),F(27,2626)),(F(1,5),F(1)),(F(1,20),F(1,10))]
ok=True
for eps,delta in params:
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r1=ctln.analyze_support(G1,list(s),eps,delta); r2=ctln.analyze_support(G2p,list(s),eps,delta)
            keys_=('is_fixed_point','positive','off_support_ok','stable','stable_fixed_point','x','charpoly')
            d=[k for k in keys_ if r1.get(k)!=r2.get(k)]
            if d:
                one = len(set(s)&{1,2})==1
                nonfp = (not r1['is_fixed_point']) and (not r2['is_fixed_point'])
                print(eps,delta,s,'differs in',d,' contains exactly one of {1,2}:',one,' non-FP in both:',nonfp)
                ok &= one and nonfp
print('every difference is on a support with exactly one twin and is non-FP in both:', ok)
