# Over all D-classes of the strict survivors: is the multiset of sizes of target-free cliques a D-invariant?
# If not, verify one witness pair exactly with ctln.analyze_support at a legal parameter.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle, collections
from itertools import combinations
n=6
data=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb'))['data']
def tfc(A):
    return tuple(sorted(len(s) for m in range(1,n+1) for s in combinations(range(n),m) if ctln.is_target_free_clique(A,list(s))))
dcls=collections.defaultdict(list)
for g in data: dcls[data[g]].append(g)
bad=[]
for c in dcls.values():
    if len(c)<2: continue
    prof={g:tfc(graph_of_idx(g)) for g in c}
    if len(set(prof.values()))>1: bad.append((c,prof))
print('non-singleton strict D-classes:',sum(1 for c in dcls.values() if len(c)>1),' with non-constant target-free-clique profile:',len(bad))
c,prof=bad[0]
gs=sorted(c,key=lambda g:prof[g]); g1,g2=gs[0],gs[-1]
print('witness class',c,'profiles',{g:prof[g] for g in c})
for eps,delta in ((F(1,4),F(1,2)),(F(1,10),F(1,8))):
    for g in (g1,g2):
        A=graph_of_idx(g)
        st=[]
        for m in range(1,n+1):
            for s in combinations(range(n),m):
                r=ctln.analyze_support(A,list(s),eps,delta)
                assert r['stable'] is not None and not r.get('BOUNDARY')
                if r['stable_fixed_point']: st.append(s)
        print(f'  eps={eps} delta={delta} q={ctln.q_of(eps,delta)} idx {g}: {show(A)}  stable supports {st}  (target-free cliques {prof[g]})')
