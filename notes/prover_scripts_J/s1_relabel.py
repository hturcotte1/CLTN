exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
n=6
best=[]
for p in itertools.permutations(range(n)):
    B = graphs.permute(G2, p)
    diff = [(i,j) for i in range(n) for j in range(n) if i!=j and G1[i][j]!=B[i][j]]
    best.append((len(diff), p, diff))
best.sort()
for b in best[:8]:
    print(b)
print('G1', show(G1))
p=best[0][1]; B=graphs.permute(G2,p); print('G2p', show(B))
