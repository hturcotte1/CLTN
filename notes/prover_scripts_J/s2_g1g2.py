exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
p=(4, 3, 5, 2, 1, 0)
G2p=[list(r) for r in graphs.permute(G2,p)]
print('G1 ', show(G1)); print('G2p', show(G2p))
# in-neighbourhoods of 1 and 2 in G1
inn=lambda A,v:{k for k in range(6) if A[k][v]}
out=lambda A,v:{k for k in range(6) if A[v][k]}
for A,name in ((G1,'G1'),(G2p,'G2p')):
    print(name, 'in(1)',inn(A,1),'in(2)',inn(A,2),'out(1)',out(A,1),'out(2)',out(A,2), 'in(4)', inn(A,4))
b1=fastregion.beta_data(G1,list(range(6))); b2=fastregion.beta_data(G2p,list(range(6)))
print('a vectors equal (same labelling):', b1['a']==b2['a'], ' d,phi,omega equal:', b1['d']==b2['d'], b1['phi']==b2['phi'], b1['omega']==b2['omega'])
N1=sp.Matrix(Nmat(G1)); N2=sp.Matrix(Nmat(G2p))
print('N2-N1 =', (N2-N1).tolist())
x=sp.Matrix([0,1,-1,0,0,0])
print('x^T N1 =', (x.T*N1).tolist(), ' (= -1 * x^T ?)', (x.T*N1 + x.T).is_zero_matrix)
print('A[1][2],A[2][1]=',G1[1][2],G1[2][1])
