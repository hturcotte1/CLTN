import json,sys
sys.path.insert(0,'/home/user/CLTN')
from sympy import Matrix, symbols, factor, Rational, eye, ones, zeros
def d6(s):
    assert s[0]=='&'; s=s[1:]
    n=ord(s[0])-63; bits=[]
    for ch in s[1:]:
        v=ord(ch)-63
        bits += [(v>>(5-i))&1 for i in range(6)]
    return [[bits[i*n+j] for j in range(n)] for i in range(n)]
def rows(A): return ' '.join(''.join(map(str,r)) for r in A)
def edges(A): return sum(map(sum,A))
z=symbols('z')
def chi_perp(A):
    n=len(A); M=Matrix(A)
    V=zeros(n,n-1)
    for i in range(n-1): V[i,i]=1; V[n-1,i]=-1
    G=(V.T*V).inv()*V.T*M*V
    return factor(G.charpoly(z).as_expr())
for s,claimed in [('&FNFL~[udo?','0011110 0011100 1101111 1110111 0011011 0100101 1100000'),
                  ('&F^U|m^QTq?','0111110 1011011 1101101 1100111 1101001 0010101 1100100')]:
    A=d6(s); print(s, rows(A)==claimed, rows(A), edges(A))
A=d6('&FNFL~[udo?'); print('FNFL 7chi_perp', factor(7*chi_perp(A)))
# Geneson paper labelling
G1p=[[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
print('G1 paper rows', rows(G1p), edges(G1p))
print('geneson.d6', rows(d6('&EZjtwWE')))
for f in ['COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json','COUNTEREXAMPLE_CANDIDATE_n6_idx1424629.json']:
    d=json.load(open('/home/user/CLTN/certificates/'+f))
    print(f, d['idx'], rows(d['A_ours']), edges(d['A_ours']), d.get('isomorphic_to_geneson'))
    print('  features', d['features'])
    print('  stable_q_limit', d['stable_q_limit'])
