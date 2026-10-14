import json
def d6(A):
    bits=[A[i][j] for i in range(len(A)) for j in range(len(A))]; bits+=[0]*((-len(bits))%6)
    return '&'+chr(63+len(A))+''.join(chr(63+int(''.join(map(str,bits[k:k+6])),2)) for k in range(0,len(bits),6))
L=[json.loads(x) for x in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
print(type(L[0]['A']), str(L[0]['A'])[:80])
def toA(a): return [[int(c) for c in r] for r in a] if isinstance(a[0],str) else a
print('d6 consistent with A in all 123:', all(d6(toA(r['A']))==r['d6'] for r in L))
