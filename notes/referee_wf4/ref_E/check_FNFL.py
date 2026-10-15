"""Referee: mechanism of &FNFL~[udo? (notes_theory.md §4.1). Exact."""
import sys, json, itertools
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E')
import sympy as sp
from reflib import *
out = {}
A = d6_decode('&FNFL~[udo?'); n = 7
cert = json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
out['matches_FNFL_certificate_rows'] = [''.join(map(str, r)) for r in A] == cert['graph']['adjacency_rows']
# closed twins 2,3: bidirected, same out/in neighbours elsewhere
tw = A[2][3] == 1 and A[3][2] == 1 and all(A[2][k] == A[3][k] and A[k][2] == A[k][3] for k in range(n) if k not in (2, 3))
out['nodes_2_3_closed_twins'] = tw
H3 = json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json'))['A_ours']
def iso(X, Y):
    m = len(X)
    for p in itertools.permutations(range(m)):
        if all(X[i][j] == Y[p[i]][p[j]] for i in range(m) for j in range(m)): return p
    return None
for d in (2, 3):
    keep = [i for i in range(n) if i != d]
    S = [[A[i][j] for j in keep] for i in keep]
    out[f'delete_{d}_iso_to_H3'] = str(iso(S, H3))
# 7 chi_perp factorisation
P = nchi_perp(A)
out['7chi_perp_factor'] = str(sp.factor(P))
out['7chi_perp_claim_ok'] = sp.expand(P - (l+1)*(l**2+l+1)*(7*l**3+14*l**2+14*l+5)) == 0
# node-blow-up family of H3 at the doubled node: map H3 node -> FNFL node 2 via iso of deleting 3
keep = [0, 1, 2, 4, 5, 6]
S = [[A[i][j] for j in keep] for i in keep]   # S = FNFL minus node 3 ; node index 2 in S is the twin node
k = sp.Symbol('k')
m = 6
D = sp.ones(m, m) - sp.eye(m) - sp.Matrix(S)
mu = [1, 1, k+1, 1, 1, 1]
Phimu = sp.expand((sp.Matrix([mu]) * (z*sp.eye(m) - D*sp.diag(*mu)).adjugate() * sp.ones(m, 1))[0])
Phik = sp.expand(z*Phimu)
claim = sp.expand(z*((6+k)*z**5 + (12+2*k)*z**4 + (18+3*k)*z**3 + (13+3*k)*z**2 + (6+3*k)*z + (1+k)))
out['Phi_k_claim_ok'] = sp.expand(Phik - claim) == 0
out['Phi_1_factor'] = str(sp.factor(Phik.subs(k, 1)))
rem = sp.rem(sp.Poly(Phik, z), sp.Poly(z**2+z+1, z))
out['Phi_k_mod_(z^2+z+1)'] = str(rem.as_expr())   # claim: (k-1) z^2 == (k-1)(-z-1)
out['Phi_k_omega_claim_ok'] = sp.expand(rem.as_expr() - sp.rem(sp.Poly((k-1)*z**2, z), sp.Poly(z**2+z+1, z)).as_expr()) == 0
# z psi_D == -z^2 and chi_L == z^2 mod z^2+z+1, with G = H3 (S), class = twin of node 2: b = 1 - S_col? use X1 with S_i = A[i][2] (i->v), T_i = A[2][i] (v->i)
Ssend = [S[i][2] for i in range(m)]; Trecv = [S[2][i] for i in range(m)]
# for a closed twin of node 2, v<->2 bidirected: S_2 = T_2 = 1
Ssend[2] = 1; Trecv[2] = 1
b = sp.Matrix([1 - s for s in Ssend]); cc = sp.Matrix([1 - t for t in Trecv])
psiD = sp.expand((sp.ones(1, m) * (z*sp.eye(m) - D).adjugate() * sp.ones(m, 1))[0])
L = D - b*sp.ones(1, m) - sp.ones(m, 1)*cc.T
chiL = sp.expand((z*sp.eye(m) - L).det())
out['zpsiD_mod'] = str(sp.rem(sp.Poly(z*psiD, z), sp.Poly(z**2+z+1, z)).as_expr())
out['chiL_mod'] = str(sp.rem(sp.Poly(chiL, z), sp.Poly(z**2+z+1, z)).as_expr())
out['X1_consistency'] = sp.expand(z*psiD + k*chiL - Phik) == 0
# no 3-cycle module: a 3-set {x,y,w} inducing a directed 3-cycle (one-way) with identical relations to every outside node
mods = []
for T in itertools.combinations(range(n), 3):
    rest = [v for v in range(n) if v not in T]
    if all(all(A[u][v] == A[T[0]][v] and A[v][u] == A[v][T[0]] for u in T) for v in rest):
        sub = [[A[i][j] for j in T] for i in T]
        e = sum(sub[i][j] for i in range(3) for j in range(3))
        oneway = all(sub[i][j] + sub[j][i] == 1 for i in range(3) for j in range(3) if i < j)
        if oneway and e == 3: mods.append(T)
out['three_cycle_modules'] = mods
print(json.dumps(out, indent=1))
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E/check_FNFL.json', 'w'), indent=1)
