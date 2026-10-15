"""Referee re-derivation of Theorem E closed forms, symbolic in a. Exact; no floats in decisions."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E')
import sympy as sp
from reflib import *

out = {}
H = d6_decode('&FNDL~[uds?')
H_notes = [[0,0,1,1,1,1,0],[0,0,1,0,1,0,0],[1,1,0,1,1,1,1],[1,1,1,0,1,1,1],[0,0,1,1,0,1,1],[0,1,0,0,1,0,1],[1,1,0,1,0,0,0]]
assert H == H_notes; out['H_matches_d6'] = True
m = 7
mu = [1, 1, 2, a, 1, 1, 1]
M = sp.diag(*mu)
A = sp.Matrix(H)
D = sp.ones(m, m) - sp.eye(m) - A
# --- (E2) D picture: Phi = mu^T adj(zI - D M) 1
X = z*sp.eye(m) - D*M
Phi = sp.expand((sp.Matrix([mu]) * X.adjugate() * sp.ones(m, 1))[0])
target = sp.expand(z*(z**2+z+1)*((a+7)*z*(z**2+z+1) + 2))
out['E2_Phi_minus_target'] = str(sp.expand(Phi - target)); assert sp.expand(Phi - target) == 0
out['E2_Phi_degree_in_a'] = sp.Poly(Phi, a).degree()
# --- A picture (W2): s(z) = mu^T (zI + A_Q)^{-1} 1, A_Q = (A+I)M - I ;  P = mu^T adj(zI+A_Q) 1 , Q = det(zI + A_Q)
AQ = (A + sp.eye(m))*M - sp.eye(m)
XA = z*sp.eye(m) + AQ
Q = sp.expand(XA.det(method='berkowitz'))
P = sp.expand((sp.Matrix([mu]) * XA.adjugate() * sp.ones(m, 1))[0])
# consistency: psi_A(z) = psi_D(z-1)  =>  P(z) = Phi(z-1)
assert sp.expand(P - Phi.subs(z, z-1)) == 0; out['P_equals_Phi_shifted'] = True
# --- (E3) cubic c z(z^2+z+1)+2 roots Re > -1/2 : y = -(z + 1/2)
c, y = sp.symbols('c y')
cub = c*z*(z**2+z+1) + 2
h = sp.Poly(sp.expand(-cub.subs(z, -y - sp.Rational(1, 2))), y)
co = h.all_coeffs()
D2 = sp.expand(co[1]*co[2] - co[0]*co[3])
out['E3_coeffs'] = [str(x) for x in co]; out['E3_Delta2'] = str(D2)
# Hurwitz cubic: co0>0, co1>0, co3>0, Delta2>0 (Delta3 = co3*Delta2). c = a+7 >= 7 > 16/3.
assert co[0] == c and sp.expand(co[1] - c/2) == 0 and sp.expand(co[3] - (3*c/8 - 2)) == 0 and sp.expand(D2 - 2*c) == 0
# --- (E4) chi(z0) != 0 : Q mod z^2 - z + 1, as element of Q(a)(z0)
Qr = red(Q, z); out['E4_Q_mod_(z^2-z+1)'] = [str(Qr[0]), str(Qr[1])]
# P(z0) = 0 identically in a
Pr = red(P, z); out['E4_P_mod'] = [str(Pr[0]), str(Pr[1])]; assert Pr == (0, 0)
# simplicity: P'(z0) != 0
dP = red(sp.diff(P, z), z); out['E4_Pprime_mod'] = [str(dP[0]), str(dP[1])]
# --- (E5) 1/s'(z0) = Q(z0)/P'(z0) (since P(z0)=0)
val = qmul(Qr, qinv(dP))
out['E5_inv_sprime'] = [str(sp.factor(val[0])), str(sp.factor(val[1]))]
R = qre(val); out['E5_R'] = str(sp.factor(R))
# --- (E6) positivity at gamma = 1/2:  ((A^T+I)M - I/2) y = 1
Mpos = (A.T + sp.eye(m))*M - sp.Rational(1, 2)*sp.eye(m)
yv = [sp.factor(sp.simplify(v)) for v in Mpos.LUsolve(sp.ones(m, 1))]
out['E6_y'] = [str(v) for v in yv]
# check against claimed
claimed = [sp.Rational(k)/(18*a+17) for k in (6,2,2,18,6,2,10)]
assert all(sp.simplify(yv[i] - claimed[i]) == 0 for i in range(m)); out['E6_matches_claim'] = True
detMpos = sp.factor(Mpos.det()); out['E6_det'] = str(detMpos)
# W' = sum mu_i y_i (total) > 1 ?
Wtot = sp.factor(sum(mu[i]*yv[i] for i in range(m))); out['E6_W_total'] = str(Wtot)
print(json.dumps(out, indent=1))
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E/check_symbolic.json', 'w'), indent=1)
