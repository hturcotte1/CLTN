"""Symbolic (sympy, exact) checks for the Lemma R / e_2 / Lemma AS audit.
No floating point. Re-run: python3 a_symbolic.py"""
import sympy as sp, random, itertools, sys
random.seed(20261014)
beta, t, lam, mu, z, c, e = sp.symbols('beta t lambda mu z c e')
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok &= bool(cond)

# ---- Remark F: K(lam beta, lam t) = (1-lam) I + lam K(beta,t), for a generic N (symbolic 0/1 pattern irrelevant)
m = 4
Nsym = sp.Matrix(m, m, lambda i, j: 0 if i == j else sp.Symbol(f'n{i}{j}'))
I, J = sp.eye(m), sp.ones(m, m)
K = lambda b, tt: I + b * Nsym + tt * J
check("Remark F identity", sp.expand(K(lam*beta, lam*t) - ((1-lam)*I + lam*K(beta, t))) == sp.zeros(m, m))

# ---- IVT root explicitly: f(lam) = lam t (lam beta - 2) - 1; lam* = (1 + sqrt(1 + beta/t))/beta
lst = (1 + sp.sqrt(1 + beta/t))/beta
f = lambda L: L*t*(L*beta - 2) - 1
check("lam* solves f = 0", sp.simplify(f(lst)) == 0)
# lam* < 1 <=> t(beta-2) > 1 (for beta>2, t>0): 1+sqrt(1+beta/t) < beta <=> 1+beta/t < (beta-1)^2 <=> beta/t < beta(beta-2)
check("lam*<1 equivalence algebra", sp.expand((beta-1)**2 - (1 + beta/t) - beta*(beta-2 - 1/t)) == 0)
# lam* > 2/beta always: sqrt(1+beta/t) > 1
# beta' = lam* beta = 1 + sqrt(1+beta/t) > 2 ; t' = lam* t; t'(beta'-2) = 1
bp, tp = sp.simplify(lst*beta), lst*t
check("t'(beta'-2) = 1", sp.simplify(tp*(bp-2) - 1) == 0)

# ---- curve scaling: t*K(2+1/t, t) = tI + (2t+1)N + t^2 J
check("t K(2+1/t,t) = L(t)", sp.expand(t*K(2+1/t, t) - (t*I + (2*t+1)*Nsym + t**2*J)) == sp.zeros(m, m))

# ---- 2x2 principal minors of K = I + beta N + tJ
def minor2(nij, nji):
    M = sp.Matrix([[1+t, beta*nij + t], [beta*nji + t, 1+t]])
    return sp.expand(M.det())
check("bidirected minor = 1+2t", sp.expand(minor2(0, 0) - (1+2*t)) == 0)
check("non-adjacent minor = (1-beta)(1+beta+2t)", sp.expand(minor2(1, 1) - (1-beta)*(1+beta+2*t)) == 0)
check("one-way minor (N_ij=1) = 1 - t(beta-2)", sp.expand(minor2(1, 0) - (1 - t*(beta-2))) == 0)
check("one-way minor (N_ji=1) = 1 - t(beta-2)", sp.expand(minor2(0, 1) - (1 - t*(beta-2))) == 0)

# ---- e_2 rewriting and the inequality chain
B, Nn, O, mm = sp.symbols('B Nn O m')
e2 = B*(1+2*t) - Nn*(beta-1)*(beta+1+2*t) + O*(1 - t*(beta-2))
Phi = mm*(mm-1) - beta*(2*Nn + O)
alt = (t*Phi + B + O - Nn*(beta**2-1)).subs(mm*(mm-1), 2*(B+Nn+O))
check("e_2 = t Phi + B + O - Nn(beta^2-1) (with m(m-1) = 2(B+Nn+O))", sp.expand(e2 - alt) == 0)
# e_2 - (B - Nn)(1+2t) = -Nn[(beta^2-2) + 2t(beta-2)] + O(1 - t(beta-2))
gap = sp.expand(e2 - (B - Nn)*(1+2*t))
check("e_2 - (B-Nn)(1+2t) = -Nn(beta^2-2+2t(beta-2)) + O(1-t(beta-2))",
      sp.expand(gap - (-Nn*(beta**2 - 2 + 2*t*(beta-2)) + O*(1 - t*(beta-2)))) == 0)
# at legal points (t>0, beta>2, t(beta-2)>1) both bracketed terms have the right sign: beta^2-2+2t(beta-2) > 2 > 0 and 1-t(beta-2) < 0,
# so e_2 <= (B-Nn)(1+2t), and e_2 > 0 forces B > Nn.  At m = 7: B+Nn+O = 21, edges = 2B+O = 21 + B - Nn.
for Bv in range(22):
    for Nv in range(22 - Bv):
        Ov = 21 - Bv - Nv
        if (2*Bv + Ov >= 22) != (Bv > Nv): check(f"edges>=22 <=> B>Nn at B={Bv},Nn={Nv}", False)
check("m=7: edges >= 22 <=> B > Nn over all (B,Nn,O) with B+Nn+O=21", True)

# ---- Lemma AS ingredients on random small graphs
def rand_graph(m, p):
    return sp.Matrix(m, m, lambda i, j: 0 if i == j else (1 if random.random() < p else 0))
for trial in range(12):
    m = random.choice([3, 4, 5])
    A = rand_graph(m, random.choice([0.3, 0.5, 0.7]))
    I, J, one = sp.eye(m), sp.ones(m, m), sp.ones(m, 1)
    chi = sp.expand((z*I + A).det())
    psi = sp.expand((one.T * (z*I + A).adjugate() * one)[0])
    lhs = sp.expand((z*I - (c*J - A)).det())
    ok1 = sp.expand(lhs - (chi - c*psi)) == 0
    V = sp.Matrix(m, m-1, lambda i, j: 1 if i == j else (-1 if i == m-1 else 0))      # rational basis of 1-perp
    Pm = V * (V.T*V).inv() * V.T                                                      # orthogonal projector onto 1-perp
    Aperp = (V.T*V).inv() * V.T * A * V                                               # matrix of P A |_(1-perp) in basis V
    ok2 = sp.expand(psi - m*(z*sp.eye(m-1) + Aperp).det()) == 0
    if not (ok1 and ok2): print("AS ingredient fail", A, ok1, ok2)
    ok &= ok1 and ok2
check("det(zI-(cJ-A)) = chi - c psi and psi = m det(zI + A|perp) on 12 random graphs", ok)
# IFT derivative: F(z,e) = psi1 - e chi1, z'(0) = chi1(z0)/psi1'(z0) = 1/s'(z0) where s = psi1/chi1
p1, c1 = sp.Function('p1'), sp.Function('c1')
s_ = p1(z)/c1(z)
P1d, C0, C1d = sp.symbols('P1d C0 C1d')
sd = sp.diff(s_, z).subs({sp.Derivative(p1(z), z): P1d, sp.Derivative(c1(z), z): C1d}).subs({p1(z): 0, c1(z): C0})
check("s'(z0) = psi1'(z0)/chi1(z0) when psi1(z0)=0", sp.simplify(sd - P1d/C0) == 0)
print("ALL PASS" if ok else "SOME FAIL")
