"""End-to-end exact spot checks of the drops against the project's independent verifier ctln.analyze_support
(the CLAUDE.md W convention) and against sympy-computed a_i.  Usage: e2e.py harness_out droppedU_file seed"""
import sys, random
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import ctln
from sympy import symbols, Poly, ZZ
from sympy.polys.matrices import DomainMatrix
b, u = symbols('b u'); R = ZZ[b]
rng = random.Random(int(sys.argv[3]))

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A
def a_polys(A):
    n = len(A)
    Nm = [[1 if (i != j and A[j][i] == 0) else 0 for j in range(n)] for i in range(n)]
    out = []
    for i in range(n):
        rows = [[R(1) if c == i else R.convert(1 if r == c else 0) + R.convert(Nm[r][c]) * R.gens[0] for c in range(n)] for r in range(n)]
        out.append(Poly(R.to_sympy(DomainMatrix(rows, (n, n), R).det()), b))
    return out
def rand_legal():
    k = rng.choice([2, 3, 4, 5, 7, 10, 20, 50, 200])
    eps = F(rng.randint(1, k - 1), k) if rng.random() < 0.5 else F(1, k)
    dmin = eps / (1 - eps)                       # q > 1  <=>  delta > eps/(1-eps)
    r = rng.choice([F(1, 1000), F(1, 50), F(1, 7), F(1, 2), F(1), F(3), F(10)]) * F(rng.randint(1, 9), 5)
    delta = dmin * (1 + r)
    assert ctln.legal(eps, delta)
    return eps, delta

H = [l for l in open(sys.argv[1])]
Ndrops = [l.split()[0] for l in H if ' df=1 ' in l and '#drop=1' in l]
Doms = [l.split()[0] for l in H if ' df=0 ' in l]
Udrops = [l.split()[0] for l in open(sys.argv[2])]
bad = 0; st = {}
def bump(k): st[k] = st.get(k, 0) + 1

# (1) Descartes certificate re-derived with sympy on 3000 real N-dropped graphs
for s in rng.sample(Ndrops, 3000):
    A = from_d6(s); n = len(A)
    ok = False
    for ap in a_polys(A):
        sh = Poly(ap.as_expr().subs(b, 2 + u), u)
        if ap.is_zero or all(c <= 0 for c in sh.all_coeffs()): ok = True; break
    if not ok: print('N-DROP NOT CERTIFIED BY SYMPY', s); bad += 1
    bump('N_sympy')

# (2) ctln at random legal points: dropped graphs never have a positive stable full-support fixed point;
#     x is proportional to a(beta) with sign(x_i) = sign(a_i) * sign(det(I - W))
for kind, pool, cnt in (('N', Ndrops, 700), ('U', Udrops, 700), ('dom', Doms, 400)):
    for s in rng.sample(pool, cnt):
        A = from_d6(s); n = len(A); ap = a_polys(A) if kind != 'dom' else None
        for _ in range(2):
            eps, delta = rand_legal()
            r = ctln.analyze_support(A, list(range(n)), eps, delta)
            bump(kind + '_pts')
            if r['x'] is None: bump(kind + '_singular'); continue
            if r['stable_fixed_point'] or (r['positive'] and r['stable'] is not False):
                print('DROPPED GRAPH WITH POSITIVE NON-UNSTABLE FP', kind, s, eps, delta, r['stable']); bad += 1
            if kind == 'dom' and r['positive']:
                print('DOMINATED GRAPH WITH x > 0', s, eps, delta); bad += 1
            if ap is not None:
                beta = 1 + delta / eps
                av = [F(p.eval(beta)) for p in ap]
                x = [F(v) for v in r['x']]
                # x = (theta/eps) a / D,  D = det(I-W)/eps^n  -> x_i a_j == x_j a_i and same sign pattern up to sign(D)
                for i in range(n):
                    for j in range(n):
                        if x[i] * av[j] != x[j] * av[i]: print('NOT PROPORTIONAL', s); bad += 1; break
                nz = [i for i in range(n) if av[i] != 0]
                if nz:
                    i0 = nz[0]; c = x[i0] / av[i0]
                    if any(x[i] != c * av[i] for i in range(n)): print('SCALE MISMATCH', s); bad += 1
                    bump(kind + '_scale_pos' if c > 0 else kind + '_scale_neg')
print(st, 'BAD', bad)
