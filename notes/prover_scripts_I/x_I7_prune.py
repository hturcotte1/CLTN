"""EXPLORATION (not a proof, not a certificate): how many digraphs does the Lemma I.7 prune remove?
Full support. h = gcd(chi_A, chi_{A+J}); KILL if some root mu of h has Re mu >= -1/2 (decided exactly by
Routh-Hurwitz on each irreducible factor f(z - 1/2); verdict False => killed, None => killed-at-boundary, counted apart).
Also reports overlap with the domination prune (Lemma D.3) and with alpha(A|perp) >= -1/2 (limit test, exact)."""
import sys, random, time
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import sympy as sp
import ctln, graphs, enum_cache
from collections import Counter
z = sp.symbols('z')
def cp(M):  # exact charpoly coefficients [1, c1, ..., cn] of det(z I - M)
    return ctln.charpoly(M)
def poly(coeffs):
    return sp.Poly([sp.Rational(c.numerator, c.denominator) for c in coeffs], z)
def some_root_re_ge_half(f):
    """f in Q[z]; True if some root has Re >= -1/2 (exact), with flag whether only via a vanishing minor."""
    boundary = False
    for g, _ in sp.factor_list(f.as_expr())[1]:
        g = sp.Poly(g, z)
        if g.degree() == 0: continue
        sh = sp.Poly(g.as_expr().subs(z, z - sp.Rational(1, 2)), z)
        co = [F(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in sh.all_coeffs()]
        if co[0] < 0: co = [-c for c in co]
        v, _ = ctln.hurwitz_stable(co)
        if v is False: return True, False
        if v is None: boundary = True
    return boundary, boundary
def analyse(A):
    n = len(A)
    Af = [[F(x) for x in row] for row in A]
    AJ = [[F(A[i][j] + 1) for j in range(n)] for i in range(n)]
    h = sp.gcd(poly(cp(Af)), poly(cp(AJ)))
    chip = (poly(cp(Af)) - poly(cp(AJ))).quo(sp.Poly(n, z))
    kill, bnd = (False, False) if h.degree() == 0 else some_root_re_ge_half(h)
    lim, lbnd = some_root_re_ge_half(chip)
    return kill, bnd, lim, lbnd, h
def run(n, keys):
    c = Counter(); t0 = time.time(); examples = []
    for v in keys:
        A = graphs.from_key(enum_cache.key_from_int(v, n), n)
        if graphs.is_complete_bidirected(A): c['complete'] += 1; continue
        kill, bnd, lim, lbnd, h = analyse(A)
        df = graphs.is_domination_free(A)
        c['total'] += 1; c['domfree'] += df
        c['I7_kill'] += kill; c['I7_kill_boundary_only'] += bnd
        c['I7_kill_and_domfree'] += kill and df
        c['limit_dead(alpha>=-1/2)'] += lim; c['limit_dead_boundary_only'] += lbnd
        c['domfree_and_limit_alive'] += df and not lim
        c['domfree_limit_alive_I7_kill'] += df and (not lim) and kill
        if df and (not lim) and kill and len(examples) < 3: examples.append((A, h.as_expr()))
    print(f"n={n}: {dict(c)}  ({time.time()-t0:.1f}s)")
    for A, h in examples: print("   alive-in-limit but I.7-killed:", A, " h =", sp.factor(h))
keys5 = enum_cache.load_keys(5, '/home/user/CLTN/certificates/digraphs_n5.keys')
run(5, keys5)
if len(sys.argv) > 1:
    keys6 = enum_cache.load_keys(6, '/home/user/CLTN/certificates/digraphs_n6.keys')
    random.seed(1); run(6, random.sample(keys6, int(sys.argv[1])))
