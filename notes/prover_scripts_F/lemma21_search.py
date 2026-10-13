"""Search for N (0/1, zero diagonal) and beta > 2 with adj(C)1 of common strict sign s and
sign(det C) != s, C = I + beta N.  Then K = C + tJ is permitted for large t (sign D = s) while C is not.
Exact: rational beta grid + exact cell sampling via sympy isolating intervals."""
import sys, random
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
import sympy as sp
import fastregion as fr
import graphs

b = sp.symbols("b")

def poly(p):  # low->high int list -> sympy Poly
    return sp.Poly(list(reversed([int(c) for c in p])) if p else [0], b)

def cells_samples(polys):
    """rational sample points in every open cell of (2, inf) cut by real roots of polys, plus roots excluded"""
    prod = sp.Poly(1, b)
    for p in polys:
        if not p.is_zero and p.degree() > 0:
            prod = prod * p
    pts = []
    if prod.degree() > 0:
        sq = sp.Poly(sp.sqf_part(prod.as_expr()), b)
        ivs = sq.intervals(eps=sp.Rational(1, 10**6))
        ends = sorted([(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), _ in ivs])
        rs = [e for e in ends if e[1] > 2]
        cuts = [sp.Rational(2)]
        for lo, hi in rs:
            cuts.append(max(lo, sp.Rational(2))); cuts.append(hi)
        cuts.append(max(cuts[-1] + 10, sp.Rational(1000)))
        for i in range(0, len(cuts) - 1, 2):
            lo, hi = cuts[i], cuts[i + 1]
            if hi > lo:
                pts.append((lo + hi) / 2)
    else:
        pts = [sp.Rational(3)]
    return pts

def scan(A):
    n = len(A)
    bd = fr.beta_data(A, list(range(n)))
    d = poly(bd["d"]); a = [poly(ai) for ai in bd["a"]]
    found = []
    for b0 in cells_samples([d] + a):
        dv = d.eval(b0); av = [ai.eval(b0) for ai in a]
        if all(v > 0 for v in av) and dv <= 0:
            found.append(("a>0,d<=0", b0, dv))
        if all(v < 0 for v in av) and dv >= 0:
            found.append(("a<0,d>=0", b0, dv))
    return found

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    tot = 0; bad = []
    for n in range(2, nmax + 1):
        for key in sorted(graphs.enumerate_digraphs(n)):
            A = graphs.from_key(key, n)
            tot += 1
            f = scan(A)
            if f:
                bad.append((n, A, f))
        print("n =", n, "graphs so far", tot, "hits so far", len(bad), flush=True)
    for x in bad[:20]:
        print(x)
