"""
Tk_table.py -- EXACT stable t-sets of the family members at a fixed beta (question (a) of the task).
For a family and integer k, Psi_k(z) = z det(zI - D - aJ) - a k sum rho_j chi_{L_j}(z) is LINEAR in a = t/beta, so the Hurwitz
minors of (-1)^deg Psi_k(-1/beta - x) are polynomials Delta_j(a).  The stable a-set is {a > a1 : all Delta_j(a) > 0}
(a1 = legal boundary); it is computed exactly from the real roots of the Delta_j (sympy real_roots: exact isolation) and one
exact rational sample per cell.  T_k(beta) := beta * inf{a0 : stable for every a > a0} (inf over the legal range; 'inf' if the
member is unstable for all large t).  Output: list of stable t-intervals (endpoints algebraic; decimals for display) and T_k.
Writes Tk_table.json.
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

G1 = [[0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0]]
FAM = dict(FAMILIES); FAM['F'] = (None, [(1, 1, 1, 1, 1, 0)], [1])


def stable_a_set(g, beta, kk):
    r = 1 / sp.Rational(beta); a1 = fl.legal_a_bound(r)
    gk = sp.expand(g.subs(fl.k, kk))
    co, mins = fl.region_minors(gk, r)
    crit = set()
    for mn in mins:
        P = sp.Poly(mn, fl.a_s)
        if P.degree() > 0:
            for rt in P.real_roots():
                if rt > a1: crit.add(rt)
    crit = sorted(crit, key=lambda v: sp.N(v, 60))
    pts = [a1] + crit
    cells = []
    for i in range(len(pts)):
        lo = pts[i]; hi = pts[i + 1] if i + 1 < len(pts) else None
        # exact rational sample strictly inside the cell
        if hi is None:
            smp = sp.Rational(int(sp.ceiling(sp.N(lo, 30))) + 1)
        else:
            smp = sp.nsimplify((sp.N(lo, 40) + sp.N(hi, 40)) / 2, rational=True, tolerance=sp.Rational(1, 10 ** 35))
            assert lo < smp < hi
        vals = [mn.subs(fl.a_s, smp) for mn in mins]
        ok = all(v > 0 for v in vals)
        cells.append((lo, hi, ok))
    # T_k: start of the final run of stable cells reaching infinity
    if not cells[-1][2]:
        Tk = None
    else:
        i = len(cells) - 1
        while i > 0 and cells[i - 1][2]: i -= 1
        Tk = cells[i][0] * sp.Rational(beta)
    return cells, Tk


def run(nm, beta, ks):
    base, classes, rhos = FAM[nm]
    A = G1 if nm == 'F' else fl.from_d6(base)
    g, j0 = fl.reduced_psi(A, classes, rhos)
    rows = []
    for kk in ks:
        cells, Tk = stable_a_set(g, beta, kk)
        st = [[str(sp.N(lo * sp.Rational(beta), 10)), (str(sp.N(hi * sp.Rational(beta), 10)) if hi is not None else 'inf')]
              for lo, hi, ok in cells if ok]
        rows.append({'k': kk, 'stable_t_intervals_display': st,
                     'T_k': ('inf (unstable for all large t)' if Tk is None else str(sp.N(Tk, 12))),
                     'T_k_is_legal_boundary': (Tk is not None and Tk == fl.legal_a_bound(1 / sp.Rational(beta)) * sp.Rational(beta))})
        print(nm, beta, kk, rows[-1]['T_k'], st, flush=True)
    return {'family': nm, 'beta': str(beta), 'legal_boundary_t': str(1 / (sp.Rational(beta) - 2)), 'rows': rows}


if __name__ == '__main__':
    out = []
    out.append(run('Q', '5/2', [0, 5, 10, 20, 24, 25, 26, 30, 35, 39, 40, 41, 50, 100]))
    out.append(run('R', '5/2', [0, 5, 10, 16, 17, 20, 25, 30, 33, 34, 50, 100]))
    out.append(run('F', '21/10', [0, 5, 10, 20, 30, 35, 38, 39, 40, 41, 45, 60]))
    json.dump(out, open('Tk_table.json', 'w'), indent=1)
