"""Exact anatomy of selected candidates: positivity set on [2, oo) from the polynomials
det(M_b), adj(M_b) 1 with M_b = b A^T + (b - 1) I (w = b adj(M_b)1 / det M_b, Lemma D2 / theory notes 3.1),
the limit polynomial psi_D(z) = chi_D - chi_{D+J}, and the stored stable interval.
Usage: python3 extremal.py d6 [d6 ...]   -> prints a JSON report (decimals marked 'display')."""
import sys, json
import sympy as sp
import lib

b, z = sp.symbols('b z')


def anatomy(d6):
    A = lib.from_d6(d6); n = len(A)
    M = sp.Matrix(n, n, lambda i, j: b * A[j][i] + (b - 1) * (1 if i == j else 0))
    det = sp.expand(M.det(method='berkowitz'))
    adj1 = [sp.expand(x) for x in M.adjugate(method='berkowitz') * sp.ones(n, 1)]
    def roots_above_2(p):
        P = sp.Poly(p, b)
        rs = [r for r in P.real_roots() if r > 2]
        return [str(sp.factor(sp.Poly(r.as_poly() if hasattr(r, 'as_poly') else P, b).as_expr())) for r in []] , rs
    info = {'d6': d6}
    # positivity on (2, oo): sign pattern between consecutive critical points
    crit = set()
    polys = [det] + adj1
    for p in polys:
        for r in sp.Poly(p, b).real_roots():
            if r > 2:
                crit.add(r)
    crit = sorted(crit, key=lambda r: sp.N(r, 50))
    pts = [sp.Integer(2)] + crit
    segs = []
    for k in range(len(pts)):
        lo = pts[k]; hi = pts[k + 1] if k + 1 < len(pts) else None
        # rational sample strictly inside
        if hi is None:
            s = sp.Integer(sp.ceiling(sp.N(lo, 30))) + 1
        else:
            a1, a2 = sp.Rational(str(sp.N(lo, 40))), sp.Rational(str(sp.N(hi, 40)))
            s = (a1 + a2) / 2
            assert lo < s < hi
        d = det.subs(b, s)
        w = [s * x.subs(b, s) / d for x in adj1]
        segs.append(dict(lo_display=str(sp.N(lo, 12)), hi_display=(str(sp.N(hi, 12)) if hi is not None else 'oo'),
                         sample=str(s), positive=all(x > 0 for x in w)))
    info['positivity_segments_on_(2,oo)'] = segs
    info['det_Mb_factored'] = str(sp.factor(det))
    info['adj1_factored'] = [str(sp.factor(x)) for x in adj1]
    D = lib.nonedge(A)
    chiD = sp.Matrix(D).charpoly(z).as_expr()
    chiDJ = sp.Matrix(n, n, lambda i, j: D[i][j] + 1).charpoly(z).as_expr()
    psi = sp.expand(chiD - chiDJ)
    info['psi_D_factored'] = str(sp.factor(psi))
    rts = sp.Poly(psi, z).nroots(n=30)
    minre = min(sp.re(r) for r in rts)
    info['psi_D_min_real_part_display'] = str(sp.N(minre, 10))
    info['beta_c_display = -1/minRe'] = str(sp.N(-1 / minre, 10)) if minre < 0 else 'oo'
    info['D_arcs (i,j): i does not send to j'] = [(i, j) for i in range(n) for j in range(n) if D[i][j]]
    info['outdeg'] = [lib.outdeg(A, i) for i in range(n)]
    info['indeg'] = [lib.indeg(A, i) for i in range(n)]
    info['A_rows'] = [''.join(map(str, r)) for r in A]
    return info


if __name__ == '__main__':
    print(json.dumps([anatomy(d) for d in sys.argv[1:]], indent=1))
