"""INDEPENDENT re-derivation of the curve verdict, not using Lemma B(d)/Orlando or the repo interpolation.
For each graph: c_k(t) = coefficients of det(x I + L(t)) by DomainMatrix Berkowitz over ZZ[t] (symbolic t);
ALL Hurwitz minors Delta_1..Delta_{m-1}(t) by DomainMatrix det over ZZ[t]. By Routh-Hurwitz (b0)+(b1) the
stable set on the curve is exactly {t > 0 : all c_k > 0 and all Delta_k > 0}, which is decided on the cells of
the arrangement of ALL their positive roots. Compared with n7/curve_test.curve_test:
  (1) c_k(t) identical to coeff_polys (interpolation check), Delta_{m-1} identical to fastregion.hurwitz_det_poly;
  (2) verdict nonempty/empty identical;
  (3) every stability change between adjacent fine cells happens at a root of c_m or Delta_{m-1} (Lemma B(d));
  (4) every fine cell lying strictly inside a coarse curve_test cell has the coarse cell's verdict.
usage: python3 e_allminors.py SHARD NSHARDS"""
import sys, json, glob, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rootisolation import dup_isolate_real_roots_list, dup_refine_real_root
from sympy.polys.domains import ZZ
from n7.curve_test import coeff_polys, curve_test, nonedge, from_d6
from fastregion import hurwitz_det_poly
D = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/'
t = sp.Symbol('t')

def poly_lo_hi(expr):
    p = sp.Poly(expr, t, domain=ZZ) if expr != 0 else None
    return [] if p is None else [int(c) for c in p.all_coeffs()[::-1]]

def pev(p, x):
    r = F(0)
    for c in reversed(p): r = r * x + c
    return r

def indep(A):
    m = len(A); N = nonedge(A)
    L = sp.Matrix(m, m, lambda i, j: (t if i == j else 0) + (2 * t + 1) * N[i][j] + t ** 2)
    dM = DomainMatrix.from_Matrix(-L); cp = [dM.domain.to_sympy(c) for c in dM.charpoly()]
    C = [poly_lo_hi(sp.expand(c)) for c in cp]
    Dl = []
    for k in range(1, m):
        H = sp.Matrix(k, k, lambda i, j: cp[2 * (i + 1) - (j + 1)] if 0 <= 2 * (i + 1) - (j + 1) <= m else 0)
        dH = DomainMatrix.from_Matrix(H) if any(isinstance(x, sp.Expr) and x.has(t) for x in H) else None
        val = dH.domain.to_sympy(dH.det()) if dH is not None else H.det()
        Dl.append(poly_lo_hi(sp.expand(val)))
    return C, Dl

def run(A):
    m = len(A); C, Dl = indep(A); P = coeff_polys(A)
    rec = {"coeffs_match": C == P}
    Delta_repo = hurwitz_det_poly(P, m - 1)
    rec["Delta_match"] = (Dl[m - 2] == Delta_repo)
    r = curve_test(A, P); rec["verdict"] = r["verdict"]
    polys = [C[k] for k in range(1, m + 1)] + Dl
    if any(not p for p in polys):
        rec["indep_verdict"] = "never(identically zero)"; rec["ok"] = rec["coeffs_match"] and rec["Delta_match"] and r["verdict"] == "never"; return rec
    nonconst = [p for p in polys if len(p) > 1]
    iso = dup_isolate_real_roots_list([[int(c) for c in reversed(p)] for p in nonconst], ZZ, inf=0, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == 0)]
    bounds = []
    for (a_, b_), idx, f in iso:
        a_, b_ = F(int(a_.numerator), int(a_.denominator)), F(int(b_.numerator), int(b_.denominator))
        while a_ <= 0:
            a_, b_ = dup_refine_real_root(f, a_, b_, ZZ, steps=1); a_, b_ = F(int(a_.numerator), int(a_.denominator)), F(int(b_.numerator), int(b_.denominator))
        bounds.append((a_, b_, f))
    cells = []; prev = F(0)
    for a_, b_, f in bounds:
        if not prev < a_: rec["ok"] = False; rec["err"] = "intervals not separated"; return rec
        cells.append((prev, a_, (prev + a_) / 2)); prev = b_
    cells.append((prev, None, prev + 1))
    stab = []
    for lo, hi, s in cells:
        vals = [pev(p, s) for p in polys]
        assert all(v != 0 for v in vals)          # sample avoids every root
        stab.append(all(v > 0 for v in vals))
    rec["indep_verdict"] = "stable_somewhere" if any(stab) else "never"
    rec["n_fine_cells"] = len(cells)
    # (3) stability changes only at roots of c_m or Delta_{m-1}
    Pcm = sp.Poly(list(reversed(C[m])), t); PD = sp.Poly(list(reversed(Dl[m - 2])), t)
    bad_change = 0
    for i in range(len(bounds)):
        if stab[i] != stab[i + 1]:
            f = sp.Poly(list(bounds[i][2]), t)
            if not (Pcm.rem(f).is_zero or PD.rem(f).is_zero): bad_change += 1
    rec["changes_not_at_cm_or_Delta"] = bad_change
    # (4) fine cells inside coarse cells agree
    mism = 0; attributed = 0
    for (lo, hi, s), st in zip(cells, stab):
        for c in r["cells"]:
            clo = F(c["lo"]); chi_ = None if c["hi"] is None else F(c["hi"])
            if clo <= lo and (chi_ is None or (hi is not None and hi <= chi_)):
                attributed += 1
                if c["stable"] is not None and c["stable"] != st: mism += 1
    rec["fine_attributed"] = attributed; rec["fine_mismatch"] = mism
    rec["ok"] = rec["coeffs_match"] and rec["Delta_match"] and rec["indep_verdict"] == r["verdict"] and bad_change == 0 and mism == 0
    return rec

def main():
    sh, nsh = int(sys.argv[1]), int(sys.argv[2])
    items = []
    for f in sorted(glob.glob(D + 'b_random_shard*.jsonl')):
        for l in open(f):
            x = json.loads(l); items.append(("rand%d" % x["idx"], [[int(ch) for ch in row] for row in x["A"]]))
    items.sort(key=lambda z: int(z[0][4:]))
    for l in open('/home/user/CLTN/n7/n7_limit_HITS.jsonl'):
        x = json.loads(l); items.append(("hit:" + x["d6"], from_d6(x["d6"])))
    for s in ["&FNFL~[udo?", "&FK@dyUfXp?", "&FMEHz\\vpv?", "&FNALh\\Uwp?", "&FM?_~LVd~?"]:
        items.append(("special:" + s, from_d6(s)))
    items.append(("K7", [[0 if i == j else 1 for j in range(7)] for i in range(7)]))
    out = open(D + f'e_allminors_shard{sh}.jsonl', 'w')
    for i in range(sh, len(items), nsh):
        name, A = items[i]; t0 = time.time()
        try:
            rec = run(A)
        except Exception as e:
            rec = {"ok": False, "err": repr(e)}
        rec["name"] = name; rec["sec"] = round(time.time() - t0, 3)
        out.write(json.dumps(rec) + "\n"); out.flush()

main()
