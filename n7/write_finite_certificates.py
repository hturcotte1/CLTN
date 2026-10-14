"""
write_finite_certificates.py -- exact certificates for finite-eps phenomena at n = 7 found by the
Lemma-M workflow (2026-10-14). No floating point.

1. certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json
   graph &FNFL~[udo? (nauty labelg canonical &FS@{`]f|~?): domination-free, 28 edges, NOT a clique,
   spectral abscissa alpha(A|perp) = -1/2 EXACTLY (filter7 'B'; 7 chi_perp = (l+1)(l^2+l+1)(7l^3+14l^2+14l+5)),
   hence never limit-stable (Lemma I.3) and not among the 123 limit candidates -- yet it carries a stable
   positive full-support fixed point at legal finite (eps, delta). Full two-parameter region (region2d) +
   pointwise exact witnesses.
2. certificates/n7_limit_hits_finite_eps_bulge_witnesses.json
   two of the 123 limit candidates whose finite-eps stable region extends to beta > beta_c (bounded above in t):
   pointwise exact witnesses outside the limit interval + region2d.
"""
import json, sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ctln import analyze_support, legal, is_clique
from graphs import is_domination_free, canonical_form, is_complete_bidirected
from fastregion import region, INF, beta_data
from region2d import region2d


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def witness(A, eps, delta):
    r = analyze_support(A, list(range(len(A))), F(eps), F(delta))
    t = (1 - F(eps)) / F(eps); beta = 1 + F(delta) / F(eps); q = F(delta) * (1 - F(eps)) / F(eps)
    return {"eps": str(F(eps)), "delta": str(F(delta)), "t": str(t), "beta": str(beta), "q": str(q),
            "legal": bool(legal(F(eps), F(delta))) if callable(legal) else bool(q > 1),
            "is_fixed_point": bool(r["is_fixed_point"]), "x": [str(v) for v in r["x"]] if r.get("x") is not None else None,
            "positive": bool(r.get("x") is not None and all(F(str(v)) > 0 for v in r["x"])),
            "stable": r["stable"], "stable_detail": r.get("stability_detail"), "stable_fixed_point": bool(r.get("stable_fixed_point")),
            "target_free_clique": bool(r.get("target_free")), "charpoly": [str(c) for c in r["charpoly"]] if r.get("charpoly") else None}


def graph_block(d6):
    A = from_d6(d6)
    bd = beta_data(A, list(range(7)))
    lim = region(A, list(range(7)), INF, bd=bd, replay=True)
    r2 = region2d(A, list(range(7)), crosscheck=True, replay=True)
    return A, {"digraph6_as_enumerated": d6, "adjacency_rows": ["".join(map(str, row)) for row in A],
               "canonical_key_int": None,
               "edges": sum(map(sum, A)), "clique": bool(is_clique(A, list(range(7)))), "domination_free": bool(is_domination_free(A)),
               "limit_region": {"verdict": lim["verdict"], "stable_intervals": lim.get("stable_intervals"), "reason": lim.get("reason")},
               "region2d": r2}


def main():
    out = {}
    # 1. the new finite-eps-only counterexample
    d6 = "&FNFL~[udo?"
    A, blk = graph_block(d6)
    blk["labelg_canonical"] = "&FS@{`]f|~?"
    blk["alpha_statement"] = "alpha(A|perp) = -1/2 exactly: 7*chi_perp(l) = (l+1)(l^2+l+1)(7l^3+14l^2+14l+5); the primitive cube roots of unity give real part -1/2, every other root has real part < -1/2 (filter7 flag 'B'; exact factorisation by the workflow refuter, re-verified by the referee)."
    blk["witnesses"] = {
        "stable_positive_legal": [witness(A, F(1, 49), F(33, 1568)), witness(A, F(1, 36), F(17, 588))],
        "unstable_same_beta_larger_t": [witness(A, F(1, 101), F(33, 3232)), witness(A, F(1, 1001), F(33, 32032)), witness(A, F(1, 61), F(51, 2989))],
    }
    out["FNFL"] = blk
    json.dump({"claim": "non-clique full support of a 7-node digraph is a stable positive fixed point at legal finite (eps, delta); the graph is never limit-stable (alpha = -1/2 exactly), so it is invisible to the eps -> 0 classification; its stable region is bounded above in t",
               "status": "CANDIDATE -- awaiting human hand-verification (CLAUDE.md section 3); nothing is built on it",
               "found_by": "workflow wf_945cf15e-c68, two independent refuters (float search heuristics, exact certification with ctln.analyze_support), each confirmed by an independent exact referee",
               "conventions": "CLAUDE.md section 2; A[i][j]=1 iff i->j; t=(1-eps)/eps; beta=1+delta/eps; legal iff t(beta-2)>1",
               "graph": blk}, open("certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json", "w"), indent=1)
    print("FNFL:", "region2d", blk["region2d"]["verdict"], "stable_cells", blk["region2d"]["stable_cells"], "| limit", blk["limit_region"]["verdict"])
    for w in blk["witnesses"]["stable_positive_legal"]: print("  S", w["eps"], w["delta"], "legal", w["legal"], "pos", w["positive"], "stable", w["stable"])
    for w in blk["witnesses"]["unstable_same_beta_larger_t"]: print("  U", w["eps"], w["delta"], "legal", w["legal"], "pos", w["positive"], "stable", w["stable"], w["stable_detail"])
    # 2. bulge witnesses on two limit candidates
    bul = {}
    for d6, pts_s, pts_u in [("&FNDL}[u|o?", [(F(1, 21), F(9, 175))], [(F(1, 61), F(27, 1525)), (F(1, 100001), F(27, 2500025))]),
                             ("&FK@dyUfXp?", [(F(1, 11), F(5561, 55000))], [(F(1, 1001), F(5561, 5005000)), (F(1, 10001), F(5561, 50005000))])]:
        A, blk = graph_block(d6)
        blk["witnesses"] = {"stable_positive_legal_beyond_limit_interval": [witness(A, e, d) for e, d in pts_s],
                            "unstable_same_beta_larger_t": [witness(A, e, d) for e, d in pts_u]}
        bul[d6] = blk
        print(d6, "limit", blk["limit_region"]["stable_intervals"], "| region2d", blk["region2d"]["verdict"], blk["region2d"]["stable_cells"])
        for w in blk["witnesses"]["stable_positive_legal_beyond_limit_interval"]: print("  S beta", w["beta"], "t", w["t"], "legal", w["legal"], "pos", w["positive"], "stable", w["stable"])
        for w in blk["witnesses"]["unstable_same_beta_larger_t"]: print("  U beta", w["beta"], "t", w["t"], "stable", w["stable"], w["stable_detail"])
    json.dump({"claim": "for these limit candidates the finite-eps stable region extends to beta above beta_c and is bounded above in t at such beta (Lemma AS: first-order term R* > 0); the t -> infinity classification understates their stable region",
               "status": "witnesses on EXISTING candidates (already awaiting hand-verification); descriptive only",
               "graphs": bul}, open("certificates/n7_limit_hits_finite_eps_bulge_witnesses.json", "w"), indent=1)


if __name__ == "__main__":
    main()
