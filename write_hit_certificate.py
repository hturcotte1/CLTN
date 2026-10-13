"""
write_hit_certificate.py --n 6 --idx <graph index>
Writes certificates/COUNTEREXAMPLE_CANDIDATE_n<n>_idx<idx>.json for a graph that
the limit sweep flagged as stable_somewhere: the limit certificate from the HITS
file plus exact pointwise witnesses from ctln.analyze_support at a rational beta
strictly inside the first stable cell and several large t. Mechanical; no
judgment. Status is always 'awaiting human hand-verification'.
"""
import argparse, glob, json
from fractions import Fraction as F

from enum_cache import key_from_int
from graphs import from_key, is_domination_free
from ctln import analyze_support, legal, is_clique, targets
from aggregate import features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--idx", type=int, required=True)
    args = ap.parse_args()
    n = args.n
    hit = None
    for p in glob.glob(f"certificates/region_n{n}_limit_shard*_HITS.jsonl"):
        for line in open(p):
            r = json.loads(line)
            if r["idx"] == args.idx:
                hit = r
    assert hit, "index not in any HITS file"
    c = hit["cert"]
    A = from_key(key_from_int(hit["key"], n), n)
    f = features(A)
    stable_cells = [cell for cell in c["cells"] if cell["positive"] and cell["off_ok"] and cell["stable"] is True]
    beta0 = F(stable_cells[0]["sample_beta"])
    witness = []
    for t in (F(100), F(1000), F(10000), F(100000)):
        eps = 1 / (t + 1); delta = (beta0 - 1) * eps
        rec = analyze_support(A, list(range(n)), eps, delta)
        witness.append({"beta": str(beta0), "t": str(t), "eps": str(eps), "delta": str(delta), "q": rec["q"],
                        "legal": legal(eps, delta), "is_fixed_point": rec["is_fixed_point"], "stable": rec["stable"],
                        "stability_detail": rec["stability_detail"], "stable_fixed_point": rec["stable_fixed_point"],
                        "target_free": rec["target_free"], "agrees_with_conjecture": rec["agrees"], "x": rec["x"], "charpoly": rec["charpoly"]})
    cert = {"status": "CLAIMED COUNTEREXAMPLE -- AWAITING HUMAN HAND-VERIFICATION (CLAUDE.md section 3). Nothing has been built on it.",
            "n": n, "idx": hit["idx"], "canonical_key": hit["key"], "A_ours": A, "support": list(range(n)),
            "why_not_target_free_clique": {"clique": is_clique(A, range(n)), "targets": targets(A, range(n))},
            "features": f, "limit_certificate": c, "stable_q_limit": hit["stable_q"],
            "pointwise_witnesses_ctln_analyze_support": witness,
            "how_to_rerun": "fastregion.region(A, range(n), fastregion.INF) reproduces limit_certificate; ctln.analyze_support(A, range(n), eps, delta) with each witness's eps, delta reproduces that witness",
            "isomorphic_to_geneson": f["isomorphic_to_geneson"]}
    path = f"certificates/COUNTEREXAMPLE_CANDIDATE_n{n}_idx{hit['idx']}.json"
    json.dump(cert, open(path, "w"), indent=1)
    print(path)
    print("A =", A)
    print("edges", f["edges"], "bidirected", f["bidirected_pairs"], "mutually non-adjacent", f["mutually_nonadjacent_pairs"],
          "|Aut|", f["automorphisms"], "domfree", f["domination_free"], "iso Geneson", f["isomorphic_to_geneson"])
    print("limit stable q:", hit["stable_q"])
    print("compressed charpoly:", f["compressed_charpoly_hi_to_lo"])
    for w in witness:
        print(f"  t={w['t']}: legal={w['legal']} q={w['q']} fp={w['is_fixed_point']} stable={w['stable']} agrees={w['agrees_with_conjecture']}")


if __name__ == "__main__":
    main()
