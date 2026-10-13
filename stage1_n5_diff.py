"""
stage1_n5_diff.py -- Stage 1: reproduce the published n=5 classification.

SOURCE  github.com/ccurto/n5-graphs-package, commit 6a24cf5e23f3632906dfb16b71809428c3f9a83d,
        file n5_digraphs.mat. Its sha256 is computed at export time and written
        into the header of data/n5_package_FP.csv (never hand-typed).
CONVENTION  the package stores sA with  sA(i,j) = 1  iff  j -> i
        (graph2net.m: W(i,j) = -1+e iff sA(i,j)=1, no transpose; plot_graph.m:
        "arrow from node j to node i if sA(i,j)=1").  Ours is A[i][j]=1 iff i->j,
        so  A = sA^T.  The graph SET is closed under reversal, so only the
        per-graph fixed-point data can detect a convention error.
DATA    FPcell{g} = every fixed-point support (stable or not) at the package
        defaults (graph2net.m: e = 0.25, d = 2e = 0.5, theta = 1), i.e. q = 3/2.
        The package records NO stability data (check_stability.m is not shipped).
        Its off-support test is non-strict (<= 0); ours is strict (< 0); exact
        ties are counted separately below.

Usage
  python3 stage1_n5_diff.py --mat <path>/n5_digraphs.mat   export data/n5_package_FP.csv, then diff
  python3 stage1_n5_diff.py                                diff against the committed CSV
Writes certificates/stage1_n5_diff.json.  The diff must be ZERO.
"""
import argparse, csv, hashlib, json, os, re, sys, time
from fractions import Fraction as F
from itertools import combinations
from multiprocessing import Pool

from ctln import analyze_support
from graphs import canonical_form, enumerate_digraphs, A000273

CSV = os.path.join("data", "n5_package_FP.csv")
EPS, DELTA = F(1, 4), F(1, 2)          # package defaults, q = 3/2
N = 5


def parse_supports(s):
    """'{ [1  2  3], [4] }' -> [(0,1,2), (3,)]  (0-indexed, sorted).  '{ [] }' -> []."""
    out = []
    for grp in re.findall(r"\[([^\]]*)\]", s):
        nodes = tuple(sorted(int(v) - 1 for v in grp.split()))
        if nodes:
            out.append(nodes)
    return sorted(out)


def export(mat_path):
    import numpy as np
    from scipy.io import loadmat
    sha = hashlib.sha256(open(mat_path, "rb").read()).hexdigest()
    M = loadmat(mat_path)
    sAcell, FPcell, labels = M["sAcell"], M["FPcell"], M["sA_label"]
    assert sAcell.shape == (1, 9608), sAcell.shape
    os.makedirs("data", exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        fh.write(f"# source: ccurto/n5-graphs-package n5_digraphs.mat sha256={sha}\n")
        fh.write("# A_rowmajor is OUR convention (A[i][j]=1 iff i->j) = transpose of the package's sA\n")
        w = csv.writer(fh)
        w.writerow(["graph", "label", "A_rowmajor", "FP_supports_1indexed"])
        for g in range(9608):
            sA = np.asarray(sAcell[0, g]).astype(int)
            assert sA.shape == (5, 5) and set(sA.flatten()) <= {0, 1} and all(sA[i, i] == 0 for i in range(5))
            A = sA.T
            fp = str(FPcell[0, g].flat[0]) if FPcell[0, g].size else "{ [] }"
            lab = str(labels[0, g].flat[0]) if labels[0, g].size else ""
            sup = parse_supports(fp)
            w.writerow([g + 1, lab, "".join(str(A[i, j]) for i in range(5) for j in range(5)),
                        ";".join(" ".join(str(v + 1) for v in s) for s in sup)])
    print(f"exported {CSV}  (mat sha256 {sha})")


def load_csv():
    rows = []
    with open(CSV) as fh:
        header = [fh.readline(), fh.readline()]
        for r in csv.DictReader(fh):
            A = [[int(r["A_rowmajor"][5 * i + j]) for j in range(5)] for i in range(5)]
            sup = sorted(tuple(sorted(int(v) - 1 for v in s.split())) for s in r["FP_supports_1indexed"].split(";") if s)
            rows.append((int(r["graph"]), r["label"], A, sup))
    return header, rows


def analyze_graph(item):
    g, A, sup_pkg = item
    fps, stable, viol, boundary, ties = [], [], [], 0, 0
    for size in range(1, N + 1):
        for sigma in combinations(range(N), size):
            rec = analyze_support(A, list(sigma), EPS, DELTA)
            if rec.get("BOUNDARY"):
                boundary += 1
            if rec["is_fixed_point"]:
                fps.append(sigma)
            elif rec["x"] is not None and rec["positive"] and not rec["off_support_ok"]:
                # could be a tie (drive == 0 exactly): the package counts <= 0 as a fixed point
                from ctln import weight_matrix
                W = weight_matrix(A, EPS, DELTA)
                x = [F(v) for v in rec["x"]]
                if all(F(1) + sum(W[k][sigma[j]] * x[j] for j in range(len(sigma))) <= 0
                       for k in range(N) if k not in sigma):
                    ties += 1
            if rec.get("stable_fixed_point"):
                stable.append(sigma)
            if rec["agrees"] is False:
                viol.append(sigma)
    return g, sorted(fps) == sup_pkg, len(fps), stable, viol, boundary, ties, sorted(fps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", default=None)
    ap.add_argument("--procs", type=int, default=3)
    args = ap.parse_args()
    if args.mat:
        export(args.mat)
    header, rows = load_csv()
    assert len(rows) == 9608, len(rows)

    # graph-set check against our enumerator
    keys = {canonical_form(A)[0] for _, _, A, _ in rows}
    ours = enumerate_digraphs(5)
    assert len(keys) == 9608 == A000273[5], "package graphs are not 9608 distinct isomorphism classes"
    assert keys == ours, "package graph set differs from our enumeration"
    print("graph set: 9608 distinct classes, identical to enumerate_digraphs(5)")

    t0 = time.time()
    with Pool(args.procs) as pool:
        res = pool.map(analyze_graph, [(g, A, sup) for g, _, A, sup in rows], chunksize=50)
    mism = [(g, ours_fps) for g, ok, _, _, _, _, _, ours_fps in res if not ok]
    n_fp = sum(r[2] for r in res)
    n_stable = sum(len(r[3]) for r in res)
    viol = [(g, v) for g, _, _, _, v, _, _, _ in res if v]
    n_boundary = sum(r[5] for r in res)
    n_ties = sum(r[6] for r in res)
    cert = {"stage": 1, "source": "ccurto/n5-graphs-package n5_digraphs.mat", "csv_header": header,
            "eps": str(EPS), "delta": str(DELTA), "q": str(DELTA * (1 - EPS) / EPS),
            "graphs": 9608, "graph_set_matches_A000273_enumeration": True,
            "fp_support_mismatches": len(mism), "mismatch_graphs": mism[:50],
            "total_fp_supports": n_fp, "total_stable_fp_supports": n_stable,
            "conjecture_violations": viol, "boundary_cases": n_boundary, "off_support_ties": n_ties,
            "seconds": round(time.time() - t0, 1)}
    os.makedirs("certificates", exist_ok=True)
    json.dump(cert, open(os.path.join("certificates", "stage1_n5_diff.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in cert.items() if k not in ("csv_header", "mismatch_graphs")}, indent=1))
    if mism or viol:
        print("*** STAGE 1 DID NOT REPRODUCE -- STOP ***")
        sys.exit(1)
    print("STAGE 1 REPRODUCED: diff of zero.")


if __name__ == "__main__":
    main()
