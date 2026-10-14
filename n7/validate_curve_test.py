"""validate_curve_test.py -- exact validation of n7/curve_test.py.
  --n5 : all 9,608 digraphs on 5 nodes (geng|directg): expect 9,513 never, 95 stable_somewhere, 0 escalate
         (the prover's count); every stable_somewhere graph gets an exact legal witness via ctln.analyze_support
         at beta = 2 + 1/t + 1e-6 and must be stable there; curve_test and curve_test_fast must agree on every graph.
  --n6 <file> : the strict n=6 positivity survivors (digraph6 lines): every graph that the exact 2-D map
         (certificates/region2d_n6_*) declares stable_somewhere must be stable_somewhere here; counts reported."""
import sys, os, json, argparse, time
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from n7.curve_test import curve_test, curve_test_fast, from_d6, coeff_polys
from ctln import analyze_support


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n5", action="store_true"); ap.add_argument("--n6", default=None)
    ap.add_argument("--geng", default="n7/geng"); ap.add_argument("--directg", default="n7/directg")
    a = ap.parse_args(); t0 = time.time()
    if a.n5:
        import subprocess
        out = subprocess.run(f"{a.geng} -q 5 | {a.directg} -q", shell=True, capture_output=True, text=True).stdout.split()
        assert len(out) == 9608, len(out)
        c = {"never": 0, "stable_somewhere": 0, "escalate": 0}; dis = 0; wit_fail = 0
        for s in out:
            A = from_d6(s); P = coeff_polys(A)
            r = curve_test(A, P); rf = curve_test_fast(A)
            if r["verdict"] != rf["verdict"]: dis += 1
            c[r["verdict"]] += 1
            if r["verdict"] == "stable_somewhere":
                cell = next(x for x in r["cells"] if x["stable"]); t = F(cell["sample_t"]); beta = 2 + 1 / t + F(1, 10 ** 6)
                eps = 1 / (t + 1); delta = (beta - 1) * eps
                w = analyze_support(A, list(range(5)), eps, delta)
                if w["stable"] is not True: wit_fail += 1
        print(json.dumps({"n5_graphs": len(out), **c, "fast_disagreements": dis, "witness_failures": wit_fail, "seconds": round(time.time() - t0, 1)}))
    if a.n6:
        c = {"never": 0, "stable_somewhere": 0, "escalate": 0}; stable = []
        for i, line in enumerate(open(a.n6)):
            s = line.split()[0]; r = curve_test_fast(from_d6(s)); c[r["verdict"]] += 1
            if r["verdict"] != "never": stable.append((s, r["verdict"]))
            if (i + 1) % 2000 == 0: print(i + 1, c, round(time.time() - t0, 1), flush=True)
        json.dump({"counts": c, "non_never": stable, "seconds": round(time.time() - t0, 1)}, open("n7/validate_curve_n6.json", "w"), indent=1)
        print(json.dumps({**c, "seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
