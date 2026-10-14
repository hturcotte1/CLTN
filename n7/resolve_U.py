"""
resolve_U.py -- exact Python resolution of the 'U' (undecided) lines of filter7p.

For each input line "<digraph6> <edges> <P|U> <twin>":
  P  -> kept unchanged (positivity somewhere certified in C by exact evaluation);
  U  -> fastregion.region(A, V, INF, positivity_only=True, replay=True) decides
        exactly whether some beta > 2 has all a_i(beta) > 0 (Lemma A.2 condition);
        kept as 'R' if yes, dropped (counted) if no.
Output: resolved survivor file (same 4-column format, flags P/R) and a stats json.
No floating point anywhere.
"""
import sys, json, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fastregion import region, INF


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dropped", default=None)
    a = ap.parse_args()
    t0 = time.time()
    st = {"P": 0, "U": 0, "U_resolved_positive": 0, "U_resolved_never": 0, "other_verdict": 0}
    drop = open(a.dropped, "w") if a.dropped else None
    with open(a.out, "w") as out:
        for path in a.inputs:
            for line in open(path):
                parts = line.split()
                if not parts:
                    continue
                s, e, f, tw = parts[0], parts[1], parts[2], parts[3] if len(parts) > 3 else "0"
                if f == "P":
                    st["P"] += 1; out.write(f"{s} {e} P {tw}\n"); continue
                assert f == "U", line
                st["U"] += 1
                A = from_d6(s); n = len(A)
                c = region(A, list(range(n)), INF, positivity_only=True, replay=True)
                if c["verdict"] == "positivity_possible":
                    st["U_resolved_positive"] += 1; out.write(f"{s} {e} R {tw}\n")
                elif c["verdict"] == "never":
                    st["U_resolved_never"] += 1
                    if drop: drop.write(f"{s} {e} U {tw} {c.get('reason')}\n")
                else:
                    st["other_verdict"] += 1
                    if drop: drop.write(f"{s} {e} U {tw} OTHER:{c['verdict']}\n")
    st["seconds"] = round(time.time() - t0, 1)
    st["kept"] = st["P"] + st["U_resolved_positive"]
    json.dump(st, open(a.out + ".stats.json", "w"), indent=1)
    print(json.dumps(st))


if __name__ == "__main__":
    main()
