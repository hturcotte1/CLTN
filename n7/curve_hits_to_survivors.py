"""curve_hits_to_survivors.py -- merge the curve-stage shard outputs into the 4-column survivor format
for n7/run_n7_region2d.py, and print the stage tally.
   python3 n7/curve_hits_to_survivors.py n7/curve n7/n7_curve_survivors.txt"""
import sys, os, json, glob
src, out = sys.argv[1], sys.argv[2]
tot = {"lines": 0, "edges_lt_22": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "errors": 0, "shards": 0}
for p in sorted(glob.glob(os.path.join(src, "curve_stats_*.json"))):
    s = json.load(open(p)); tot["shards"] += 1
    for k in ("lines", "edges_lt_22", "never", "stable_somewhere", "escalate", "errors"):
        tot[k] += s.get(k, 0)
n = 0
with open(out, "w") as fh:
    for p in sorted(glob.glob(os.path.join(src, "curve_hits_*.jsonl"))):
        for line in open(p):
            r = json.loads(line)
            fh.write(f"{r['d6']} {r['edges']} {r['verdict'][0].upper()} 0\n"); n += 1
tot["written"] = n
json.dump(tot, open(os.path.join(src, "curve_TOTAL.json"), "w"), indent=1)
print(json.dumps(tot))
