#!/bin/sh
# Unattended continuation once the three n=6 limit shards have written their tallies:
#  1. collect positivity survivors (records carrying "cert") into certificates/n6_survivors.idx
#  2. run the two-parameter certification on them in 3 shards
#  3. run the eps-slices on them (cross-checks) in a 4th process, sequentially
cd /home/user/CLTN
while [ ! -f certificates/region_n6_limit_shard0of3_tally.json ] || [ ! -f certificates/region_n6_limit_shard1of3_tally.json ] || [ ! -f certificates/region_n6_limit_shard2of3_tally.json ]; do sleep 60; done
python3 - <<'PY'
import json
idx = []
for s in range(3):
    for line in open(f"certificates/region_n6_limit_shard{s}of3.jsonl"):
        r = json.loads(line)
        if "cert" in r:
            idx.append(r["idx"])
idx.sort()
open("certificates/n6_survivors.idx", "w").write("\n".join(map(str, idx)) + "\n")
print("survivors:", len(idx))
PY
for s in 0 1 2; do
  nohup python3 run_region2d_sweep.py --n 6 --from certificates/region_n6_limit_shard0of3.jsonl certificates/region_n6_limit_shard1of3.jsonl certificates/region_n6_limit_shard2of3.jsonl --shard $s/3 > /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/region2d_n6_shard$s.log 2>&1 &
done
for e in 1/100 1/10 1/4 1/2 1/1000; do
  python3 run_region_sweep.py --n 6 --mode slice --eps $e --only certificates/n6_survivors.idx --crosscheck-every 5 > /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/slice_n6_$(echo $e | tr / -).log 2>&1
done
echo ALL_DONE > /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/after_limit.done
