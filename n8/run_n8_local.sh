#!/bin/bash
# run_n8_local.sh -- the n = 8 limit prefilter on this machine, in resumable shards at low priority.
#   bash n8/run_n8_local.sh <workers> <K>        e.g. bash n8/run_n8_local.sh 3 400
# Shard k of K: geng 8 15:28 k/K (underlying graphs with >= 15 edges; fewer cannot give >= 29 arcs)
#   | directg -e29:56 (digraphs with >= 29 arcs: the e_2 edge bound, proved and Lean-checked, says a stable
#     full support on 8 nodes needs at least 8*7/2 + 1 = 29 arcs; domination-free is tested by the filter)
#   | filter8 (overflow-guarded exact spectral filter; keeps S = alpha < -1/2 strictly, B = boundary, O = overflow)
# Output per shard: n8/parts/partK.txt (S and O lines only; B lines are counted in the .stats file) and
# n8/parts/partK.stats; n8/parts/partK.done marks completion. Re-running skips finished shards.
# NOTE: geng's res/mod classes are very uneven; most shards are tiny and a few dozen hold almost everything.
cd "$(dirname "$0")/.." || exit 1
W=${1:-3}; K=${2:-400}
run_shard() {
  local k=$1
  nice -n 19 sh -c "n7/geng -q 8 15:28 $k/$K | n7/directg -q -e29:56 | n7/filter8 2> n8/parts/part$k.stats | awk '\$3!=\"B\"' > n8/parts/part$k.txt.tmp" \
    && mv "n8/parts/part$k.txt.tmp" "n8/parts/part$k.txt" && touch "n8/parts/part$k.done"
}
running=0
for k in $(seq 0 $((K-1))); do
  [ -f "n8/parts/part$k.done" ] && continue
  while [ "$running" -ge "$W" ]; do wait -n; running=$((running-1)); done
  run_shard "$k" &
  running=$((running+1))
  echo "$(date -u +%H:%M:%S) started shard $k" >> n8/parts/progress.log
done
wait
echo "$(date -u +%H:%M:%S) ALL_SHARDS_DONE" >> n8/parts/progress.log
