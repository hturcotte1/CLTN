#!/bin/sh
# run_n8_local.sh -- the n = 8 limit prefilter on this machine, in resumable shards at low priority.
#   sh n8/run_n8_local.sh <workers> <K>        e.g. sh n8/run_n8_local.sh 4 400
# Shard k of K: geng 8 15:28 k/K (underlying graphs with >= 15 edges; fewer cannot give >= 29 arcs)
#   | directg -e29:56 (digraphs with >= 29 arcs: the e_2 edge bound, proved and Lean-checked, says a stable
#     full support on 8 nodes needs at least 8*7/2 + 1 = 29 arcs; domination-free is tested by the filter)
#   | filter8 (overflow-guarded exact spectral filter; keeps S = alpha < -1/2 strictly, B = boundary, O = overflow)
# Output per shard: n8/parts/partK.txt (S and O lines only; B lines are counted in the .stats file) and
# n8/parts/partK.stats ("total ... domfree ... survivors ... boundary ... overflow ..."); n8/parts/partK.done marks completion.
cd "$(dirname "$0")/.." || exit 1
W=${1:-4}; K=${2:-400}
run_shard() {
  k=$1
  [ -f "n8/parts/part$k.done" ] && return
  nice -n 19 sh -c "n7/geng -q 8 15:28 $k/$K | n7/directg -q -e29:56 | n7/filter8 2> n8/parts/part$k.stats | awk '\$3!=\"B\"' > n8/parts/part$k.txt.tmp" && mv "n8/parts/part$k.txt.tmp" "n8/parts/part$k.txt" && touch "n8/parts/part$k.done"
}
k=0
while [ $k -lt $K ]; do
  while [ "$(jobs -p | wc -l)" -ge "$W" ]; do sleep 5; done
  run_shard $k &
  k=$((k+1))
done
wait
echo ALL_SHARDS_DONE
