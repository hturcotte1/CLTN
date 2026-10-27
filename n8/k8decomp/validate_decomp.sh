#!/bin/bash
# validate_decomp.sh n workdir -- (V1) the decomposition by bidirected subgraph produces exactly the same set
# of isomorphism classes as the reference enumeration directg -e(C+1):(2C) on K_n (canonical forms by
# labelg, multiset comparison, duplicate check), and (V2) filter8 applied to both streams gives the same
# set of (canonical digraph, arcs, flag) lines.  Exit 0 iff both hold.
set -uo pipefail
n=$1; W=$2/val$n; mkdir -p "$W"; C=$((n*(n-1)/2)); D=$(cd "$(dirname "$0")/../.." && pwd); K=$D/n8/k8decomp
python3 -c "import sys; sys.path.insert(0,'$K'); from g6tools import complete_g6; print(complete_g6($n))" > "$W/Kn.g6"
python3 "$K/gen_batches.py" $n 1:$C 1 "$W/batches" > "$W/gen.log"
s=$(date +%s); "$D/n7/directg" -q -e$((C+1)):$((2*C)) "$W/Kn.g6" > "$W/ref.d6"; t1=$(( $(date +%s) - s ))
s=$(date +%s); bash "$K/gen_batch.sh" "$W/batches/batch0.txt" 2> "$W/dec.genlog" > "$W/dec.d6" || { echo "gen_batch failed"; exit 1; }; t2=$(( $(date +%s) - s ))
"$D/n7/labelg" -q "$W/ref.d6" | sort > "$W/ref.can"; "$D/n7/labelg" -q "$W/dec.d6" | sort > "$W/dec.can"
nr=$(wc -l < "$W/ref.can"); nd=$(wc -l < "$W/dec.can"); nu=$(sort -u "$W/dec.can" | wc -l)
ok=0
if cmp -s "$W/ref.can" "$W/dec.can" && [ "$nd" = "$nu" ]; then echo "V1 n=$n: EQUAL sets, $nr classes, no duplicates (reference ${t1}s, decomposition ${t2}s)"; else echo "V1 n=$n: DIFFER ref=$nr dec=$nd unique=$nu"; ok=1; fi
"$D/n7/filter8" < "$W/ref.d6" 2> "$W/ref.stats" > "$W/ref.out"; "$D/n7/filter8" < "$W/dec.d6" 2> "$W/dec.stats" > "$W/dec.out"
paste -d' ' <(cut -d' ' -f1 "$W/ref.out" | "$D/n7/labelg" -q) <(cut -d' ' -f2,3 "$W/ref.out") | sort > "$W/ref.fcan"
paste -d' ' <(cut -d' ' -f1 "$W/dec.out" | "$D/n7/labelg" -q) <(cut -d' ' -f2,3 "$W/dec.out") | sort > "$W/dec.fcan"
if cmp -s "$W/ref.fcan" "$W/dec.fcan"; then echo "V2 n=$n: EQUAL filter outputs, $(wc -l < "$W/ref.fcan") lines; stats ref: $(cat "$W/ref.stats") | dec: $(cat "$W/dec.stats")"; else echo "V2 n=$n: DIFFER"; diff "$W/ref.fcan" "$W/dec.fcan" | head -n 5; ok=1; fi
exit $ok
