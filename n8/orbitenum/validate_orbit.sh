#!/bin/bash
# validate_orbit.sh <graph6 G> <min_arcs> <workdir> -- (V1) orbitenum's output is exactly the set of
# isomorphism classes that nauty's directg -e<min_arcs>:<2m> lists for G (canonical forms by labelg, multiset
# comparison, duplicate check) and (V2) filter8 gives the same (canonical digraph, arcs, flag) lines on both
# streams.  Exit 0 iff both hold.
set -uo pipefail
G=$1; lo=$2; W=$3; mkdir -p "$W"; D=$(cd "$(dirname "$0")/../.." && pwd)
m=$(python3 -c "import sys; sys.path.insert(0,'$D/n8/k8decomp'); from g6tools import edges_g6; print(len(edges_g6(sys.argv[1])))" "$G")
s=$(date +%s); printf '%s\n' "$G" | "$D/n7/directg" -q -e$lo:$((2*m)) > "$W/ref.d6"; t1=$(( $(date +%s) - s ))
s=$(date +%s); "$D/n8/orbitenum/orbitenum" "$G" "$lo" > "$W/orb.d6" 2> "$W/orb.log" || { echo "orbitenum failed"; exit 1; }; t2=$(( $(date +%s) - s ))
"$D/n7/labelg" -q "$W/ref.d6" | sort > "$W/ref.can"; "$D/n7/labelg" -q "$W/orb.d6" | sort > "$W/orb.can"
nr=$(wc -l < "$W/ref.can"); no=$(wc -l < "$W/orb.can"); nu=$(sort -u "$W/orb.can" | wc -l); ok=0
if cmp -s "$W/ref.can" "$W/orb.can" && [ "$no" = "$nu" ]; then echo "V1 $G lo=$lo: EQUAL sets, $nr classes, no duplicates (directg ${t1}s, orbitenum ${t2}s)"; else echo "V1 $G lo=$lo: DIFFER ref=$nr orb=$no unique=$nu"; ok=1; fi
"$D/n7/filter8" < "$W/ref.d6" 2> "$W/ref.stats" > "$W/ref.out"; "$D/n7/filter8" < "$W/orb.d6" 2> "$W/orb.stats" > "$W/orb.out"
paste -d' ' <(cut -d' ' -f1 "$W/ref.out" | "$D/n7/labelg" -q) <(cut -d' ' -f2,3 "$W/ref.out") | sort > "$W/ref.fcan"
paste -d' ' <(cut -d' ' -f1 "$W/orb.out" | "$D/n7/labelg" -q) <(cut -d' ' -f2,3 "$W/orb.out") | sort > "$W/orb.fcan"
if cmp -s "$W/ref.fcan" "$W/orb.fcan"; then echo "V2 $G lo=$lo: EQUAL filter outputs, $(wc -l < "$W/ref.fcan") lines; $(cat "$W/orb.stats")"; else echo "V2 $G lo=$lo: DIFFER"; ok=1; fi
exit $ok
