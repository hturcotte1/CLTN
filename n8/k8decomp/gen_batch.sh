#!/bin/bash
# gen_batch.sh <batchfile> [tag] -- digraph6 stream of every digraph on K_n whose bidirected subgraph is one
# of the batch's graphs B: for each line "B H e" of the batch file,  printf H | directg -q -o | addarcs B.
# directg -o lists the orientations of H (one per orbit of Aut(H) = Aut(B)); addarcs ORs in the bidirected
# arcs of B.  Per-B line counts go to stderr (addarcs).  Any failure stops the stream with exit 1.
# The optional tag is ignored: the n = 8 driver passes its job<j>.g6 path so that its running-job guard
# (a grep of the command lines for job<j>.g6) sees this pipeline.
set -o pipefail
D=$(cd "$(dirname "$0")/../.." && pwd)
while read -r B H e; do
  printf '%s\n' "$H" | "$D/n7/directg" -q -o | "$D/n8/k8decomp/addarcs" "$B" || { echo "gen_batch: FAILED on B=$B" >&2; exit 1; }
done < "$1"
