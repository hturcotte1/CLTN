#!/bin/sh
D=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_prefilter_algebra
cd $D
for k in 00 01 02 03; do
  ( ./harness < n7U_part$k | python3 indep_pos.py > n7U_indep_$k.txt 2> n7U_indep_$k.err ) &
done
wait
echo done > indep.done
