# n = 8, eps -> 0 limit map -- runbook (needs rented compute)

Same pipeline as the n = 7 limit stage (notes/REPORT_n6.md, n7/run_n7_limit.py), with the
overflow-guarded filter `n7/filter8.c` (byte-identical to `filter7` on all n <= 7 data;
overflow -> flag 'O', kept). Measured on this machine: ~0.9 M digraphs/s per core, directg-bound;
0 overflows on 155 M eight-node digraphs.

| quantity | value |
|---|---|
| digraphs on 8 nodes (A000273(8)) | 1,793,359,192,848 |
| filter8 + directg rate | ~0.9 M / s / core |
| full pass | 500-800 CPU-hours (e.g. 96 vCPUs: 6-9 h; 4 cores here: ~6 days) |
| expected strict survivors 'S' | 1e5 - 5e5 (extrapolated), analysed by fastregion.region(INF) at ~50 ms each |
| boundary 'B' lines | ~1e8 (~2 GB); never limit-stable (exact), can be counted instead of stored |

```sh
sh n7/setup_remote.sh; gcc -O2 -o n7/filter8 n7/filter8.c
K=96; mkdir -p n8
for k in $(seq 0 $((K-1))); do
  (n7/geng -q 8 $k/$K | n7/directg -q | n7/filter8 2> n8/part$k.stats | awk '$3!="B"' > n8/part$k.txt) &
done; wait
cat n8/part*.stats     # totals must sum to 1793359192848 digraphs
# then: python3 n7/run_n7_limit.py adapted to n = 8 (support range(8)) over the 'S' and 'O' lines
```
The 'S' lines then get the exact limit region (fastregion), and every non-complete limit-stable
graph is a candidate for hand-verification, as at n = 7. Note geng's res/mod classes are uneven
at n = 7 (some empty, some ~2e8 digraphs); use K >= 10x the core count for load balance.
