# n = 7, finite eps: the complete two-parameter map — runbook

Goal: for every digraph on 7 nodes decide exactly whether the full support carries a
stable positive fixed point at SOME legal (eps, delta), and for the ones that do, the
exact region in the (beta, t) plane (beta = 1 + delta/eps, t = (1-eps)/eps, legal iff
t(beta-2) > 1). Theorem F reduces the whole conjecture at n = 7 to this full-support
question plus the already certified maps for n <= 6.

Every stage is exact (integers, Fractions, exact root isolation; no floating point in
any decision path), every stage is a *necessary*-condition prune except the last, and
every stage was validated to diff zero against the n = 6 answer before use.

## Measured costs (this machine, one core)

| stage | tool | per graph | input | output | validation |
|---|---|---|---|---|---|
| 0. enumerate | `geng 7 r/K \| directg` (nauty 2.8.8) | ~5 us | — | 882,033,440 digraph6 lines = A000273(7) | count reproduced twice |
| 1. prefilter | `n7/filter7p` (C, exact) | ~3 us | all | 9,632,133 kept: 8,845,537 P + 786,596 U; domination-free 117,232,922 | all n=6: kept set == strict survivor set, diff 0 |
| 2. resolve U | `n7/resolve_U.py` | ~19 ms per U line | U lines | ~9.8% of U lines positive | n=6: 38/427 positive, set diff 0 |
| 3. cheap prune | `n7/cheap_prune.py` | 64 ms | all survivors (~8.9 M) | ~12.4% survive (~1.1 M) | n=6: 2,538/28,350 survive, all 4 stable graphs survive; 124 n=7 limit-stable graphs survive; 72,239 random exact points on killed graphs infeasible |
| 4. CAD | `n7/run_n7_region2d.py` (region2d) | 10.1 s | cheap survivors | verdict + exact cells per graph; HITS file | n=5, n=6 complete maps; 321-case projection check |
| 5. aggregate | (to write: `aggregate_n7.py`) | — | shard outputs | summary, equivalence classes, report | — |

Totals: stage 1 ~0.8 CPU-h (49 min wall on 4 cores); stage 2 ~4 CPU-h; stage 3 ~160 CPU-h;
stage 4 ~3,100 CPU-h (± the sample's ±0.3% on the 12.4%). About 3,300 CPU-hours in all.
On 4 cores: ~35 days. On one 96-vCPU machine: ~35 h. On four: ~9 h.
Memory: ~250 MB per Python process. Disk: stage-1 output 193 MB raw.

## Commands

```sh
sh n7/setup_remote.sh                              # sympy, nauty, filter build, A000273(6) + Geneson checks, test_ctln.py
K=96                                               # shards = cores
mkdir -p n7/pos n7/cheap n7/region2d
# stage 0+1 (all shards in parallel; each shard is a balanced res/mod slice of geng)
for k in $(seq 0 $((K-1))); do
  (n7/geng -q 7 $k/$K | n7/directg -q | n7/filter7p > n7/pos/part$k.txt 2> n7/pos/part$k.stats) &
done; wait
cat n7/pos/part*.stats | awk '{t+=$2; d+=$4; n+=$6; p+=$8; u+=$10} END {print "total",t,"domfree",d,"dropN",n,"P",p,"U",u}'
#   must print total 882033440 domfree 117232922 dropN 107600789 P 8845537 U 786596
# stage 2
for k in $(seq 0 $((K-1))); do
  (python3 n7/resolve_U.py n7/pos/part$k.txt --out n7/pos/part${k}_resolved.txt --dropped n7/pos/part${k}_droppedU.txt) &
done; wait
cat n7/pos/part*_resolved.txt > n7/n7_pos_survivors.txt
# stage 3 (resumable: rerun the same command after an interruption)
for k in $(seq 0 $((K-1))); do
  (python3 n7/cheap_prune.py n7/n7_pos_survivors.txt n7/cheap/shard${k}of$K.out --shard $k/$K --resume) &
done; wait
cat n7/cheap/shard*.out | awk '$2=="survives"{print $1" 0 C 0"}' > n7/n7_cheap_survivors.txt
# stage 4 (resumable; every non-'never' non-complete graph is appended to *_HITS.jsonl immediately)
for k in $(seq 0 $((K-1))); do
  (python3 n7/run_n7_region2d.py --from n7/n7_cheap_survivors.txt --shard $k/$K --out n7/region2d --resume) &
done; wait
cat n7/region2d/*_tally.json | python3 -c "import sys,json; ..."   # or aggregate_n7.py
```

Stop rules (CLAUDE.md §3): any line in a `*_HITS.jsonl` that is not already one of the 123
limit counterexamples (certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl) is a NEW
finite-eps counterexample invisible in the limit; the run stops there and the certificate
goes to a human. K7 must come out `stable_somewhere` in stage 4 (the driver exits 3 otherwise).

## What the outcome means

* 0 new hits: the finite-eps map at n = 7 equals the limit map; the 123 (+K7) limit-stable
  graphs are the complete list of non-clique stable full supports at n = 7, and the exact
  stable region of each is its CAD certificate. Combined with Theorem F and the n <= 6 maps
  this closes the conjecture at n = 7 for every (eps, delta).
* new hits: graphs stable only at finite eps (conditional stability); each one needs the
  pointwise exact witnesses and hand verification before anything is built on it.
