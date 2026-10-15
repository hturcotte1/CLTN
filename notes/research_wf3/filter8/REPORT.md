# filter8: overflow-safe n=8 version of n7/filter7.c

Scratch dir: /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/filter8/
Repository untouched.

## (a) Changes (filter8.c)
* Every long long op in the integer Faddeev-LeVerrier charpoly, every __int128 op in the
  (mu - 1/2) shift (binomials, 2^sh scaling, accumulation) and every product/difference in
  the Bareiss det128 (Hurwitz minors) is guarded with __builtin_{mul,add,sub}_overflow.
  Bareiss divisions are exact and |df/prev| <= |df|, so they cannot overflow.
* On any overflow the graph is KEPT and printed with flag 'O' (undecided); nothing aborts.
  stderr stats gain "overflow N (charpoly a shift b hurwitz c)".
* decode_d6 additionally checks the line is long enough for n*n bits.
* -DI128_IS_INT64 test build forces i128 to 64 bits so the guards must fire.
* Output format and every verdict are unchanged when no overflow occurs.

## (b) Validation (both binaries built with gcc -O2)
| input | lines | filter7 | filter8 | diff |
|---|---|---|---|---|
| all n=6 digraphs (geng 6 \| directg) | 1,540,944 | 3,594 surv (3,566 B, 28 S) | same | byte-identical, 0 overflow |
| n=7 geng 0/400 head 2,000,000 (task slice; sparse, 0 survivors) | 2,000,000 | 0 surv | 0 surv | identical, 0 overflow |
| n=7 geng 14:21 3/400 head 2,000,000 (dense) | 2,000,000 | 5 B | 5 B | byte-identical, 0 overflow |
| n=7 full published survivor list (n7_survivors_all.txt.gz, 1,889 S + 636,153 B) | 638,042 | 638,042 | 638,042 | byte-identical to each other AND to the published file, 0 overflow |
Overflow path: filter8_ovftest (64-bit i128) on the 3,000,007-graph n=8 sample: 10,387 'O',
1,044 'B', 0 'S'; every non-'O' line equals the real build's line, and every graph the real
build kept is either kept with the same flag or 'O' (no graph is lost, none changes verdict).

## (c) n=8 measurements (single filter process; machine load ~6 on 4 cores from the n7 sweep,
so all rates are pessimistic)
* geng -q 8 0/20000 | directg -q | filter8: 909 undirected graphs -> 107,955,318 digraphs,
  252 s wall for the pipeline (directg alone 119 s). 14,407,464 domination-free (13.3%),
  0 survivors, 0 overflow. (This slice is the sparse end of geng's order.)
* filter8 alone on pre-generated files:
  - 30M digraphs, geng 18:28 7/300 head (sparse orientations): 15.8 s = 1.90M graphs/s, 0 surv.
  - 30M digraphs, same with directg -e28:56: 41.2 s = 0.73M/s; 4,653 B, 0 S, 0 O.
  - 13.9M digraphs, geng 23:28 3/100, directg -e38:56: 11.5 s = 1.21M/s; 318 B, 0 S, 0 O.
  - 3,000,007 random dense digraphs (gen_dense.py, + K8 controls): 2.83 s = 1.06M/s;
    292,156 domfree, 3,280 B, 2,034 S (1,965 are K8 duplicates; 70 distinct S), 0 O.
  Fractions on the ~155M n=8 digraphs seen: S 2034/155M (sample-biased), B ~5e-5, O = 0 exactly.
* Exact Python cross-check (crosscheck8.py): independent big-int F2 (fastregion.icharpoly +
  Fraction shift + fastregion.idet Hurwitz) on 50 random distinct 'S' lines, all 'O' lines
  (there are none) and 20 'B' lines: 0 mismatches. fastregion.region(A, range(8), INF) on the
  same 70: K8 -> stable_somewhere; 68 -> never; and ONE non-clique graph
  &GLjlNZ}NzH\o (41 arcs, canonical &GKYM\R~Ryn^w) -> stable_somewhere, limit interval
  beta in (2, root of a degree-10 minpoly in (16/7, 23/10)), q in (1, ~1.29).
  ctln.analyze_support confirms stable at exact legal points eps=1/1000 and 1/100, beta=21/10.
  STATUS: CANDIDATE pending human hand-verification (CANDIDATE_n8_limit_GLjlNZ.json); not an
  established counterexample; found in a random (non-canonical) sample, not an enumeration.
  Max |g| coefficient seen at n=8: 5376 (~2^13), so the Hurwitz minors are ~2^13*7 bits at
  most in practice -- far below 2^127; the 0% overflow rate is expected to hold everywhere.

## (d) Full n=8 cost estimate (1,793,359,192,848 digraphs)
* filter8 only, at the measured 0.73-1.9M/s (population is mostly sparse, so ~1.3M/s weighted):
  ~1.4e6 CPU-s = ~380 CPU-h (range 260-680 CPU-h).
* Whole pipeline (geng + directg + filter8): directg generates ~0.9M/s, so ~0.6M digraphs per
  CPU-s overall -> ~2.9e6 CPU-s = ~800 CPU-h (optimistic ~500 CPU-h on an unloaded machine).
  On this 4-core box that is 5-9 days wall once the n=7 sweep finishes; rented compute
  (e.g. 64 cores) brings it to ~half a day. Output volume: B lines extrapolate to ~1e8
  (n=6: 3,566; n=7: 636,153), i.e. ~2 GB; the downstream limit stage must be run on S lines
  first and B lines (alpha = -1/2 exactly) handled as at n=7.
* Expected 'S' lines: n=6 28, n=7 1,889 (x67 while the population grew x572); extrapolating
  the S fraction (1.8e-5 -> 2.1e-6 -> ~2.5e-7) gives ~4.5e5, extrapolating the count ratio
  gives ~1.3e5: expect 1e5-5e5 'S' lines, each then needing fastregion.region(INF)
  (~0.01-0.1 s each, so negligible next to the filter).
