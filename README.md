# CTLN target-free clique project

Mapping the parameter gap where the target-free clique conjecture fails.

## Start here
1. `python3 test_ctln.py`  -> must print "20 passed, 0 failed"
2. Read `CLAUDE.md` (problem, conventions, invariants, staged plan)
3. Open Claude Code in this folder with Fable 5.1, paste `RUN_PROMPT.md`

## Layout
    ctln.py           exact pointwise verifier (no floats anywhere)
    symbolic.py       exact parameter-region analysis -- the centre of the project
    graphs.py         canonical forms + enumeration, validated vs OEIS A000273
    test_ctln.py      20 tests against independently known facts (incl. Stage 0/1 reproductions)
    run_stage_a.py    Stage A cheap-filter driver, shardable
    CLAUDE.md         the operating document
    RUN_PROMPT.md     what to paste into Claude Code
    ANOMALIES.md      append-only surprise log
    .claude/agents/   prover (opus) / symbolic-worker (sonnet) /
                      sweeper (haiku) / literature (sonnet)
    certificates/     output (see certificates/INDEX.md)
    n7/               n = 7 pipeline: C prefilters, Lemma-R curve test, sharded drivers, runbooks
    notes/            lemmas, reports (REPORT_n6.md, REPORT_n7.md), audit and referee records

## Requires
    pip install sympy
Optional but much faster for n=6 enumeration: nauty (`directg`, `geng`).
If you swap in nauty, re-run the A000273 check first.

## Results so far (2026-10-15)
* n <= 6: complete exact map (notes/REPORT_n6.md). Exactly three non-clique graphs on 6 nodes have a
  stable full support; the conjecture holds at every legal parameter with q >= q_6* (a sextic root,
  about 1.1664) and fails for every q in (1, q_6*).
* n = 7: complete exact map at every legal parameter (notes/REPORT_n7.md): exactly 124 non-clique
  graphs, 123 visible in the eps -> 0 limit and one only at finite eps. All 124 are candidates pending
  human hand-verification.
* All n >= 6: counterexamples exist for every n (Theorem F), including an infinite family invisible in
  the limit (Theorem E); the critical q is unbounded in n (Theorem U); and at one fixed legal parameter
  the stable non-clique support size is unbounded (Theorem A). All theorems refereed by independent
  agents; 87 supporting lemmas machine-checked in Lean 4 + Mathlib (notes/lean/).
* Start with notes/FINAL_SUMMARY.md.
