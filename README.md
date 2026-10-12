# CTLN target-free clique project

Mapping the parameter gap where the target-free clique conjecture fails.

## Start here
1. `python3 test_ctln.py`  -> must print "12 passed, 0 failed"
2. Read `CLAUDE.md` (problem, conventions, invariants, staged plan)
3. Open Claude Code in this folder with Fable 5.1, paste `RUN_PROMPT.md`

## Layout
    ctln.py           exact pointwise verifier (no floats anywhere)
    symbolic.py       exact parameter-region analysis -- the centre of the project
    graphs.py         canonical forms + enumeration, validated vs OEIS A000273
    test_ctln.py      12 tests against independently known facts
    run_stage_a.py    Stage A cheap-filter driver, shardable
    CLAUDE.md         the operating document
    RUN_PROMPT.md     what to paste into Claude Code
    ANOMALIES.md      append-only surprise log
    .claude/agents/   prover (opus) / symbolic-worker (sonnet) /
                      sweeper (haiku) / literature (sonnet)
    certificates/     output

## Requires
    pip install sympy
Optional but much faster for n=6 enumeration: nauty (`directg`, `geng`).
If you swap in nauty, re-run the A000273 check first.
