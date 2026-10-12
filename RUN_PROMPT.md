# Run prompt

Open a Claude Code session in this folder with **Fable 5.1** as the
orchestrator. The four subagents are already defined in `.claude/agents/`
(`prover` → Opus 5.5, `symbolic-worker` → Sonnet 5.5, `sweeper` → Haiku 5.5,
`literature` → Sonnet 5.5). Paste everything below the line.

---

You are running a long autonomous research project. Read `CLAUDE.md` in full
before doing anything else — it has the problem statement, the fixed
conventions, the hard invariants, the measured scale, and the staged plan.
Then run `python3 test_ctln.py`. All 12 tests must pass before you touch
anything. If any fails, stop and report.

**The problem.** A directed graph determines a nonlinear dynamical system (a
CTLN). The target-free clique conjecture said its stable patterns are exactly
the target-free cliques. Geneson disproved it in July 2026 with a six-node
counterexample sitting just above the legal parameter minimum, while proving it
holds once the parameter `q` exceeds roughly `n−2`. **Nobody knows what happens
in the gap between.** Your job is to map that gap for `n = 6`, then `n = 7`.

**What already works.** `ctln.py` is an exact-arithmetic pointwise verifier,
validated against twelve independently known facts. `symbolic.py` computes exact
parameter regions with `q` kept symbolic — this is where the real answer lives,
because critical values are frequently irrational and no finite set of rational
samples will find a thin stable interval. `graphs.py` enumerates digraphs up to
isomorphism and reproduces OEIS A000273 exactly through `n=5` (9,608, matching
Curto's published package). `run_stage_a.py` is the cheap filter driver; a
smoke test at `n=4` returns exactly one stable full support, the complete
bidirected graph, precisely as the conjecture predicts.

**Your first objective is the validation gate.** Two things whose correct
answer is already known:

*Stage 0* — send `literature` for arXiv:2607.21396. Extract Geneson's explicit
six-node graph and its parameters, note which convention the paper uses and
whether you had to transpose, encode it, and confirm `ctln.py` flags it. Add it
as a permanent test. This is the sharpest test the pipeline will ever face and
it is not yet encoded.

*Stage 1* — reproduce the published `n=5` classification from
`ccurto/n5-graphs-package` with a diff of zero. Add it as a permanent
regression test.

**If both reproduce clean, keep going** — start Stage 2, the full-support sweep
at `n=6`. Shard `run_stage_a.py` across parallel `sweeper` agents, then escalate
the undecided set to `symbolic-worker` agents for Stage 3. Use `prover` for the
reduction lemmas, anything symbolic that resists, and the Lean formalization.
Work for as long as you have useful work.

**Two things stop the run cold.** First, any claimed counterexample — a
non-target-free support with a nonempty stable q-interval inside the legal
region. Write the full certificate, stop, and say so. Do not try to confirm it
yourself and do not build anything on top of it. Second, a failure to reproduce
Stage 0 or Stage 1. That is the finding; report it, do not work around it.

**Invariants, non-negotiable.** No floating point in any decision path —
exactness is the entire point, because the question is always which side of a
boundary a case sits on. Stability via Hurwitz determinants on the exact
characteristic polynomial, never a numerical eigenvalue routine. A vanishing
Hurwitz minor means *undecided*, never unstable and never a counterexample
(this bug was already found and fixed once; do not reintroduce it). Irrational
critical values get serialized exactly, never as decimals. No claim without a
certificate file that can be re-run. `test_ctln.py` passes before every commit.
Never use `graphs.dominates()` in a production sweep until `prover` has replaced
it with the paper's real definition and proved the preservation lemma — an
unsound prune hides counterexamples silently.

**Be ambitious about method.** Reformulate the stability condition if a better
basis exists. Train a classifier over graph features to steer the sweep toward
likely counterexamples — as a search heuristic only, with every hit still
requiring an exact certificate. Borrow from matrix analysis, numerical range,
oriented matroid theory, anything that helps. Attack a weaker statement first if
it opens the door. Look at `n=7` structure early if it hints at `n=6`. Hunt for
the counterexamples' structural signature rather than finding them one at a
time. If this problem stalls entirely, the other thirteen open questions in
Curto–Morrison §7.1 (arXiv:2301.12638) are adjacent and several are untouched.
Nothing is off the table on *how* you attack this. The invariants above are not
part of that latitude — a creative search that produces an uncertified result
has produced nothing.

**Keep progress countable.** "Processed 412,000 of 1,540,944; 11,803 undecided;
0 disagreements" is a status. "Made good progress" is not. Append every surprise
to `ANOMALIES.md` rather than resolving it quietly — a slow case, a critical
value of odd algebraic degree, a count that didn't match. The result usually
lives in the anomalies.

**Don't put a model in a deterministic loop.** Verification over a million cases
is a `for` loop, not a million inference calls. Reserve inference for the
lemmas, the resistant cases, the anomaly interpretation, the conjecture, and the
Lean work.

Have `literature` run a scoop check early and periodically — Geneson published
in July and mapping this gap is his obvious next move. If he's already on it,
that changes what we do.

Start now. Report at the validation gate, then continue.
