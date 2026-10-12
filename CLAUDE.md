# CTLN target-free clique project

Map the parameter gap in which the target-free clique conjecture fails.

---

## 1. The problem

A **CTLN** turns a directed graph into a dynamical system. For a graph `G` on
`n` nodes and two positive rational knobs `eps`, `delta`:

```
dx_i/dt = -x_i + [ sum_j W[i][j] x_j + theta ]_+
```

**Curto–Morrison Open Question 10 (the target-free clique conjecture):** the
supports of stable fixed points are exactly the target-free cliques.

**Geneson, arXiv:2607.21396 (Jul 2026)** disproved it with a 6-node graph at
`q ≈ 29/25`, and proved it *does* hold for

```
q  >=  n - 2 - (n-3)*eps/2        where  q := delta*(1-eps)/eps
```

A CTLN is only legal for `q > 1`. **The open region is the gap between `q`
slightly above 1 and `q` around `n-2`.**

**Our goal:** settle that gap for `n = 6`, then `n = 7`. Three outcomes are all
real results — finding more counterexamples and characterizing what they share;
pushing the proven-safe threshold down from `n-2` toward 1; or a complete
certified map showing where the boundary actually is.

---

## 2. Conventions — DO NOT DRIFT

Long runs on this problem die the same way: a sign or transpose convention
flips a thousand cases in, and nothing afterward means anything. These are
fixed. If you think one is wrong, log it and stop — don't change it.

```
A[i][j] == 1   iff   there is an edge  i -> j      (row i = who i SENDS to)

W[i][j] = 0           if i == j
W[i][j] = -1 + eps    if A[j][i] == 1    (j sends to i -> WEAK inhibition)
W[i][j] = -1 - delta  otherwise          (no j->i edge -> STRONG inhibition)
```

Note the transpose: `W[i][j]` is the influence **of j on i**, so it reads
`A[j][i]`. This is the most common place to go wrong.

**Closed form, cross-checked in both `ctln.py` and `symbolic.py`:**

```
I - W_sigma  =  -delta*I  +  (1+delta)*J  -  (eps+delta)*A_sigma
```

Every entry is a polynomial in `(eps, delta)` with integer coefficients, which
is what makes exact and symbolic work possible.

**`sigma` supports a stable fixed point iff all three hold:**

1. `(I - W_sigma) x = theta*1` has a solution with `x > 0` strictly
2. for every `k` outside `sigma`: `theta + sum_{j in sigma} W[k][j] x_j < 0`
3. all eigenvalues of `-I + W_sigma` have real part `< 0`

**Graph notions.** `sigma` is a *clique* if every distinct pair has edges both
ways. A node `k` outside `sigma` is a *target* if every member sends to `k`.
`sigma` is *target-free* if it is a clique with no target.

---

## 3. Hard invariants

Not style preferences. Violating any one invalidates everything downstream.

- **No floating point in any decision path.** The interesting region sits near
  a boundary where a float solver cannot tell `-1e-12` from `+1e-12`. Use
  `fractions.Fraction` or symbolic parameters.
- **Stability decided by Hurwitz determinants on the exact characteristic
  polynomial.** Never a numerical eigenvalue routine.
- **Boundary cases are three-valued.** A vanishing Hurwitz minor means
  *undecided* — not unstable, not a counterexample. Collapsing undecided into
  either bucket makes the search chase phantoms. This bug was already found and
  fixed once; do not reintroduce it.
- **Critical values are often irrational.** An early test produced
  `sqrt(2261)/14 - 2`. Serialize exactly. A decimal in a certificate makes it
  worthless, because the entire question is which side of a boundary you are on.
- **No claim without a re-runnable certificate file.** No result exists until it
  round-trips.
- **`test_ctln.py` passes before any commit.** All 12, every time.
- **A claimed counterexample stops the run** and goes to a human for
  hand-verification before anything is built on it.

---

## 4. The full-support simplification — read this before designing a sweep

When `sigma` is the whole vertex set there are no outside nodes, so **condition
2 is vacuous** and there is nothing that can be a target. "Target-free clique"
then collapses to just "clique", which for the full vertex set means the
complete bidirected graph. So:

> **Restricted to full support, the conjecture says: the full vertex set is a
> stable fixed point if and only if the graph is complete bidirected.**

That is sharply testable over graphs alone — no support enumeration, no
off-support condition. And it is exactly the form of Geneson's counterexample
(a stable fixed point with non-clique full support), so this is simultaneously
the cleanest attack surface and the place the known counterexample lives.

**Start here.** The size-5-support case within `n=6` needs condition 2 back and
should come second.

---

## 5. Measured scale and the two-stage architecture

Numbers from this machine, single-threaded:

| step | rate | implication |
|---|---|---|
| enumeration, n≤5 | 3 s for all 9,608 | validated, see §6 |
| pointwise exact (`ctln.analyze_support`), n=6 | **221 cases/sec** | full 1,540,944 sweep at one q ≈ **2 hours** |
| undecided (BOUNDARY) rate at eps=1/2, q=3/2 | **2.8%** | ≈ 43,000 escalations |
| symbolic `q_region`, n=5 full support | **1.4 s/case** | heavier at n=6 — cannot run on everything |

So the architecture is forced, and it is the right one:

**Stage A — cheap exact filter.** Run `ctln.analyze_support` over the whole
space at several `q` values spread across the legal gap. Exact rational, no
floats, ~2 hours per q value, trivially parallel. Most cases are decided here.

**Stage B — symbolic escalation.** Run `symbolic.q_region` only on the
survivors: everything undecided, everything that disagrees, and everything whose
verdict changes between sampled q values. This is where the real answer comes
from, because a case stable on a thin interval between two algebraic endpoints
is invisible to any finite set of samples.

**Stage C — human and prover.** Cases that resist, the reduction lemmas, the
structural conjecture, the Lean formalization.

A sweep of 1,495 pairs at sampled q found zero decided counterexamples and 84
undecided cases. **The counterexamples are in the undecided pile.** Design
accordingly: Stage A is a filter, not an answer.

---

## 6. Prune hard, but prove every prune

Already proven, so skip: **symmetric graphs**, **oriented graphs**, **supports
of size ≤ 4**. At `n=6` the live cases are *mixed* graphs with supports of size
5 or 6.

Curto's domination theorems (arXiv:2510.05098) would let you delete dominated
nodes and search only irreducible graphs — a large further cut.
**`graphs.dominates()` is a placeholder written from a secondhand summary and
must not be used in a production sweep until the real definition is pulled from
the paper and the preservation lemma is proved.** An unsound prune silently
hides counterexamples and nobody finds out.

**Validation anchor, already passing:** the count of digraphs on `n` nodes up to
isomorphism is OEIS A000273 — `1, 1, 3, 16, 218, 9608, 1540944`. Our enumerator
reproduces `n ≤ 5` exactly, and the 9,608 matches the figure in Curto's
published n=5 package. If any replacement enumerator fails this check, stop.

---

## 7. Stages and definition of done

Each stage has a countable output. Do not start the next until the current one
is complete and checked.

**Stage 0 — reproduce Geneson's counterexample.** Pull arXiv:2607.21396, extract
the explicit 6-node graph and parameters, confirm the verifier flags it, add it
as a permanent test. *Done when:* it round-trips and matches the paper's claim.

**Stage 1 — reproduce the published n=5 classification.** Diff of zero against
`ccurto/n5-graphs-package`. *Done when:* every graph agrees. **If you cannot
recover a published result, nothing at n=6 is trustworthy.**

**Stage 2 — full-support sweep at n=6.** Stage A filter over all 1,540,944
canonical digraphs at several q values. *Done when:* every graph has a
certificate and the undecided set is enumerated.

**Stage 3 — symbolic regions on the escalation set.** Exact q-intervals for
every survivor. *Done when:* each has `always` / `never` / explicit region.

**Stage 4 — size-5 supports within n=6.** Condition 2 comes back. Same shape.

**Stage 5 — aggregate.** Is there a `q*` below which counterexamples appear?
What do they share structurally? *Done when:* the n=6 gap is mapped and written.

**Stage 6 (stretch) — n=7 and the general theorem.** Not the definition of
success. A complete certified n=6 map is a real result on its own.

---

## 8. Latitude

Be aggressive about method. Explicitly on the table: reformulating the stability
condition in a basis that makes the structure visible; learning a classifier
over graph features to steer the sweep toward likely counterexamples (as a
*search heuristic* only — every hit still needs an exact certificate); borrowing
machinery from matrix analysis, numerical range, or oriented matroid theory;
attacking a weaker statement first; looking at n=7 structure early for hints
about n=6; characterizing the counterexample's signature rather than hunting
one at a time; and if this problem stalls, the other 13 open questions in
Curto–Morrison §7.1 (arXiv:2301.12638) are nearby and some are untouched —
"proper sources die" (Q13) needs no data at all.

The latitude is about **method**, never about the invariants in §3 or about
what counts as verified. A creative search that produces an uncertified result
has produced nothing.

---

## 9. Anomaly log

Keep `ANOMALIES.md`, append-only. Anything surprising goes in: a case far slower
than its neighbours, a critical value of unexpected algebraic degree, a count
that didn't match, a reduction behaving oddly. Do not tidy these away. The
result usually lives in the anomalies, and a run that silently smooths over
surprises produces a clean-looking report about nothing.

---

## 10. Where inference is actually needed

Running the verifier over a million cases is deterministic code. It should be a
loop, not a million model calls. Reserve inference for: proving the reduction
lemmas, symbolic cases that resist, interpreting the anomaly log, forming the
structural conjecture, the Lean formalization, and deciding what to do when
something doesn't fit. Batch verification, file wrangling and format conversion
run as plain scripts with no model in the loop.

**Agent roster** (`.claude/agents/`):

| agent | model | job |
|---|---|---|
| `prover` | Opus 5.5 | lemmas, resistant cases, the conjecture, Lean |
| `symbolic-worker` | Sonnet 5.5 | Stage B symbolic batches |
| `sweeper` | Haiku 5.5 | Stage A mechanical sweeps, collation, counts |
| `literature` | Sonnet 5.5 | paper retrieval, verbatim definitions, scoop checks |

Orchestrator runs on Fable 5.1 and should stay out of the volume work.

---

## 10b. Network access — read before Stage 0

Outbound network from the **shell** is usually allowlisted: package registries
and GitHub work, most other sites do not. `curl`/`wget` to arxiv.org will often
be refused before it reaches the site, and every shell route goes through the
same allowlist, so retrying with a different command will not help.

- **Papers:** use `WebFetch` / `WebSearch`, not curl. `arxiv.org/abs/...` reads
  fine; try `arxiv.org/pdf/...` for full text.
- **Repos:** `git clone` or `gh` for public GitHub repos.
- **Python packages:** `pip install sympy` (add `--break-system-packages` if pip
  refuses).

If a fetch is refused, **say so and ask for the file to be attached** rather
than reconstructing a paper's contents from citations or memory. For Stage 0
specifically, a wrong guess at Geneson's graph is worse than no graph at all.

---

## 11. Files

- `ctln.py` — exact pointwise verifier. Read the module docstring first.
- `symbolic.py` — exact parameter-region analysis. The centre of the project.
- `graphs.py` — canonical forms and enumeration, validated against A000273.
- `test_ctln.py` — 12 tests against independently known facts. Run first.
- `ANOMALIES.md` — append-only surprise log.
- `certificates/` — one file per stage.

---

## 12. People to contact early

**Jesse Geneson** — published the counterexample in July; mapping the gap is his
natural next move, so find out now whether he is already on it. **Anne Shiu**,
**Carina Curto** — the surrounding program. Send the n=5 reproduction as
evidence the pipeline is real.
