# Lemma F: the target-free clique conjecture reduces to full supports

Author: prover agent, 2026-10-13. Status of each item is labelled
**[PROVED]** (complete proof below), **[CITED]** (taken from a paper, not
re-derived here), or **[SANITY]** (an exact computation that supports a claim
but does not prove it). Nothing here changes the conventions of CLAUDE.md §2.

Sanity scripts (exact arithmetic only, `python3 -I`), all in
`/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverF/`:
`check_FGH.py` (every lemma below, on random cases; prints `ALL SANITY CHECKS PASSED`),
`small_k_region2d.py` / `.json` (k = 3, 4 CAD sweep), `undecided_rate.py`,
`undecided_kinds.py`, `undecided_factor.py`, `undecided_probe.py`, `lemma21_search.py`.

---

## 0. Setting and notation

Graph `G` on `[n]`, `A[i][j] = 1` iff `i -> j`. For `i != j` write

    c_ij := -W_ij = 1 - eps   if j -> i        (A[j][i] = 1)
                  = 1 + delta if j -/-> i

so `(I - W)_ii = 1` and `(I - W)_ij = c_ij` for `i != j`. Parameters: `eps > 0`,
`delta > 0`, `theta > 0`. **Legal** means `q := delta(1-eps)/eps > 1`.

* Legal forces `eps < 1` (if `eps >= 1` then `q <= 0`).
* Identity used repeatedly: `(1+delta)(1-eps) - 1 = delta(1-eps) - eps = eps(q-1)`,
  so legal iff `(1+delta)(1-eps) > 1`.

For a nonempty support `sigma` write `W_sigma = (W_ij)_{i,j in sigma}`. The
row-`i` fixed-point equation and the drive at an outside node `k` are

    (FP_i)     x_i + sum_{j in sigma, j != i} c_ij x_j = theta          (i in sigma)
    drive_k(x) := theta + sum_{j in sigma} W_kj x_j = theta - sum_{j in sigma} c_kj x_j   (k not in sigma)

`sigma` **supports a stable fixed point** (written SFP) iff
(1) `(I - W_sigma) x = theta 1` has a solution `x > 0`;
(2) `drive_k(x) < 0` for every `k` not in `sigma`;
(3) every eigenvalue of `-I + W_sigma` has real part `< 0`.

**Observation 0 [PROVED].** If (3) holds then `det(I - W_sigma) > 0`. In
particular `I - W_sigma` is nonsingular, the `x` in (1) is unique, and (2) is
unambiguous. If `I - W_sigma` is singular, (3) fails and `sigma` is not an SFP.

*Proof.* `I - W_sigma = -(-I + W_sigma)`, so its eigenvalues have positive real
part. It is a real matrix, so its non-real eigenvalues come in conjugate pairs
with product `|mu|^2 > 0`, and its real eigenvalues are positive. The
determinant is the product of the eigenvalues. ∎

**Dynamical meaning (remark, [PROVED]).** If (1), (2), (3) hold, then
`x* = (x on sigma, 0 off sigma)` is an equilibrium of the CTLN: for `i in sigma`
the argument of `[.]_+` equals `x_i > 0` by (FP_i), and for `k` outside it is
`drive_k < 0`. All these arguments are bounded away from 0 near `x*`, so near
`x*` the vector field is affine with Jacobian `[[-I + W_sigma, *], [0, -I]]`.
Its spectrum is `spec(-I + W_sigma)` together with `-1`, so (3) gives local
asymptotic stability. This is the "core certificate lemma" proposed for Lean.

Notation for this note. For a graph `H` on `k` nodes and a legal point
`p = (eps, delta)`:

* `FS(H, p)`: the full support `[k]` of `H` satisfies (1) and (3) at `p`.
  Condition (2) is vacuous for the full support.
* `FS(k, p)`: every `k`-node graph `H` with `FS(H, p)` is complete bidirected.
* `Conj(n, p)`: for every graph `G` on `n` nodes and every nonempty `sigma`,
  `sigma` is an SFP of `G` at `p` iff `sigma` is a target-free clique.

By Lemma G these do not depend on `theta`.

---

## 1. Lemma G: `theta` only rescales [PROVED]

**Lemma G.** Fix `G`, `sigma`, `eps`, `delta`. For every `theta > 0`, each of
(1), (2) and (3) has the same truth value at `theta` as at `theta = 1`. The
same holds for every boundary event: singular `I - W_sigma`, some `x_i = 0`,
some `drive_k = 0`, or a vanishing characteristic-polynomial coefficient or
Hurwitz minor. Hence `theta = 1` loses no generality.

*Proof.* `W` does not involve `theta`. The map `x -> y = x/theta` is a bijection
from the solutions of `(I - W_sigma) x = theta 1` to the solutions of
`(I - W_sigma) y = 1`. Because `theta > 0`, `x_i > 0` iff `y_i > 0`, and
`x_i = 0` iff `y_i = 0`. Next, `drive_k(x; theta) = theta (1 + sum_j W_kj y_j)
= theta drive_k(y; 1)`, which has the same sign. Condition (3), the
characteristic polynomial of `-I + W_sigma` and its Hurwitz minors do not
involve `theta`, and neither does singularity of `I - W_sigma`. ∎

*Dynamical version.* Since `[theta u]_+ = theta [u]_+` for `theta > 0`, the map
`x(t) -> x(t)/theta` carries trajectories of the CTLN with input `theta` to
trajectories with input 1. It is linear, keeps the orthant and keeps supports,
so equilibria, their supports and their stability correspond exactly.

[SANITY] `check_FGH.py::check_G`: on 300 random cases, every field of
`ctln.analyze_support` except `x` and `theta` is identical for `theta = 1` and
a random rational `theta`. `x` scales exactly. So `fastregion.py`,
`region2d.py` and the sweeps lose nothing by fixing `theta = 1`.

---

## 2. Lemma R: conditions (1) and (3) only see the induced subgraph [PROVED]

**Lemma R (restriction).** Let `sigma = {s_1 < ... < s_k}` and let `H = G|sigma`
be the graph on `[k]` with `A^H[a][b] = A[s_a][s_b]`. Then for all `eps, delta`:

    W^H[a][b] = W[s_a][s_b]      for all a, b in [k],

that is, `W^H = W_sigma` under the relabelling `a <-> s_a`. Therefore:

* (1) for `sigma` in `G` iff (1) for the full support of `H`. It is the same
  linear system, with the same solutions.
* (3) for `sigma` in `G` iff (3) for the full support of `H`. It is the same
  matrix.
* (1) and (3) together for `sigma` in `G` iff `FS(H, p)`.
* Condition (2) for `sigma` in `G` reads `W_kj` for `k` outside `sigma` and `j`
  inside, which depends on the edges `j -> k` from `sigma` to the outside. Those
  are not edges of `H`, so (2) is **not** a property of `H`.

*Proof.* Both diagonals are 0. For `a != b`: `W^H[a][b] = -1 + eps` iff
`A^H[b][a] = 1` iff `A[s_b][s_a] = 1` iff `W[s_a][s_b] = -1 + eps`. Otherwise
both entries are `-1 - delta`. ∎

**Relabelling invariance [PROVED].** Let `pi` be a permutation of `[k]` and
`H'` the graph with `A'[pi(a)][pi(b)] = A^H[a][b]`. The same entrywise argument
gives `W'[pi(a)][pi(b)] = W^H[a][b]`, so `I - W' = P (I - W^H) P^T` with `P` the
permutation matrix. Solutions correspond via `x' = P x`, which keeps
positivity, and the characteristic polynomials are equal. Complete
bidirectedness is an isomorphism invariant. So `FS(H, p)` depends only on the
isomorphism class of `H`, and sweeping one canonical representative per class
(the A000273 enumeration) loses nothing.

[SANITY] `check_FGH.py::check_F`, `check_perm`: `ctln.submatrix(W, sigma) ==
ctln.weight_matrix(G|sigma)` on 400 random `(G, sigma, p)`. The `x`,
`positive`, `stable`, `charpoly` and `BOUNDARY` fields agree. Relabelling
invariance holds on 200 cases.

---

## 3. Lemma R2: the clique rule, proved directly [PROVED]

This is Curto–Morrison Rule 2 (cliques) and Curto–Geneson–Morrison 2019 [CITED].
Here is a self-contained proof in our notation.

**Lemma R2.** Let `0 < eps < 1`, `delta > 0`, `theta > 0`, and let `sigma` be a
clique of size `m`. Legality is not needed.

(a) `I - W_sigma = eps I + (1-eps) J`. Its unique solution is
    `x = c 1` with `c = theta / (eps + m(1-eps)) > 0`, so (1) holds.
(b) `-I + W_sigma = -(eps I + (1-eps) J)` has eigenvalues `-eps` (multiplicity
    `m-1`, on the complement of `1`) and `-(eps + m(1-eps))` (on `1`). So (3)
    holds.
(c) Let `k` be outside `sigma`, and let `r` be the number of `j in sigma` with
    `j -/-> k`. Then `drive_k = theta [eps - r(eps+delta)] / (eps + m(1-eps))`.
    This is `> 0` if `r = 0` (k is a target) and `<= -theta delta/(eps + m(1-eps)) < 0` if `r >= 1`.
(d) Hence a clique `sigma` is an SFP iff it is target-free.

*Proof.* (a) For `i != j` in `sigma` both edges exist, so `c_ij = 1 - eps`.
Then `(eps I + (1-eps)J) c1 = c(eps + m(1-eps)) 1 = theta 1`. The matrix is
nonsingular by (b), so the solution is unique.
(b) `J` has eigenvalue `m` on `1` and `0` on the `(m-1)`-dimensional complement.
(c) `drive_k = theta + c[(m-r)(-1+eps) + r(-1-delta)]`. Multiply by
`(eps + m(1-eps))/theta = 1/c`: you get
`eps + m - m eps - m + (m-r) eps - r delta = eps - r(eps + delta)`.
If `r >= 1` this is at most `eps - (eps + delta) = -delta`.
(d) follows from (a), (b) and (c). ∎

[SANITY] `check_FGH.py::check_R2_symbolic` checks (a), (b) and (c)
symbolically in `(eps, delta, theta)` for `m <= 7`. `check_R2_numeric` checks
(d) with `ctln` on 300 random cliques inside random graphs.

---

## 4. Lemma T: what an outside node sees, given only (1) [PROVED]

**Lemma T.** Let `x > 0` solve (FP) for `sigma` in `G`, and let `k` be outside
`sigma`. Write `S := {j in sigma : j -> k}` and
`In(i) := {j in sigma \ {i} : j -> i}`. Only `eps, delta > 0` and `eps < 1` are
used.

(a) **Monotone.** `drive_k` depends on `k` only through `S`. If `S' ⊆ S`, then
    `drive(S') = drive(S) - (eps + delta) sum_{j in S \ S'} x_j`, which is
    `<= drive(S)`, with strict inequality when `S' != S`.
(b) **Inside-out.** If some `i in sigma` has `i notin S` and `S ⊆ In(i)`, then
    `drive_k <= -delta x_i < 0`. In particular `S = ∅` always gives (2) at `k`,
    and so does `S = In(i)` for any `i in sigma`.
(c) **Outside-in.** If some `i in sigma` has `i in S` and `In(i) ⊆ S`, then
    `drive_k >= eps x_i > 0`, so (2) fails. In particular a target (`S = sigma`)
    always violates (2), for every `sigma`, clique or not.

*Proof.* (a) `c_kj` is `1 - eps` for `j in S` and `1 + delta` otherwise.
For (b) and (c), subtract (FP_i) from the definition of the drive:

    drive_k = (1 - c_ki) x_i + sum_{j in sigma, j != i} (c_ij - c_kj) x_j.        (*)

(b) `i notin S` gives `c_ki = 1 + delta`, so the first term is `-delta x_i`. A
summand with `c_ij - c_kj > 0` needs `c_ij = 1 + delta` and `c_kj = 1 - eps`,
that is `j -/-> i` and `j -> k`, so `j in S \ In(i) = ∅`. Hence every other
term is `<= 0`.
(c) `i in S` gives `c_ki = 1 - eps`, so the first term is `eps x_i`. A summand
with `c_ij - c_kj < 0` needs `j in In(i) \ S = ∅`. Hence every other term is
`>= 0`. ∎

(b) and (c) are the `sigma`-relative graphical-domination configurations
("`i` dominates `k`" and "`k` dominates `i`"). I believe they are the
inside-out and outside-in cases in Curto–Geneson–Morrison 2019, but I have not
checked that paper's wording, and the proof above does not depend on it. The
same subtraction applied to two nodes `j, k` inside `sigma` gives the
inside-inside case, used only in §8 (Lemma D'):

**Lemma D' [PROVED].** Let `j, k in sigma` with `j -> k`, `k -/-> j`, and
`l -> j => l -> k` for every `l in sigma \ {j,k}`. Then (FP) has no solution
`x > 0`. *Proof:* (FP_j) minus (FP_k) gives
`0 = eps x_j + delta x_k + sum_{l != j,k} (c_jl - c_kl) x_l`, and every summand
is `>= 0`. ∎
(This is the `sigma`-relative form of Curto arXiv:2510.05098 Def. 1.1 /
Cor. 3.6, as recorded in `graphs.py`. The full domination lemma belongs to
`notes/lemmas.md` and is not claimed here.)

[SANITY] `check_FGH.py::check_E`: on 1,287 random positive cases, (a), (b) and
`S = ∅` hold for every in-set `S`, and they agree with
`ctln.analyze_support(...)["off_support_ok"]` on the padded graph.
`check_H_strong` checks (b) and (c) on 398 cases.

---

## 5. Lemma K2: two nodes [PROVED]

**Lemma K2.** At every legal `(eps, delta)`, no non-complete 2-node graph
satisfies (3) on its full support. So `FS(2, p)` holds for every legal `p`.

*Proof.* `I - W = [[1, c_12], [c_21, 1]]` and
`det(lam I - (-I + W)) = (lam + 1)^2 - c_12 c_21`, so the eigenvalues are
`-1 ± sqrt(c_12 c_21)`, where `c_12 c_21 > 0`. A non-complete graph has some
`c = 1 + delta`, and the other `c` is `>= 1 - eps > 0`. Therefore
`c_12 c_21 >= (1+delta)(1-eps) = 1 + eps(q-1) > 1`, so `-1 + sqrt(c_12 c_21) > 0`. ∎

(For the one-way edge, `det(C + tJ) = -(t(beta-2) - 1)`. The legal boundary
`q = 1` is exactly where this 2-node determinant vanishes. That is why
`region2d.py` raises "common factor between D and L" on `k = 2`; see §8.)

---

## 6. Theorem F: the reduction [PROVED]

**Theorem F (pointwise).** Fix a legal `p = (eps, delta)` (and any `theta > 0`,
by Lemma G) and `n >= 1`. The following are equivalent:

(i)   `Conj(n, p)`.
(ii)  For every `n`-node graph `G`, every SFP of `G` at `p` is a clique.
(iii) `FS(k, p)` for every `k in {2, ..., n}`. That is, no non-complete graph on
      `k <= n` nodes has a full support satisfying (1) and (3) at `p`.

(`k = 1` is vacuous: the 1-node graph is complete.)

*Proof.*

**(i) <=> (ii).** The direction "target-free clique => SFP" holds at every `p`
by Lemma R2(d). The direction "SFP => target-free clique" splits into "SFP =>
clique", which is (ii), and "SFP clique => target-free", which is Lemma R2(d)
or Lemma T(c). So (i) and (ii) are equivalent at every `p`.

**(F1) (iii) => (ii).** Let `sigma` be an SFP of an `n`-node `G` at `p`. Case
split:

* *`sigma` is a clique.* Nothing to prove for (ii). For (i), Lemma R2(d) also
  shows it is target-free. **This branch is the only place condition (2) is
  used.** The SFP property supplies (2) at every outside node, and (2) at a
  target is impossible (R2(c) or T(c)).
* *`sigma` is not a clique.* Then `|sigma| = k` with `2 <= k <= n`, and
  `H := G|sigma` is not complete bidirected. By Lemma R, (1) and (3) for
  `sigma` in `G` give `FS(H, p)`, which contradicts `FS(k, p)`. **This branch
  uses only (1) and (3), never (2).** So this case cannot occur.

**(F2) (ii) => (iii), by padding.** Suppose `H` is a non-complete graph on
`k <= n` nodes with `FS(H, p)`. Build `G` on `n` nodes: `G|[k] = H`, and the
`n - k` new nodes receive **no** edges from `[k]`. Their other edges (among
themselves and from them into `[k]`) are arbitrary, since they do not affect
any condition for `sigma = [k]`. By Lemma R, `sigma = [k]` satisfies (1) and
(3) in `G`. Each new node has `S = ∅`, so Lemma T(b) gives (2). So `sigma` is
an SFP of `G` and is not a clique, and (ii) fails. ∎

**Correction to the task statement.** "Condition (2) is never needed in (F1)"
is **false as literally stated**. Take the graph `0 -> 1` and `sigma = {0}`.
This `sigma` is a clique and satisfies (1) and (3), but node 1 is a target and
`drive_1 = eps theta > 0`, so `{0}` is not an SFP. Without (2), (F1) would
wrongly count `{0}` as a stable support. What is true:

* the hypothesis `FS(k, p)` involves no condition (2);
* the non-clique branch never uses (2);
* (2) enters only to exclude targets of cliques, through Lemma R2 or Lemma T(c).

**Corollary F1 (uniform form, the task's (F1)).** Let `P` be any set of legal
parameters. If `FS(k, p)` holds for all `2 <= k <= n` and all `p in P`, then
`Conj(n, p)` holds for every `p in P`. With `P` = all legal parameters: if for
every `k <= n`, every `k`-node graph whose full support satisfies (1) and (3)
at some legal parameter is complete bidirected, then the target-free clique
conjecture holds for every `n`-node graph at every legal parameter. The
converse also holds pointwise, by (F2).

**Corollary F2 (where a full-support counterexample transplants).** Let `H` be
a non-complete `k`-node graph with `FS(H, p)`. Let `G` be any graph with
`G|sigma ≅ H`. Then `sigma` is a non-clique SFP of `G` at `p`, so `G` is a
counterexample at `p`, iff (2) holds at every outside node. Whether (2) holds
at `k` depends only on `S_k = {j in sigma : j -> k}`. The admissible in-sets
`S` form a down-set (T(a)) that:

* always contains `∅`;
* always contains every open in-neighbourhood `In(i)` for `i in sigma`, and
  every subset of one (T(b));
* never contains any superset of a closed in-neighbourhood `{i} ∪ In(i)`, in
  particular never `S = sigma` (T(c)).

Everything else is decided by the exact inequality of Lemma H(4) below. In
particular `G = H` (vacuous (2)) is a counterexample, and so is every padding
of `H` to any `n >= k`.

**Corollary F3 (monotonicity).** Let `Bad(n) := {p legal : Conj(n, p) fails}`
and `S_H := {p legal : FS(H, p)}`. Then

    Bad(n) = union over k = 2..n, over non-complete H on k nodes, of S_H,

and `Bad(n) ⊆ Bad(n+1)`. A counterexample on `k` nodes is a counterexample on
every `n >= k`.

---

## 7. Lemma O: counterexample sets are open [PROVED]

**Lemma O.** For fixed `(G, sigma)`, the set of legal `p` at which `sigma` is
an SFP of `G` is open in the `(eps, delta)` plane. So is `S_H` for each `H`.

*Proof.* The legal region `{eps > 0, delta > 0, delta(1-eps) > eps}` is open.
At an SFP point, Observation 0 gives `det(I - W_sigma) > 0`. That determinant
is a polynomial in `(eps, delta)`, so it is nonzero nearby. There
`x(p) = theta (I - W_sigma(p))^{-1} 1` is a rational, hence continuous,
function, so `x > 0` and `drive_k < 0` persist on a neighbourhood. Condition
(3) is equivalent (Routh–Hurwitz) to all Hurwitz determinants of the monic
characteristic polynomial being `> 0`. Those are strict polynomial
inequalities, so they persist too. ∎

Consequences:

* **(O1) Rational witnesses [PROVED].** If `S_H` is nonempty, it is a nonempty
  open set and contains points with rational `(eps, delta)`. Every
  counterexample can therefore be certified at a rational parameter by
  `ctln.analyze_support` alone. Algebraic numbers are needed to describe the
  *boundaries* of the region, never to prove that a counterexample exists.
* **(O2) A counterexample point is never "undecided" [PROVED].** At every
  `p in S_H`:
  - the characteristic polynomial of `-I + W` is monic with all roots in the
    open left half-plane, so all its coefficients are `> 0`;
  - all its Hurwitz minors are `> 0` (Routh–Hurwitz necessity);
  - `I - W` is nonsingular.

  So `ctln.analyze_support` returns `stable = True`, `BOUNDARY = False` and
  `agrees = False` there: a *decided* counterexample. An undecided record means
  `p` lies on the zero set of a coefficient, minor or determinant polynomial,
  and such a point is never in `S_H`. If that polynomial is not identically
  zero, its zero set is a curve with empty interior. If it *is* identically
  zero, the graph is stable nowhere. So, read pointwise, CLAUDE.md §5's "the
  counterexamples are in the undecided pile" is not right. Counterexamples sit
  in open regions that a finite set of samples can step over. They are not at
  the samples where a minor happens to vanish. The three-valued convention of
  §3 stays correct and conservative: at the exact point, a vanishing minor
  means (3) fails there, and "undecided" refers to the neighbourhood. I am not
  proposing to change it.
* **(O3) Dense-set results extend [PROVED].** Each nondegeneracy condition (a
  principal minor of `I - W`, or a column-replacement determinant) is a
  polynomial in `(eps, delta)`. None of these polynomials is identically zero:
  at `(eps, delta) = (1, -1)` every off-diagonal weight is 0, so `I - W = I` and
  each such determinant equals 1. So the degenerate parameters form a closed
  set with empty interior. If a theorem rules out non-clique SFPs at all
  *nondegenerate* parameters in an open set `U`, then by Lemma O it rules them
  out on all of `U`. This is what allows a literature result stated "for
  nondegenerate CTLNs" to be used as a statement about `FS(k, p)`.

**[SANITY] An empirical consequence of (O2) that matters for Stage B sizing.**
Take the same 3,000 random 6-node graphs, full support (`undecided_rate.py`):

| point `(eps, delta)` | `(beta, t)` | undecided |
|---|---|---|
| `(1/2, 3/2)` | `(4, 1)` | **59 / 3000** |
| `(1/2, 3/2 + 1/997)` | | 0 / 3000 |
| `(500/1001, 1501/1000)` | | 0 / 3000 |
| `(1/10, 19/9)` | | 0 / 3000 |
| `(1/100, 29/2500)` | | 0 / 3000 |

At `(1/2, 3/2)`:

* In all 59 undecided cases the first vanishing quantity is the cubic
  characteristic-polynomial coefficient `c_3`. That is the sum of the 20
  principal `3x3` minors of `K = I + 4N + J`, which take only 13 integer
  values at this point, all multiples of 4, so cancellation is frequent
  (`undecided_kinds.py`; `undecided_factor.py` shows each vanishing `c_3` as a
  polynomial in `(beta, t)` that is not identically zero).
* 56 of the 59 already **fail positivity**, so "not an SFP" is decided by (1)
  alone. `analyze_support` still reports them as BOUNDARY because it tests
  stability independently.

So the 2.8% undecided rate in CLAUDE.md §5 and `ANOMALIES.md` is an artifact of
the special sample point `(beta, t) = (4, 1)`. It should not be used to size
Stage B. Suggested ANOMALIES entry; I have not edited that file.

---

## 8. What this means for n = 6

Hypothesis (iii) for `k <= 4`:

| `k` | status of `FS(k, p)` at every legal `p` | source |
|---|---|---|
| 1 | vacuous | |
| 2 | holds | Lemma K2 [PROVED] |
| 3 | holds | Geneson arXiv:2607.21396 Thm 4.1 with `m = 3` gives `q < 1` [CITED] |
| 4 | holds | Curto–Geneson–Morrison 2024 (AAM 154, 102652) [CITED] |

About the citations:

* Geneson's Theorem 4.1 is stated for "K permitted and stable", which is
  exactly our (1) plus (3) for a support of size `m`, at every legal parameter,
  with no nondegeneracy assumption. I have **not** re-derived its proof.
* For `k = 4` I rely on Geneson's summary "known for … supports of size at
  most four". I have **not** read CGM 2024, so I cannot confirm its exact
  hypotheses. If it is stated only for nondegenerate CTLNs, (O3) closes the
  gap. If it carries a parameter restriction, it does not.
* [SANITY] `small_k_region2d.py` ran `region2d.py` (its CAD was not audited by
  me) on every digraph with `k = 3` and `k = 4`, full support:
  - all 15 non-complete 3-node graphs: `never`;
  - 214 of the 217 non-complete 4-node graphs: `never`;
  - the other 3 raised "common factor between L and Delta" inside `region2d`.
    Each of those 3 has a dominated node (checked), so Lemma D' kills
    positivity at every parameter.

  Net: no non-complete graph with `k <= 4` is ever stable and positive.
  `k = 2` is skipped by `region2d` (the `D = -L` coincidence of §5) and is
  covered by Lemma K2.

Geneson's Theorem 4.1 also gives `FS(k, p)` whenever
`q >= k - 2 - (k-3) eps/2` [CITED]: `q >= 3 - eps` for `k = 5`, and
`q >= 4 - 3eps/2` for `k = 6`.

**Conclusion [PROVED, modulo the two citations for k = 3, 4].** For every legal `p`,

    Conj(6, p)  <=>  FS(5, p) and FS(6, p),

and `Bad(6) = Bad_FS(5) ∪ Bad_FS(6)`, where `Bad_FS(k)` is the union of `S_H`
over non-complete `k`-node `H`. `Bad_FS(5) ⊆ {q < 3 - eps}` and
`Bad_FS(6) ⊆ {q < 4 - 3eps/2}`. `Bad_FS(6)` is nonempty: it contains Geneson's
points.

Practical consequences:

1. **For the truth of the conjecture, Stage 4 (size-5 supports inside `n = 6`)
   reduces to the `n = 5` full-support question on 9,608 graphs. No condition
   (2) is needed.**
   - If `FS(5, p)` holds at every legal `p` (which the 2-D CAD described in
     ANOMALIES would establish), then no graph of any size has a non-clique
     SFP of size 5, by Lemma R, and Stage 4 is empty.
   - If some 5-node `H` has `FS(H, p)`, then its padding is a 6-node
     counterexample at the same `p` (F2).
   - Stage 4 with condition (2) is needed only to *list* the `(G, sigma)` pairs
     that realise a counterexample. Corollary F2 describes that list: a down-set
     of admissible in-sets, decided by Lemma H(4).
2. The `n = 6` map is the union of the open 2-D regions `S_H` over
   non-complete `H` on 6 nodes (and on 5 nodes, if any). By (O1), every
   nonempty one has a rational witness.
3. Sweeping canonical representatives is enough (relabelling invariance, §2).

---

## 9. Lemma H: the rank-one structure of positivity, strengthened [PROVED]

Setting: any support `sigma` of size `m` (the full support is `sigma = [n]`).
Let `N` be the receiver-first nonedge matrix of `G|sigma`:
`N_ab = 1` iff `a != b` and `s_b -/-> s_a`, that is, `N = J - I - (A|sigma)^T`.
Define

    beta = 1 + delta/eps,   t = (1-eps)/eps,   C = I + beta N,   K = C + t J,   u := adj(C) 1.

**(H0)** `I - W_sigma = eps K`.
*Proof.* `eps K = eps I + (eps + delta) N + (1-eps) J`. Diagonal:
`eps + 1 - eps = 1`. Off-diagonal with `N = 0` (an edge `s_b -> s_a`):
`1 - eps = c_ab`. Off-diagonal with `N = 1`: `eps + delta + 1 - eps = 1 + delta = c_ab`. ∎

**(H1)** `adj(K) 1 = adj(C) 1 = u` for all `C` and `t`. No invertibility needed.
*Proof.* By Cramer, `(adj(M) 1)_i = det(M with column i replaced by 1)`. In `K`
with column `i` replaced by `1`, every other column is `C_j + t 1`. Subtracting
`t` times column `i` (which is `1`) from each of them leaves the determinant
unchanged and gives `C` with column `i` replaced by `1`. ∎

**(H2)** `det K = det C + t 1^T u`. This is the matrix determinant lemma.
*Proof.* The determinant is multilinear in the columns `C_j + t 1`. Every term
using `t 1` twice vanishes. The terms using it once sum to
`t sum_i det(C with column i -> 1) = t 1^T u`. ∎

**(H3) Positivity.** The following are equivalent:

* `I - W_sigma` is nonsingular and its solution `x` of (1) is `> 0`;
* every `u_i` is nonzero, and all `u_i` have the same sign `s`.

When they hold: `sign det C = sign det K = s`, and `x = (theta/eps) u / det K`.

*Proof.* The new ingredient is `C u = det(C) 1`, which is the adjugate identity.
`C` is entrywise nonnegative with unit diagonal, since `beta > 0` and
`N >= 0`. If `u > 0`, then `(Cu)_i = u_i + beta sum_j N_ij u_j >= u_i > 0`, so
`det C > 0`. Likewise `u < 0` gives `det C < 0`.

(<=) Since `t > 0` and `1^T u` has sign `s`, (H2) gives
`sign det K = sign det C = s`. So `K` is nonsingular, and
`x = (theta/eps) K^{-1} 1 = (theta/eps) adj(K) 1 / det K = (theta/eps) u / det K > 0`
by (H1).

(=>) `u = adj(K) 1 = det K · K^{-1} 1 = det K · (eps/theta) x`, which is
strictly one-signed with the sign of `det K`. ∎

What this means:

* **The task's form.** "All `u_i` nonzero with common sign equal to
  `sign det K`" is equivalent to the above. The clause "equal to `sign det K`"
  is automatic.
* **Geneson's Lemma 2.1 is confirmed.** "K permitted iff C permitted" holds in
  both directions, by (H3) applied with `t > 0` and with `t = 0`.
* **Permission does not depend on `t`.** It depends on `(eps, delta)` only
  through `beta`.
* **Under (3).** `det K > 0` (Observation 0 and H0), so
  `(1) and (3) <=> u > 0 entrywise and (3)`.
* [SANITY] `lemma21_search.py` found no one-signed `u` with a mismatched
  `det C` sign over all 9,845 digraphs on 2 to 5 nodes, before this proof was
  found. (H3) makes that search unnecessary.

**(H4) Condition (2) in beta-form.** If (1) holds and `k` is outside `sigma`, then

    drive_k = (theta / det K) · ( det C - beta · sum_{j in sigma, j -/-> k} u_j ),

so (2) at `k` holds iff `beta · sum_{j in sigma, j -/-> k} |u_j| > |det C|`.

*Proof.* Write `W_kj = -(1-eps) - (eps+delta) N_kj`, where `N_kj = 1` iff
`j -/-> k`. Substitute `x = (theta/eps) u / det K`:

    drive_k = theta - (theta/det K)[t 1^T u + beta sum_j N_kj u_j]
            = (theta/det K)[det C - beta sum_j N_kj u_j],

using (H2). By (H3), `u`, `det C` and `det K` all share one sign `s`. Multiply
through by `s`. ∎

**Consequences of H4:**

* (1) and (2) together depend on `(eps, delta)` only through `beta`. Only (3)
  depends on `t`. This is exactly the structure `fastregion.py` uses, and here
  it is proved.
* A target has an empty sum, so (2) fails since `|det C| > 0`. This agrees with
  Lemma T(c).
* `S = ∅` gives `beta 1^T u > det C`. That also follows from `(Cu)_i <= beta 1^T u`
  with strict inequality, because `u_i < beta u_i`. This agrees with
  Lemma T(b).

**(H5) Stability (remark).** Apply (H2) to `mu I + C` in place of `C`:
`det(mu I + K) = det(mu I + C) + t 1^T adj(mu I + C) 1`. This is the "linear in
`t`" structure recorded in ANOMALIES 2026-10-13 and used by `fastregion.py`.
Hurwitz stability of this polynomial in `mu` is equivalent to (3), since
`-I + W_sigma = -eps K` and `eps > 0` only rescales the roots.

[SANITY] `check_FGH.py`:

* `check_H0`: 300 cases.
* `check_H_symbolic`: (H1) and (H2) for a generic symbolic `C`, `m <= 4`.
* `check_H_numeric`: 500 cases.
* `check_H_strong`: (H3), (H4) and the beta-form of (2), 398 positive cases.

---

## 10. Status summary

**[PROVED] here:**

* Observation 0
* Lemmas G, R (with relabelling), R2, T, D', K2, O, H0–H4
* Theorem F: (i) <=> (ii) <=> (iii) pointwise, with corollaries F1, F2, F3
* the `n = 6` reduction, modulo the `k = 3, 4` citations

**[CITED], not re-derived:**

* Geneson Thm 4.1 (`k = 3`, and the high-`q` range for all `k`)
* Curto–Geneson–Morrison 2024 (`k = 4`); exact hypotheses not checked
* the identification of T(b), T(c) with CGM 2019's domination cases

**[SANITY], not proofs:**

* all `check_FGH.py` checks
* the `k = 3, 4` `region2d` sweep (CAD not audited by me)
* the undecided-rate measurements

**Not done:**

* a self-contained proof for `k = 3` or `k = 4`. A partial `k = 3` route:
  dominated graphs fail (1) by D'; symmetric non-complete graphs fail (3)
  because a non-adjacent pair gives a non-positive-definite principal `2x2`
  block `[[1, 1+delta], [1+delta, 1]]`; some graphs fail (3) because
  `c_2 = sum over pairs of the principal 2x2 minors of K` is `< 0`. Not every
  domination-free 3-node case is closed this way, so it is not claimed.
* no Lean formalisation yet.

**Flags for the orchestrator:**

1. The task's "condition (2) is never needed in (F1)" is false as stated (§6).
2. The 2.8% undecided rate is an artifact of the sample point `(1/2, 3/2)`
   (§7). 56 of 59 such records already fail positivity.
3. The 3 four-node graphs on which `region2d.py` raises ArithmeticError are all
   dominated. That is a `region2d` limitation when `L` shares a factor with a
   Hurwitz factor, not a mathematical gap.
