# Is sup_n q_n* infinite? (workflow wf8, label `prover_q`)

Date 2026-10-15. Conventions are those of CLAUDE.md §2: `A[i][j] = 1` iff `i -> j`; `D := J - I - A` (`D[i][j] = 1` iff
`i != j` and `i` does not send to `j`); `beta = 1 + delta/eps`, `r = 1/beta`, `gamma = 1 - r`, `t = (1-eps)/eps`,
`q = delta(1-eps)/eps = (beta-1) t/(t+1)`; legal iff `t(beta-2) > 1`. The repository was not modified. Every file named
below is in `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/` (below: `PQ/`),
and each script re-runs standalone with `python3 <script>`. The scripts import the repo read-only.

Labels: **PROVED** means a proof is written here, or the statement is cited from a repo note that marks it proved and
refereed. **CERTIFIED** means an exact integer/Fraction computation, with the script and output file named.
**HEURISTIC** means floating point used only to steer a search. **OPEN** means unresolved.

> **Hand-verification status (CLAUDE.md §3).** Every graph named here is a **CANDIDATE counterexample pending human
> hand-verification**. The smallest explicit ones (n = 28 and n = 41, §5.3, adjacency rows in
> `certificates_explicit.json`) are the natural first targets. Theorem U is a proof. It cites only lemmas that the repo
> marks PROVED and refereed (`notes/lemma_I_spectral.md` I.1, I.3, I.4, (R); `notes/research_wf3/theory/notes_theory.md`
> D1, D2, X1, X2). It does not use the fact that any particular graph is a counterexample.

---

## 0. Answer and summary

**Answer: yes. sup_n q_n* = infinity. More strongly, q_n* -> infinity as n -> infinity (Theorem U, PROVED).** For every
integer N >= 3 there is n0(N) such that for every n >= n0(N) some non-clique n-node digraph has a stable positive
full-support fixed point at legal (eps, delta) with q arbitrarily close to N - 1. Hence q_n* >= N - 1 for all
n >= n0(N). This is consistent with Geneson's upper bound q_n* <= n - 2, since n0(N) is far larger than N.

| # | statement | status | where |
|---|---|---|---|
| 1 | Lemma T: for an in-tree D with root rho, `chi_{D - e_rho 1^T}(z) = z^m Lam(1/z)` and `chi_D(z) = z^m (1 - Gam(1/z))` | PROVED (determinant lemma); CERTIFIED on 40 random in-trees + 5 forests | §2, `check_lemmas.py` |
| 2 | Lemma Q: `psi_D(z) = z^(m-1) Q(1/z)`, `Q = [(m+1-|N|)(1-Gam) + (|N|-1) Lam]/(1-x)` | PROVED (Sherman–Morrison); CERTIFIED (same) | §2 |
| 3 | Arms forest: any `p = z^d + z^(d-1) + c_2 z^(d-2) + ... + c_d` with integers `c_2 >= 1`, `c_j >= 0` (forest constraint) is realised as `chi_L = z^(m-d) p`, `m = 2p(1) - 1` | PROVED | §3 |
| 4 | Lemma P: closed-form positivity of the arms-forest core | PROVED; CERTIFIED vs direct exact solves (12 cases) | §3.2, `check_lemmas.py` |
| 5 | Lemma S (rounding): `p~` = coefficientwise-rounded `prod_k ((z+s)^2 + W^2 k^2)`, `s = 1/(2N)`: all roots in `-3s/2 < Re z < -s/2` once `W > W0(N) = 2N 4^N (N+1)^(2N-2)` | PROVED (Rouché); spot-CERTIFIED at W = W0(N), N = 3..8 | §4, `proof_ingredients.py` |
| 6 | **Theorem U**: `q_n* -> infinity` | **PROVED** (modulo the cited repo lemmas) | §4 |
| 7 | Instances of Theorem U: q -> 2, 3, 4 (N = 3, 4, 5), every n >= 1812, 6,214,995, 1,590,755,528 | CERTIFIED (reduced exact certificate, all k >= k0) | §5.1, `theorem_U_instances.json` |
| 8 | 65-node core `c = [1,1,11,7,10,3,1]`: q -> 3, 4, 9/2, **23/5** for every n >= 219, 657, 2475, **4909** | CERTIFIED (reduced; explicit Lemma T/Q/X1/positivity checks on the 65-node core and k = 1, 2) | §5.2, `certificates_forest_small.json` |
| 9 | Explicit graphs, checked directly by `ctln.analyze_support`: n = 28 at q = 20000/10001, n = 41 at q = 22500/10001 | CERTIFIED (E3 is `ctln` alone, no lemma used; E1/E2 use only D1/D2) | §5.3, `certificates_explicit.json` |
| 10 | Tool validation against known values: q_6* bracket (561/481, 7/6); F_1 (n=7) beta_c in (2.2746, 2.2748); plastic positivity edge of F_2; n=7 cores of Theorems Q, Q2 | CERTIFIED, ALL PASS | §6, `validate_small.json` |
| 10b | Route (a), circulant (normal, out-regular) cores alone: `alpha(A|perp) >= -1/2`, never limit-stable at legal beta | PROVED | §7 |
| 11 | Route (b), a universal bound `q_n* <= C` or `alpha(A|perp) >= -1 + c` for all non-complete positive digraphs | REFUTED by Theorem U: alpha(A|perp) < -1 + 1/N occurs for every N | §7 |
| 12 | Growth rate of q_n* in n; an explicit n0(N) | OPEN (our n0(N) is super-exponential and only partly explicit) | §7 |

---

## 1. Tools cited (PROVED in the repo, refereed)

- **D1** (notes_theory §1, from Lemma I.1): `A|perp = -I - D|perp`. So `alpha(A|perp) < -gamma` iff every root of
  `psi_D(z) := 1^T adj(zI - D) 1 = m chi_{D|perp}(z)` has `Re z > -r`.
- **D2**: `w = (A^T + gamma I)^{-1} 1 > 0` iff `u := (D^T + rI)^{-1} 1 > 0` and `U := 1^T u > 1`, and then `w = u/(U-1)`.
- **X1** (clique-class extension `G[S,T;k]`, k bidirected new nodes, `b = 1 - S`, `c = 1 - T`):
  `n chi_{D'|perp} = z^(k-1) Phi_k` with `Phi_k = z psi_D + k chi_L`, where `L = D - b 1^T - 1 c^T`.
- **X2** (T = all): `G[S,all;k]` is positive at beta iff `w_G(beta) > 0` and `sigma_S := sum_{i in S} w_G,i > 1`.
  This does not depend on k.
- **(R)** root continuity for monic polynomials (lemma_I_spectral §1).
- **Limit rule** (lemma_I_spectral §9, I.3 + I.4). Let beta > 2 be rational, `w(gamma) > 0` and `alpha(A|perp) < -gamma`.
  Then there is T such that the full support is a stable fixed point at the legal points
  `(eps, delta) = (1/(1+t), (beta-1)/(1+t))` for every t > T. These points have `q = (beta-1) t/(t+1) -> beta - 1`.

## 2. Lemmas T and Q (PROVED here)

Let D be a 0/1 zero-diagonal m x m matrix with a root rho. Every node u != rho has exactly one D-out-neighbour, and
iterating it reaches rho in `ell(u)` steps (`ell(rho) = 0`). Let `N := N^+(rho)` be the D-out-neighbourhood of rho, and
`x := 1/z`. Define
`Lam(x) = 1 + sum_{v not in N} x^(1+ell(v))` (rho itself contributes x) and `Gam(x) = sum_{w in N} x^(1+ell(w))`.

**Lemma T.** `chi_{D - e_rho 1^T}(z) = z^m Lam(x)` and `chi_D(z) = z^m (1 - Gam(x))`.

*Proof.* Write `D = D0 + e_rho 1_N^T`, where D0 is D with row rho zeroed. D0 is nilpotent (a forest of paths into rho), so
`det(zI - D0) = z^m` and `R0 := (zI - D0)^{-1} = sum_k x^(k+1) D0^k`, with `(D0^k)[u][v] = 1` iff v is the k-th successor
of u. Hence `(R0 e_rho)_u = x^(1+ell(u))`. By the matrix determinant lemma,
`chi_D = z^m (1 - 1_N^T R0 e_rho) = z^m (1 - Gam)`. Similarly `D - e_rho 1^T = D0 - e_rho 1_{N^c}^T` (rho is not in N), so
`chi_{D - e_rho 1^T} = z^m (1 + 1_{N^c}^T R0 e_rho) = z^m Lam`. ∎

**Lemma Q.** `psi_D(z) = z^(m-1) Q(x)`, where `Q(x) = [(m+1-|N|)(1 - Gam) + (|N|-1) Lam] / (1 - x)` (a polynomial).

*Proof.* Sherman–Morrison gives
`1^T (zI-D)^{-1} 1 = 1^T R0 1 + (1^T R0 e_rho)(1_N^T R0 1)/(1 - Gam)`.
The ingredients are `1^T R0 e_rho = sum_u x^(1+ell(u)) = Lam - 1 + Gam`,
`(R0 1)_u = x (1 - x^(ell(u)+1))/(1-x)`, `1^T R0 1 = x (m - Lam + 1 - Gam)/(1-x)` and `1_N^T R0 1 = x (|N| - Gam)/(1-x)`.
Multiply by `chi_D = z^m (1-Gam)`. The numerator divided by x is
`(1-Gam)(m+1-Lam-Gam) + (Lam-1+Gam)(|N|-Gam)`. The terms `-(1-Gam)(Lam+Gam) - Gam(Lam-1+Gam)` sum to `-Lam`, so the
numerator equals `(m+1-|N|)(1-Gam) + (|N|-1) Lam`. ∎

CERTIFIED: `check_lemmas.py` -> `check_lemmas.json` checks both identities exactly on 40 random in-trees (m = 4..25) and
on 5 arms forests.

## 3. The arms forest (PROVED)

### 3.1 Construction

Take nodes c0, c1, c2 with D-arcs `c0 -> c1 -> c2 -> c0`. Each **arm** a adds two nodes a, z_a with D-arcs `c2 -> a`,
`a -> z_a` and `z_a -> target(a)`. An arm of depth 1 targets c2, c1 or c0 (its "exit"). An arm of depth k >= 2 targets
the first arm of depth k-1 with the same exit. Then D is an in-tree rooted at `rho = c2` with
`N = {c0} ∪ {arms}`, and Lemma T gives
`Lam(x) = 1 + x + x^2 + sum_a x^len(a)`, where `len(a) = 2 depth(a) + (0, 1, 2)` for exits (c2, c1, c0).

So for integers `c_2 >= 1` and `c_j >= 0`, with `c_2 - 1` arms of length 2 and `c_j` arms of length j, we get
`chi_{D - e_c2 1^T}(z) = z^(m-d) p(z)` with `p(z) = z^d + z^(d-1) + c_2 z^(d-2) + ... + c_d`. The forest constraint is
that an arm of length j >= 5 needs some arm of length j-2. The core has `m = 3 + 2(p(1) - 2) = 2p(1) - 1` nodes.
(`forest_exact.build_core`.)

**Sink class.** `G_k := core[S, all; k]` with `S` = all nodes except c2. That is, k new nodes, bidirected among
themselves, which send to every core node and receive from every core node except c2. Then `b = e_c2` and `c = 0`, so
`L = D - e_c2 1^T`. By X1 and Lemmas T and Q,

    n chi_{D'|perp} = z^(k-1) Phi_k,   Phi_k = z^(m-d) (G + k p),   G(z) := z^d Q(1/z)   (degree d, leading coeff m).

Note that `c_1 = 1` is forced: `trace L = -1`. So the d nonzero limit roots sum to -1, and the ceiling of a degree-d core
is `beta* < d`. This is the "many-class ceiling" of notes_theory §5.5, but the core size is now unbounded.

### 3.2 Lemma P (positivity, closed form; PROVED)

Write the system `(D^T + rI) u = 1` row by row as "D-in-mass of i + r u_i = 1".
- Node z_a has D-in-set {a}, so `u_{z_a} = (1 - u_a)/r`.
- Arm a has D-in-set `{c2} ∪ {z_b : target(b) = a}`.

Put `lam := 1 - (1 - u_c2)/r`. Induction from the leaves gives `u_a = 1 - lam f_a` and `u_{z_a} = lam f_a / r`, where
`f_leaf = 1` and `f_a = 1 + r^-2 sum_{children b} f_b`. Summing f over the depth-1 arms of each exit gives

    F2 = c_2 - 1,   F1 = sum_{odd j >= 3} c_j r^(3-j),   F0 = sum_{even j >= 4} c_j r^(4-j).

The three cycle rows are then linear in `e_i := 1 - u_ci` and lam:
`e2 = r(1 - lam)`, `e1 = r(1-e2) + lam F2/r`, `e0 = r(1-e1) + lam F1/r`, `e2 = r(1-e0) + lam F0/r`. The unique solution is

    lam = r^3 (1 - r + r^2) / Den,   Den = F0 - r F1 + r^2 F2 + r^2 + r^5,

and sympy confirms the residual is 0. The core is positive and the class is activated (X2) iff
`lam > 0`, `lam * fmax < 1`, `u_c0, u_c1 > 0` and `0 < u_c2 < 1`. Here `fmax <= max(F0, F1, 1)`, and `U > 1` is
automatic because every arm contributes `u_a + u_{z_a} >= 1`. CERTIFIED: the closed form equals the direct exact solve
in 12 cases (`check_lemmas.json`), and the 65-node core's direct solve agrees (`certificates_forest_small.json`, C4).

## 4. Theorem U (PROVED)

**Theorem U.** For every integer N >= 3 there is n0(N) such that for every n >= n0(N) there is a non-clique digraph on
n nodes whose full support is a stable positive fixed point at legal points `(eps, delta)` with
`q = (N-1) t/(t+1)` for all large t. Hence `q_n* >= N - 1` for n >= n0(N), and `lim_n q_n* = sup_n q_n* = infinity`.

*Proof.* Fix N >= 3, `s = 1/(2N)` and `r = 1/N`.

**Step 1 (target polynomial).** `p_W(z) = prod_{k=1..N} ((z+s)^2 + W^2 k^2)`, with W a positive integer. It has degree
d = 2N, `c_0 = 1`, `c_1 = 2Ns = 1`, and every coefficient is positive. Let `p~` have `c~_0 = c~_1 = 1` and
`c~_j = nint(c_j)` for j >= 2, so `|c~_j - c_j| <= 1/2`. For j >= 2 we have `c_j >= 1`: for j = 2i the term
`e_i(W^2, ..., W^2 N^2) >= 1` appears, and for j = 2i+1 the term `e_i(a)(N-i)/N >= C(N,i)(N-i)/N = C(N-1,i) >= 1` appears (`a_k = W^2 k^2`). So
`c~_j >= 1` for every j >= 2, the arms forest of §3.1 exists, and the forest constraint holds.

**Step 2 (Lemma S, Rouché).** Scale `z = W zeta`. Then `P(zeta) := W^-2N p_W(W zeta) = prod((zeta + s/W)^2 + k^2)`, with
the 2N simple roots `zeta_{k,±} = -s/W ± ik`, pairwise at distance >= 1. The rounded polynomial is `P + E` with
`E(zeta) = sum_{j>=2} (c~_j - c_j) W^-j zeta^(2N-j)`, so `|E| <= W^-2 (N+1)^(2N-2)` on `|zeta| <= N+1`. On the circle
`|zeta - zeta_{k,±}| = rho := s/(2W)` we have `|P| >= rho 2^-(2N-1)`. Rouché applies once
`W > W0(N) := 2N 4^N (N+1)^(2N-2)`. Then each disc contains exactly one root, the 2N disjoint discs contain all of them,
and in z every root of `p~` satisfies `|z - (-s ± iWk)| < s/2`. (The exact spot check at W = W0(N) itself also passes,
`proof_ingredients.json`.) In particular `-3/(4N) < Re z < -1/(4N)`.

**Step 3 (positivity at beta = N).** Take W > W0(N) and let it grow. `c~_2N >= W^2N (N!)^2 - 1/2` enters F0 with weight
`r^(4-2N)`. Every odd coefficient `c~_j` (j <= 2N-1) is `O(W^(2N-2))` and `F2 = O(W^2)`, so `F1/F0, F2/F0, 1/F0 -> 0`.
Hence `Den/F0 -> 1`, `lam -> 0+` and `lam F0 -> r^3(1-r+r^2) < 1`, while `lam F1, lam F2 -> 0`. Then
`u_c2 -> 1-r`, which lies in (0,1); `u_c1 -> 1 - r + r^2 > 0`; and `u_c0 -> 1 - r + r^2 - r^3 > 0`. So for all large W,
Lemma P gives `w_core(N) > 0` and `sigma_S > 1`.
Fix one such `W > W0(N)`. The core has `m = 2 p~(1) - 1` nodes.

**Step 4 (sink class).** By X2, `G_k` is positive at beta = N for every k >= 1. `Phi_k/(m+k)` is monic and converges
coefficientwise to `chi_L = z^(m-d) p~` as k -> infinity. The roots of `chi_L` have `Re > -3/(4N)` (or are 0). By (R)
with radius `< 1/(4N)`, there is k0 such that for every k >= k0 every root of `Phi_k` has `Re > -1/N = -r`. The extra
factor `z^(k-1)` contributes roots with Re 0. By D1, `alpha(G_k|perp) < -gamma(N)`.

**Step 5.** The limit rule (§1) at the rational beta = N > 2 gives an SFP at legal points with `q = (N-1) t/(t+1)` for
all `t > T`. G_k is not a clique, because its core has non-edges. Take `n0(N) = m + k0`. ∎

Remarks.
- W0(N) and the core are explicit. k0 comes from (R) and is not explicit in general. For the instances in §5 it is
  computed exactly.
- At fixed beta, the proof uses Lemma I.3 only. The open "beta -> infinity jointly with t" case (Remark I.3.5) is not
  needed, because each N is handled at the fixed rational beta = N.
- The construction reduces "large q" to a question about integer polynomials: a polynomial with c_1 = 1, nonnegative
  integer coefficients and all roots in `Re z > -1/beta`. Positivity is automatic for "top-heavy" coefficients. This is
  the reformulation that removes the fixed-core ceiling `beta < m`: the core grows with N.

## 5. Certified instances

### 5.1 Theorem U instances (CERTIFIED, reduced certificate; `theorem_U.py` -> `theorem_U_instances.json`; re-run identical, `rerun/`)

Certificate contents: (C1) exact Hurwitz test that every root of `p~` has `Re > -1/beta*`; (C2) Lemma P closed form; and
(C3) every Hurwitz minor of `(-1)^d (G + k p~)(-1/beta* - x)` is a polynomial in k with positive leading coefficient and
no real root >= k0 (exact real-root isolation). Together these mean "for all k >= k0".

| N | W | beta* | q -> | d | core m | k0 | n0 | time |
|---|---|---|---|---|---|---|---|---|
| 3 | 1 | 3 | 2 | 6 | 257 | 1,555 | 1,812 | 0.3 s |
| 4 | 2 | 4 | 3 | 8 | 442,233 | 5,772,762 | 6,214,995 | 2 s |
| 5 | 2 | 5 | 4 | 10 | 44,039,851 | 1,546,715,677 | 1,590,755,528 | 54 s |
| 6 | 2 | 6 | 5 | 12 | 6,323,344,851 | not computed | | C1, C2 pass (`theorem_U_N6_partial.json`); W = 1 fails C1; the sympy k0 step was stopped at the time budget. Existence of k0 follows from §4 Step 4. |

### 5.2 A 65-node core (CERTIFIED; `certify_forest.py` -> `certificates_forest_small.json`)

The core is `c = [1,1,11,7,10,3,1]` (d = 6, m = 65). C1–C3 pass at each beta*, and at beta* = 5 the explicit checks
also pass: Lemma T and Lemma Q on the 65 x 65 matrix, X1 on the explicit G_1 and G_2, and the direct exact positivity
solve.

| beta* | 4 | 5 | 11/2 | 28/5 |
|---|---|---|---|---|
| q -> | 3 | 4 | 9/2 | **23/5 = 4.6** |
| k0 | 154 | 592 | 2,410 | 4,844 |
| n0 = 65 + k0 | 219 | 657 | 2,475 | 4,909 |

**So q_n* >= 23/5 for every n >= 4909.** This beats Theorem Y's 3.8 (n = 8 + 5k >= 1823) and its family limit of 3.96.
The ceiling for d = 6 is beta* < 6.

### 5.3 Explicit graphs checked by ctln (CERTIFIED; `certify_explicit.py` -> `certificates_explicit.json`)

The `ctln` check (E3) uses no lemma at all; (E1) and (E2) use only D1 and D2. Each record contains: (E0) the explicit psi of the n x n non-edge matrix equals
`z^(k-1) z^(m-d)(G + kp)`; (E1) the explicit limit test; (E2) explicit positivity; (E3) `ctln.analyze_support` on the full
support.

| c | k | n | beta | t | eps | delta | q | result |
|---|---|---|---|---|---|---|---|---|
| [1,1,5,3,1] | 9 | 28 | 3 | 10^4 | 1/10001 | 2/10001 | 20000/10001 | legal, x > 0 (min 110011/4235068), all Hurwitz minors positive, not a clique |
| [1,1,4,2,1] | 26 | 41 | 13/4 | 10^4 | 1/10001 | 9/40004 | 22500/10001 | legal, x > 0, all Hurwitz minors positive, not a clique |

Adjacency rows are stored in the JSON. These two graphs are the recommended hand-verification targets.

**Negative control** (`negative_control.py` -> `negative_control.json`, CERTIFIED). On the explicit graphs the reduced
threshold is exactly sharp. At k = k0 - 1 (n = 27, resp. 40) the explicit limit test fails with an exact negative
Hurwitz minor (minor 2, resp. 3); at k0 and k0 + 1 it passes; positivity holds at all three k, as X2 predicts. This is an
independent consistency check between the reduced polynomial `G + k p` (Lemmas T, Q, X1) and the explicit n-node graphs. They were found
by `x_small_forest.py` (HEURISTIC float ranking, then exact `certify`).

## 6. Validation against known numbers (CERTIFIED; `validate_small.py` -> `validate_small.json`, ALL PASS)

These use the same generic primitives as above (`exact_core.py`: psi via integer charpolys, shifted Hurwitz test,
exact positivity).
- V1: n = 6, Geneson's G1. Limit-stable and positive at q = 561/481; not limit-stable at q = 7/6. This brackets
  q_6* = 1.1664.
- V2: n = 7, `F_1`. beta_c lies in (2.2746, 2.2748), matching Theorem F's 2.2747.
- V3: n = 8, `F_2`. Positive at 2.3247 but not at 2.3248, a sign change of `b^3 - 3b^2 + 2b - 1` (the plastic number + 1,
  i.e. q_7* = 1.3247). It is still limit-stable at 2.3248.
- V4: the n = 7 core `&F^S{^ReL|?` (Theorem Q) has `chi_L = z^3 (z^4+z^3+3z^2+z+1)` and is positive at beta = 3 and 10.
- V5: the n = 7 core `&F^U|m^QTq?` (Theorem Q2) has `4 chi = z^3 (2z^2+z+2)^2`.

## 7. Open questions and the other routes

- **Route (b), a universal bound: REFUTED.** By D1, the graphs of Theorem U have `alpha(A|perp) < -1 + 1/N`. So no
  bound `alpha >= -1 + c` holds over non-complete digraphs with positive w.
- **Route (a), circulant cores K_m minus a pattern: PROVED to fail on their own.** Let the non-edge matrix D be nonzero,
  normal and out-regular (`D1 = d 1`); every circulant pattern qualifies. Then `alpha(A|perp) >= -1/2`, so the core is
  never limit-stable at a legal beta (I.3(d)). *Proof.* D is unitarily diagonalisable, so the real parts of its
  eigenvalues are the eigenvalues of `H = (D + D^T)/2`. A nonzero D gives a 2 x 2 principal block `[[0,w],[w,0]]` of H
  with `w >= 1/2`, so by interlacing some eigenvalue zeta of D has `Re zeta <= -1/2` (notes_theory §5.4, normal case).
  Since 1 is an eigenvector of D (and of D^T, by normality), `1^perp` is D-invariant. The excluded eigenvalue is
  `d >= 0`, so zeta lies in `spec(D|perp)`. By D1, `-1 - zeta` is an eigenvalue of `A|perp`, so
  `alpha(A|perp) >= -1 - Re zeta >= -1/2`. ∎ So symmetric patterns on K_m need a class or another symmetry-breaking device. The
  in-tree reformulation made this route unnecessary: for a sink class over an in-tree the limit polynomial is free (§3.1).
- **OPEN: the growth rate.** Our n0(N) is at least `m >= 2 p~(1) - 1 >= 2 (N!)^2 W^(2N)`, which is super-exponential in N.
  Certified points: q -> 2 at n = 19 (Theorem Q2), q -> 4.6 at n = 4909 (§5.2), q -> 4 at n ≈ 1.6e9 (§5.1). Is
  `q_n* >= c log n`, or `>= n^c`? Upper bound: `n - 2` (Geneson). A sharper construction would need polynomials with
  c_1 = 1 and roots in `Re > -1/N` whose value p(1) (hence core size) is much smaller. Integer polynomials with small
  coefficients and roots in a thin left strip are the bottleneck. Also open: an explicit k0 bound (quantitative (R)) for
  general N.
- **HEURISTIC record (previous attempt, floats only, not certificates):** `arms.log`, `grow3.log` and `hz3.log` hold the
  H/Z-design and arm-growth searches that led to the arms forest. `x_*.py` are search scripts; their outputs are not
  claims.

## 8. Items for ANOMALIES.md (not appended: the repository was not to be modified)

1. The trace "many-class ceiling" `beta < m` (notes_theory §5.5) is real but not an obstruction. In-tree cores realise
   every admissible limit polynomial, and the ceiling grows with the core. The obstruction lives in the degree
   (`beta* < d`, since `c_1 = 1` is forced by `trace L = -1`), not in the graph.
2. Positivity, the binding constraint for every small-n record (plastic number, 2.7693, 3.023), becomes free in the arms
   forest: `lam -> 0` when the deepest coefficient dominates. All the difficulty moves to the spectral side.
3. The k0 thresholds are 2.4–75 times the core size (154/65 up to 4844/65; 6, 13 and 35 for N = 3, 4, 5). Counterexamples with large
   q are therefore mostly one huge clique class, a near-complete graph. This is consistent with I.5 (`alpha = -1` only for the
   complete graph).

## 9. Files (PQ/)

- Exact libraries: `forest_exact.py` (construction, Lemma T check, Hurwitz strip test, direct positivity) and
  `exact_core.py` (generic psi / limit / positivity / ctln witness).
- Certificate scripts and outputs:
  - `check_lemmas.py` -> `check_lemmas.json` (Lemmas T, Q, P)
  - `certify_forest.py` -> `certificates_forest_small.json` (§5.2)
  - `theorem_U.py` -> `theorem_U_instances.json` (§5.1)
  - `theorem_U_N.py N` -> `theorem_U_instance_N<N>.json` (N = 6 was not finished; partial C1/C2 record in
    `theorem_U_N6_partial.json`)
  - `proof_ingredients.py` -> `proof_ingredients.json`: Lemma S strip and positivity at W0(N), N = 3..8, all pass
  - `certify_explicit.py` -> `certificates_explicit.json` (§5.3)
  - `validate_small.py` -> `validate_small.json` (§6)
  - `negative_control.py` -> `negative_control.json` (§5.3)
- Re-run logs: `rerun/` (check_lemmas, certify_forest and theorem_U outputs, all reproduced identically).
- HEURISTIC (floats steer only): `heur_lib.py`, `hz_lib.py`, `forest.py`, `x_*.py`, `hlevel.py`, `hopt.py`, `*.log`.
