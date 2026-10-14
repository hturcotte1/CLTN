# Lemma M: disproof, what is true instead, and what the n=7 finite-eps map can use

Prover notes, 2026-10-14. Conventions are exactly CLAUDE.md §2 and notes/lemmas.md §0. Nothing in the
repository was modified. Every script named below is in this directory and re-runs standalone with
`python3 <script>` (they import the repo read-only from /home/user/CLTN).

Labels: **PROVED** (proof here), **CERTIFIED** (exact rational computation, script named), **HEURISTIC**
(floating point, used only to search; never a decision), **OPEN**.

## 0. Summary

| statement | status |
|---|---|
| M_s (T(beta) upward closed), pure stability, as stated in the task | **FALSE**, CERTIFIED at m = 3 (§2.1) |
| M_w (T(beta) nonempty => alpha <= -(1-1/beta)), pure stability | **FALSE**, CERTIFIED at m = 3 (§2.1) |
| M_s and M_w with the positivity hypothesis a(beta) > 0 (domination-free, SFP) | **FALSE**, CERTIFIED at n = 7 (§2.2) |
| graph level, pure stability: legal-stable somewhere => alpha <= -1/2 | **FALSE**, CERTIFIED at n = 6 (§2.3), dominated graphs |
| graph level, domination-free: legal-stable somewhere => alpha < -1/2 | **FALSE**, CERTIFIED at n = 7 (§2.4: `&FNALh\Uwp?`, `&FM?_~LVd~?`; neither is an SFP anywhere). Holds by exact exhaustion for n <= 6 (n <= 5: only K_n; n = 6: 28 graphs, all alpha < -1/2) |
| graph level with positivity (what the application needs), M_graph^+: SFP somewhere legal => alpha <= -1/2 | true for n <= 6 (existing exact maps); n = 7 **OPEN**, no counterexample in the exact searches of §2.4 |
| Identity (E): stability depends only on nu_A(c) = min Re spec(cJ - A) | PROVED (§1) |
| Lemma R (ray reduction): legal-stable somewhere <=> stable somewhere on the curve t(beta-2) = 1 | PROVED (§3), implemented exactly (curve1d.py), validated |
| Strip decomposition: full-support SFP is a union of 1-D problems | PROVED (§3.3) |
| Lemma AS: exact large-t asymptotics near the limit boundary; sign criterion for M near (beta_c, inf) | PROVED (§4); 3 graphs CERTIFIED against the repo engine |
| e_2 closed form; legal stability => more than half of the ordered pairs are edges | PROVED (§5) |
| In-/out-regular graphs: M_s holds | PROVED (§5) |

Bottom line for the n=7 map: **Lemma M cannot be used, but it is not needed.** Lemma R gives an exact,
sound 1-D prefilter: about 0.035-0.045 s per 7-node graph in pure Python under load. On 16,000 random n=7
positivity survivors it eliminated all 16,000. Survivors go to region2d. Total cost is about 16-40 CPU-hours (§6).

---

## 1. One identity that organises everything (PROVED)

Full support, `A[i][j] = 1` iff `i -> j`, `N = J - I - A^T`, `C = I + beta N`, `K = C + t J`. Put

    c := 1 + t/beta  (> 1),      gamma := 1 - 1/beta  (in (1/2, 1) for beta > 2).

**(E)** `K = beta (c J - A^T - gamma I)`.
*Proof.* `K = I + beta(J - I - A^T) + tJ = (1-beta) I - beta A^T + (beta + t) J`. Factor out beta. ∎

Since `J = J^T`, `spec(cJ - A^T) = spec(cJ - A)`. Hence

    K positive stable  <=>  nu_A(c) > gamma,     nu_A(c) := min Re spec(cJ - A).          (1)

`nu_A` depends only on the graph and on c. By Lemma C / C' / I (notes/lemmas.md), as c -> inf one eigenvalue
of cJ - A is real and tends to +inf, and the other m-1 converge to spec(-A|perp); so
`nu_A(inf) := lim nu_A(c) = -alpha`, `alpha = alpha(A|perp)`.

**Legality in these coordinates.** `t = beta(c-1)`, `beta = 1/(1-gamma)`, so `t(beta-2) = (c-1)(2gamma-1)/(1-gamma)^2`.
Thus `t(beta-2) > 1 <=> c(2 gamma - 1) > gamma^2`, and for c > 1:

    legal  <=>  gamma > nu_P2(c) := c - sqrt(c^2 - c),                                        (2)

where P2 is the 2-node digraph with one edge: `cJ - A_P2 = [[c, c-1],[c, c]]` has eigenvalues `c ± sqrt(c(c-1))`.
(`nu_P2` decreases from 1 to 1/2 on (1, inf); the root `c + sqrt(c^2-c) > 1` never binds.) So

    legal and stable at (c, gamma)   <=>   nu_P2(c) < gamma < nu_A(c).                          (3)

The legal boundary of the CTLN is literally the stability boundary of a one-way edge.
In this language M_s at beta says `{c > gamma^2/(2gamma-1) : nu_A(c) > gamma}` is upward closed, and M_w says
that if it is nonempty then `nu_A(inf) >= gamma`. Both are statements about the monotonicity of the single
function `nu_A`, and that is exactly what fails below.

---

## 2. Counterexamples (CERTIFIED)

### 2.1 Pure stability, m = 3 (`certify_pure_stability_cex.py`)

    A = [[0,1,1],
         [0,0,1],
         [1,0,0]]          beta = 22/5 throughout

| t | eps, delta | q | legal | stable (ctln.analyze_support) |
|---|---|---|---|---|
| 11/25 | 25/36, 85/36 | 187/180 | yes | **True** (all Hurwitz minors positive) |
| 1 | 1/2, 17/10 | 17/10 | yes | False |
| 100 | 1/101, 17/505 | 340/101 | yes | False |

`charpoly(A|perp) = z^2 + (4/3) z + 2/3`, so `alpha = -2/3 > -(1 - 5/22) = -17/22`. M_s and M_w both fail.
Positivity fails (`x_1 < 0`; node 1 is dominated by node 2), so this is not a CTLN counterexample of any kind.
Mechanism: `nu_A(c)` decreases monotonically from 1 to 2/3, like `nu_P2`.

### 2.2 With positivity, domination-free, n = 7 (`certify_M_violation.py`, `crosscheck_sympy.py`, `fastregion_refine.py`)

Graph `&FK@dyUfXp?` (digraph6, bit (i,j) row-major = 1 iff i -> j), domination-free, not complete:

    A = [[0,0,1,1,0,0,0],
         [0,0,0,0,1,1,0],
         [0,1,0,1,1,1,1],
         [0,1,0,0,1,0,1],
         [1,0,1,0,0,1,1],
         [1,0,1,1,0,0,1],
         [1,1,0,0,0,1,0]]          beta = 10561/5000 throughout

- t = 10: eps = 1/11, delta = 5561/55000, q = 5561/5500 > 1. ctln.analyze_support: x > 0, **stable**
  (all Hurwitz minors positive), stable_fixed_point = True.
- t = 1000: eps = 1/1001, delta = 5561/5005000, legal, x > 0, **not stable** (Hurwitz minor 6 negative).
- Limit: Hurwitz minor signs of `omega(.; beta)` are `[+,+,+,+,-,-]`, and `omega(0) = s != 0`. So omega is
  not Hurwitz and has no root on iR (b3), hence `alpha > -(1 - 1/beta)` strictly.

So T(beta) is nonempty, bounded above, and M_w fails, all with a(beta) > 0.

Three independent exact confirmations:
1. ctln.analyze_support (above);
2. a separate sympy implementation (adjugate, det(mu I + K), Hurwitz minors) gives the same signs at t = 10, 1000, 10^5;
3. the repo engine `fastregion.region` at t = 10, 30, 100, 1000, 10^4 returns the SFP beta-interval. Its right end,
   refined to rational isolating intervals of width < 1e-13, is **strictly above** the t = INF right end beta_c
   at every t (decimals for orientation only):

| t | right end of SFP beta-interval |
|---|---|
| INF | 2.1121280067227 (= beta_c) |
| 10 | 2.1130764373736 |
| 30 | 2.1124079688000 |
| 100 | 2.1122073871008 |
| 1000 | 2.1121357542928 |
| 10^4 | 2.1121287795343 |

Second instance: `&FMEHz\vpv?`, right end above beta_c at t = 30, 100, 1000, 10^4.

> **FLAG for hand-verification.** The t = 10 point is a legal, finite-eps, non-clique stable fixed point, i.e. a
> conjecture violation. The graph is already in `certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl`, and
> so is `&FMEHz\vpv?`. This is not a new counterexample graph; it is a new finite-eps witness for an existing
> candidate. Nothing below is built on it. The disproof of Lemma M as stated rests on §2.1 alone.

### 2.3 Graph level, pure stability, n = 6 (`certify_graph_level_pure.py`)

`&EODUpt{` and `&EIKDw\{` are stable at beta = 9/4, t = 11/2 (eps = 2/13, delta = 5/26, q = 55/52, legal).
`omega(.; 2)` has Hurwitz minor signs `[+,+,+,-,-]` and `s(2) = 4`, so `alpha > -1/2` strictly: they are limit-stable
at no beta. Both have dominated nodes, so positivity fails.
HEURISTIC count (float grid): 24 such graphs at n = 6. Some have `alpha` near -0.27 with a positive margin, so pure
legal stability can sit far from limit stability. Positivity, or domination-freeness, is essential.

### 2.4 Domination-free is not enough: positivity is needed at the same beta (CERTIFIED; `domfree_exact_all.py`, `check_hits_exact.py`, `certify_all.py` C4)

`domfree_exact_all.py` gives an EXACT verdict for every graph: the domination test, then the curve test of §3,
then the exact Hurwitz test of `omega(.;2)` for each graph that is legal-stable somewhere.

- n = 3, 4, 5, all graphs: **the only domination-free graph that is legal-stable anywhere (pure stability) is the
  complete graph** (n=3: 10 dominated, 5 never, 1 = K3; n=4: 171/46/1; n=5: 8174/1433/1).
- n = 6, exhaustive and exact over all 1,540,944 graphs (`domfree_exact_n6_{a,b}.{log,jsonl}`):
  - 1,336,566 graphs are dominated;
  - of the 204,378 domination-free graphs, 204,350 are never legal-stable and 0 escalate;
  - exactly **28** are legal-stable somewhere, **all with alpha < -1/2**, each with a single unbounded stable
    t-interval on the curve. They are the same 28 that the float pre-selection `n6_domfree_legal_stable.py` found,
    and they include K6 and the three known limit-stable graphs.

  So at n = 6, domination-free + legal-stable => limit-stable somewhere. That is stronger than what the existing
  n=6 map needed.
- **n = 7: the domination-free graph-level statement is FALSE.** Two counterexamples were found by a targeted search
  (`targeted_near_boundary.py`). The search uses the first-order mechanism of §4: alpha slightly above -1/2 and
  R > 1/8. Floats only selected the candidates; every decision below is exact.

  `&FNALh\Uwp?` (24 edges, domination-free, a 'U' line of the positivity filter):

      [[0,0,1,1,1,1,0],[0,0,0,1,0,0,0],[1,1,0,1,1,0,1],[0,0,1,0,1,1,1],
       [0,1,0,1,0,1,1],[0,1,1,1,0,0,0],[1,1,0,0,0,1,0]]

  - stable at the legal point t = 19/2, beta = 2 + 2/19 + 10^-6 (ctln.analyze_support);
  - `omega(.;2)` minor signs `[+,+,+,+,-,-]`, s(2) = 143, so alpha > -1/2 strictly;
  - region2d: positivity set empty, so never an SFP.

  `&FM?_~LVd~?` (26 edges, domination-free, a 'U' line resolved to 'R'):

      [[0,0,1,1,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,1,1,1],[1,1,1,0,0,1,1],
       [0,1,0,1,0,1,1],[1,1,0,0,1,0,1],[1,1,1,1,1,1,0]]

  - stable at the legal point t = 21/2, beta = 2 + 2/21 + 10^-6;
  - `omega(.;2)` minor signs `[+,+,+,+,-,-]`, s(2) = 11, so alpha > -1/2 strictly;
  - its positivity set is the single cell with endpoints in [16/7, 51/22], roughly beta in (2.29, 2.32), while legal
    stability lives near beta of about 2.1;
  - region2d (exact, positivity included): never an SFP.

  On the curve, each of these graphs has a single **bounded** stable t-interval: about (9, 10) and (10, 11).
  Every legal-stable domination-free graph found at n <= 6, and every n=7 limit-stable graph, has an unbounded one.

So the only conjecture that survives every test is the coupled one:

> **M_graph^+ (OPEN).** If the full support is an SFP at some legal finite (eps, delta), i.e. positive AND stable
> at the same (beta, t), then alpha(A|perp) <= -1/2.

Any proof must couple positivity and stability at the same beta.
- Stability plus legality alone fails at n = 3 (§2.1) and n = 6 (§2.3).
- Adding domination-freeness still fails at n = 7 (above).
- Positivity at a different beta fails too (`&FM?_~LVd~?`).

**Exact search for finite-eps-only counterexamples in the most likely window**
(`alpha_hist.py`, `curve_window.py`, `region2d_on_hits.py`). These are graphs that are SFP at a legal finite point
but never limit-stable.
- Floats (selection only) picked every dense (edges >= 22) n=7 positivity line with alpha in (-1/2, -0.3): 34,176
  graphs.
- The exact curve test returned `never` for 34,174 of them. The 2 that are legal-stable are exactly the graphs above.
- region2d says neither is an SFP anywhere.

So there is no counterexample to M_graph^+ in this window. Outside it (alpha >= -0.3, about 4.0M dense lines) is not
searched; §6 gives the cost of searching it.

---

## 3. Lemma R: reduction of legal stability to one curve (PROVED)

**Remark F** (notes/lemmas.md §7): `K(lambda beta, lambda t) = (1-lambda) I + lambda K(beta, t)`, so stability at
(beta, t) gives stability at (lambda beta, lambda t) for every lambda in (0, 1]. In (c, gamma) coordinates this says
the stable set is the hypograph `{gamma < nu_A(c)}`, which (1) also gives directly.

### 3.1 Statement

**(R)** The full support is stable at some legal (eps, delta) **iff** there is t > 0 with `K(2 + 1/t, t)` positive stable.

*Proof.*
- (<=) Stability is open in (beta, t), since eigenvalues depend continuously on the entries. So `K(2 + 1/t + eta, t)`
  is stable for small eta > 0, and that point is legal.
- (=>) Let (beta, t) be legal and stable, so beta > 2. Put `f(lambda) := lambda t (lambda beta - 2) - 1`. Then
  `f(2/beta) = -1 < 0` and `f(1) = t(beta-2) - 1 > 0`. By the IVT there is lambda* in (2/beta, 1) with
  `f(lambda*) = 0`. Then `t' := lambda* t` and `beta' := lambda* beta` satisfy `t'(beta' - 2) = 1`, so
  `beta' = 2 + 1/t'`. Remark F makes `K(beta', t')` stable. ∎

**(R+)** If the full support is SFP at a legal (beta0, t0), the curve point from (R) has `beta' = lambda* beta0 < beta0`.
So: SFP somewhere => there is t' with `2 + 1/t' < sup(Pos)` and `K(2 + 1/t', t')` stable, where Pos is the
positivity set in beta.

### 3.2 Exact implementation (`curve1d.py`)

Along the curve, `t K(2 + 1/t, t) = L(t) := t I + (2t+1) N + t^2 J`, an integer polynomial matrix. For t > 0 it has the
same stability as K.

- `p(mu, t) := det(mu I + L(t)) = sum_k c_k(t) mu^(m-k)`, with `deg_t c_k <= k+1`.
  *Proof:* `c_k` is the sum of the k x k principal minors. For each principal block, (T2) gives
  `det(L_S) = det(tI + (2t+1) N_S) + t^2 1^T adj(tI + (2t+1) N_S) 1`, of degrees `<= k` and `<= k+1`.
- The coefficients are obtained by exact interpolation, with one check point. Delta_{m-1}(t) comes from
  `fastregion.hurwitz_det_poly`.
- The critical set is the positive real roots of `c_m(t) = det L(t)` and of `Delta_{m-1}(t)`, isolated exactly. The
  proof of Lemma B(d) applies word for word with t as the parameter: on each open cell no root of p crosses iR.
- Each cell is decided by `ctln.hurwitz_stable` at an exact rational sample; if the minor vanishes at all three
  samples, the cell is `escalate`.
- If some `c_k` (k >= 1) vanishes identically in t, the verdict is `never`. By (b1) this is an identity over the
  whole curve, not a sampled boundary value. At n = 5 it occurred exactly for the 12 tournaments: on the curve
  e_2 = B(1+2t) - Nn(beta-1)(beta+1+2t), which vanishes identically iff B = Nn = 0 (§5).

Validation:
- n = 5, all 9,608 graphs: 9,513 never, 95 stable_somewhere, 0 escalate. This equals the float graph-level count
  (95), with 0 disagreements (`validate_curve_n5.py`). Each of the 95 has an exact legal witness checked by
  ctln.analyze_support at `beta = 2 + 1/t + 10^-6`.
- n = 7, all 124 graphs whose limit verdict is stable_somewhere: all return stable_somewhere, as (R) requires
  (`curve_on_limit_hits.py`).
- n = 7, 6,000 random positivity survivors, 1,500 from each of n7_pos_part{0..3}: **6,000 never**, 0.043-0.050 s
  per graph under full machine load (`run_curve_sample.py`, `curve_sample_p*.log`).
- n = 7, 10,000 further random survivors with edges >= 22: **10,000 never**, 0.035 s per graph
  (`curve_sample_dense.py`, `curve_dense_s101.log`).

**Exact fast path** (`curve1d.coeff_never`, `curve1d.curve_test_fast`). By (b1), if no t > 0 has all coefficients
c_1..c_m of p(., t) positive, the curve is never stable. That set is decided exactly from the positive roots of the
low-degree c_k (degree <= k+1).
- Speed and power: on 1,500 random dense n=7 survivors, the coefficient test alone decides `never` for 84.6% of
  graphs, at 0.0088 s per graph (`curve_coeff_prefilter.py`). It agrees with the full test on 300 of them.
- Equivalence: `curve_test_fast` and `curve_test` give identical verdicts on all 9,608 graphs with 5 nodes
  (`validate_fast_n5.py`).
- Expected average cost on dense n=7 lines: about 0.009 + 0.15 x 0.035, i.e. roughly 0.014 s per graph.

### 3.3 Strip decomposition: SFP itself is a union of 1-D problems (PROVED)

Let `(beta_l, beta_r)` be a connected component of the positivity set inside (2, inf); beta_r may be inf.
Positivity is t-free (Lemma A), so SFP in this strip means legal and stable there. On the ray of slope
`tau = t/beta`:
- the stable beta-set is `(0, beta_s(tau))` (Remark F);
- the legal beta-set is `(beta_leg(tau), inf)` with `beta_leg(tau) = 1 + sqrt(1 + 1/tau)`.

Hence SFP somewhere in the strip **iff**

- (c1) some curve point `K(2 + 1/t, t)` with `2 + 1/t in [beta_l, beta_r)` is stable, **or**
- (c2) `K(beta_l, t)` is stable for some legal `t > 1/(beta_l - 2)` (vacuous if beta_l = 2).

*Proof.* On ray tau the SFP set is `(max(beta_l, beta_leg), min(beta_r, beta_s))`.
- If `beta_leg >= beta_l`, it is nonempty iff `beta_leg < beta_r` and `beta_leg < beta_s`, which is (c1).
- If `beta_leg < beta_l`, it is nonempty iff `beta_l < beta_s`, which is (c2) at `t = tau beta_l`; by Remark F and
  openness, stability at beta_l transfers to `beta_l + eta`. ∎

(c2) at an irrational algebraic beta_l is not implemented. Survivors of (R) go to region2d, which is exact anyway.

---

## 4. Lemma AS: exact asymptotics near the limit boundary (PROVED; explains §2.2)

Let `chi(z) := det(zI + A)`, `psi(z) := 1^T adj(zI + A) 1`, and `s^(z) := psi/chi = 1^T (zI + A)^-1 1`.

- By (T2), `det(zI - (cJ - A)) = chi(z) - c psi(z)` identically.
- By the argument of C'(a), `psi(z) = m det(zI_(m-1) + A|perp)`.
- Let `h = gcd(chi, psi)`, `chi = h chi1`, `psi = h psi1`. The roots of h are eigenvalues of cJ - A for every c
  ("fixed"). The other eigenvalues ("moving") are the roots of `chi1 - c psi1`.

**Lemma AS.** Let z0 be a simple root of psi1. Then there are eps0 > 0 and an analytic z(e), |e| < eps0, with
`z(0) = z0`, such that for every `c > 1/eps0`, `z(1/c)` is an eigenvalue of cJ - A, and

    z(1/c) = z0 + 1 / (c s^'(z0)) + O(c^-2).

*Proof.* Put e = 1/c. The equation `chi1 - c psi1 = 0` becomes `F(z, e) := psi1(z) - e chi1(z) = 0`. We have
`F(z0, 0) = 0`, and `F_z(z0, 0) = psi1'(z0) != 0` because z0 is simple. The implicit function theorem gives an
analytic z(e) with `z'(0) = chi1(z0)/psi1'(z0)`. Here `chi1(z0) != 0` since gcd(chi1, psi1) = 1. Finally
`s^'(z0) = psi1'(z0)/chi1(z0)`, because `psi1(z0) = 0`. ∎

Note that z0, s^ and R below **do not depend on beta**. Beta enters only through gamma and through c = 1 + t/beta.

**Theorem AS (local form of M).** Assume:
- `nu_inf := -alpha` is in (1/2, 1), so `beta_c = 1/(1 - nu_inf) > 2`;
- every root z0 of psi with real part nu_inf has multiplicity one in psi and is not a root of chi. Then z0 is not
  a root of h, it is a simple root of psi1, and it is a moving eigenvalue.

Let `R* := min Re(1/s^'(z0))` over those roots. Then `nu_A(c) = nu_inf + R*/c + O(c^-2)`.
- **(a)** If `R* > 0`, there is eta > 0 such that for every beta in (beta_c, beta_c + eta), T(beta) is nonempty and
  bounded above. If moreover `a(beta_c) > 0`, these points are SFP. Both forms of M fail.
- **(b)** If `R* < 0`, `nu_A` is increasing for all large c, so `T(beta) ∩ (T0, inf)` is upward closed for beta
  near beta_c. M holds near `(beta_c, inf)`.

The boundary of the stable region near `(beta_c, inf)` is `beta = beta_c + beta_c^3 R*/t + O(t^-2)`.

*Proof.*
1. For large c the other bounded eigenvalues stay near roots with real part > nu_inf, and the large one goes to
   +inf. So the minimum real part is attained on the branches z(1/c) of Lemma AS. Conjugate branches have equal
   real parts. Hence `nu_A(c) = Re z(1/c) = nu_inf + R*/c + O(c^-2)`.
2. (a) Pick c2 so large that `nu_A(c) >= nu_inf + R*/(2c)` for c >= c2, and that c2 is legal for every
   gamma <= (1 + nu_inf)/2. This is possible because `gamma^2/(2gamma-1)` is bounded for gamma in
   `[nu_inf, (1 + nu_inf)/2]`.
3. For gamma in `(nu_inf, nu_inf + R*/(2 c2))`: stable at c2 by (1), and unstable for all large c because
   `nu_A(c) -> nu_inf < gamma`. The positivity set is open, so it contains a neighbourhood of beta_c when
   `a(beta_c) > 0`.
4. (b) `d/dc Re z(1/c) = -(1/c^2)(R* + O(1/c)) > 0` for large c.
5. The boundary formula follows from `nu_A(c) = gamma`, using `gamma_c - gamma ≈ (beta - beta_c)/beta_c^2` and
   `c ≈ t/beta_c`. ∎

CERTIFIED check against the repo engine for `&FK@dyUfXp?`:
- `R* = 8.1996e-4`, so the prediction is `beta_c^3 R* = 0.0077260`;
- the exact endpoints give `t (beta_right(t) - beta_c) = 0.0077476` at t = 10^3 and `0.0077281` at t = 10^4.

The good-sign graph `&FMC`j\BTt?` (R* < 0) has right ends below beta_c, increasing to it (t = 100 ... 10^5); M
holds there. That is the n=6 situation (Geneson class: R* < 0).

HEURISTIC table (`n7_first_order_table.json`, not a certificate) for the 124 n=7 limit-stable graphs:

| class | count |
|---|---|
| positive at beta_c, R* > 0: stable region extends past beta_c, T(beta) bounded there | 39 |
| positive at beta_c, R* < 0 | 54 |
| not positive at beta_c | 30 |
| K7 (fixed critical eigenvalue) | 1 |

Stage 6b's 2-D maps should therefore show stable cells **bounded above in t** for about 39 graphs. That is a
checkable prediction, and it differs from n = 6.

---

## 5. Smaller proved facts

**e_2 closed form.** Count unordered pairs {i, j}: B bidirected, Nn non-adjacent, O one-way. The 2x2 principal
minors of K are:

| pair type | minor |
|---|---|
| bidirected | 1 + 2t |
| non-adjacent | (1-beta)(1+beta+2t) |
| one-way | 1 - t(beta-2), negative exactly when legal |

So

    e_2(K) = B(1+2t) - Nn(beta-1)(beta+1+2t) + O(1 - t(beta-2)) = t Phi(beta) + B + O - Nn(beta^2 - 1),
    Phi(beta) = m(m-1) - beta e,   e = 2 Nn + O (ordered non-edges).

Exact check: 300 random cases, `check_e2.py`.
- Stability needs `e_2 > 0` (b1). At a legal point this forces `B > Nn`, i.e. more than half of the ordered pairs
  are edges. That is exactly the trace condition C'(f) that limit stability at beta -> 2+ requires.
- If `Nn(beta^2 - 1) >= B + O`, it forces `Phi(beta) > 0`, which is C'(f) at the same beta.
- This is a proved finite-t analogue of C'(f), the open item in notes/lemmas.md §9. On the n=5 data all 95
  legal-stable graphs satisfy B > Nn.
- Cheap exact prefilter: at n = 7, B > Nn is `edges >= 22`, which keeps 27-67% of the positivity-survivor lines,
  depending on the part file.

**Regular graphs.**
- If A has constant in-degree d, then `1^T A = d 1^T`, so `s^(z) = m/(z + d)`. The only moving eigenvalue is
  `cm - d`, and all others are fixed. Hence `nu_A` is nondecreasing and M_s holds.
- Constant out-degree: the same with `A 1 = d 1`.

**Real crossings** (Lemma A.2): with s(beta) > 0, which holds on the positivity set, D = d + ts increases in t, so a
real eigenvalue of K crosses 0 only upward. Every M failure under positivity is an imaginary-pair crossing, as in §2.2.

**On the suggested tools.**
- Crossing direction, `lambda' = (w^T 1)(1^T v)/(w^T v)`: the pointwise claim "at a legal point with a(beta) > 0,
  every imaginary-axis crossing moves into Re > 0 as t grows" is false. HEURISTIC evidence: 18 downward legal
  crossings at beta > 2 among 3,000 n=7 positivity survivors x 50 beta values (`crossings.py`, `scan_cross_d6.py`).
  Theorem AS identifies the obstruction as the sign of `Re(1/s^'(z0))`; neither positivity nor legality fixes it.
- `M(eps) = J + eps B(beta)` is `I - W`; it is the same family as K at fixed beta.
- Rantzer convex directions and Bialas segment tests: §2.1 and §2.2 show that the rank-one direction J (same vector
  1 on both sides) is not a stability-monotone direction, even under positivity and domination-freeness. Such
  theory can certify a given segment but cannot yield M.

---

## 6. Consequences for the n=7 finite-eps map

1. **Do not use Lemma M.** It is false. Its graph-level repair with domination-freeness is also false at n = 7. The
   only surviving form, M_graph^+ (§2.4), is unproved.
2. **Sound replacement pipeline. Every step is exact and proved.**
   - (i) Positivity survivors: the existing exact filter, `n7_pos_all_PU.txt.gz`, has 9,632,133 lines
     (8,845,537 P and 786,596 U, with the U lines being resolved).
   - (ii) `edges >= 22`, which is exact by §5 (legal stability forces B > Nn). This leaves 4,045,243 lines.
   - (iii) The Lemma-R curve test. A `never` verdict is a complete certificate: no legal (eps, delta) makes the full
     support stable, so there is no SFP. It stores the critical roots and the cell samples.
   - (iv) region2d only on the survivors of (iii). By (R) these are exactly the graphs that are legal-stable somewhere.
     They include the 124 limit-stable graphs and graphs like the two in §2.4.
3. **Cost.** Measured, pure Python, on this 4-core machine under full load: 0.035 s per graph on dense lines
   (10,000 random dense lines, all `never`), and 0.045 s on unfiltered lines.
   - 4.05M x 0.035 s is about 40 CPU-hours with `curve_test`;
   - about 16 CPU-hours with `curve_test_fast` (§3.2), i.e. well under an hour on a rented 64-core machine and
     about 4-5 hours on this 4-core box;
   - plus about 10 s per survivor of (iii) for region2d;
   - the full 2-D map on every survivor (about 9.6M x 10 s, roughly 1,100 CPU-days) is not needed.
4. This settles M_graph^+ at n = 7 as a by-product. It holds iff every region2d SFP hit has alpha <= -1/2, which is one
   exact Hurwitz test per hit. A hit with alpha > -1/2 would be a counterexample graph **invisible to the limit
   sweep**. That is the class the finite-eps map exists to find. STOP and hand-verify if one appears.

---

## 7. Gaps (named)

- **G1.** M_graph^+ is OPEN in general and at n = 7. No proof idea survived. §2 shows that any proof must couple
  positivity and stability at the same beta: stability plus legality, plus domination-freeness, plus positivity at
  another beta, each still fails.
- **G2.** Closed: the n = 6 exhaustive exact domination-free run finished (§2.4) with 28 hits and 0 escalations.
  The cost is pure Python, about 40 min per half on this machine.
- **G3.** In the curve test a cell is `escalate` if a minor vanishes at three samples. None occurred on all 9,608
  n=5 graphs, the 124 n=7 hits, the 16,000 random n=7 survivors, or the alpha-window search. (c2) of §3.3 at an
  irrational beta_l is not implemented; region2d covers it.
- **G4.** All heuristic scans (`scan_*.py`, `first_order*.py`, `R_survey.py`, `n7_first_order_table.py`,
  `alpha_hist.py`, the float stage of `targeted_near_boundary.py`) are searches only; their counts are not
  certificates. The disproofs rest only on `certify_all.py`, which re-runs every claim C1-C5 exactly and passes.
