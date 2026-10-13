# Lemma J: when do two graphs have the same full-support data D(A)?

Prover notes, 2026-10-13. Conventions exactly CLAUDE.md §2 and notes/lemma_I_spectral.md §0:
`A[i][j] = 1` iff `i -> j`; `N := J - I - A^T` (`N[i][j] = 1` iff `i != j` and `j -/-> i`);
`C := I + beta N`, `K := C + tJ`, `a := adj(C) 1`, `d := det C`, `phi(mu) := det(mu I + C)`,
`omega(mu) := 1^T adj(mu I + C) 1 = det(mu I + C + J) - det(mu I + C)`.
`D(A) := (d, phi, omega, multiset{a_i})`, all polynomials in `beta` (and `mu`), exactly as computed by
`fastregion.beta_data(A, range(n))`. `chi_X(z) := det(zI - X)`. For a matrix `X` and vectors `v_1..v_r`,
`Kr(X; v_1..v_r)` is the smallest `X`-invariant subspace containing them (= span of all `X^k v_j`).
"Relabelling" = simultaneous permutation of rows and columns (`graphs.permute`).

Labels: **PROVED** (proof here; linear algebra only), **CHECKED** (exact integer/rational computation,
script named; scripts are in /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/ and must be copied into the repo to be
re-runnable from it -- I was told to modify no file but this one), **OBSERVED** (data, no proof).
No floating point anywhere: numpy is used only with int64 arrays on integers below 10^4.

---

## 0. Summary

| # | statement | status |
|---|---|---|
| 1 | G2 is G1 with ONE edge moved: relabelled, `G2 = G1 - (1->4) + (2->4)`, where 1, 2 are bidirected in-twins | CHECKED (s1, s2) |
| J.0 | common invariant subspace containing 1 + cospectral quotients ⇒ same D, same `a` as a vector, same full-support fixed point and Jacobian charpoly for all (eps, delta) | PROVED |
| J.1 | rank-one moves `N' = N + y x^T` with `x ⊥ Kr(N; 1, y)` preserve D; in-twin switching is the case `x = e_u - e_v` | PROVED |
| J.2 | **`D(A) = D(B)` iff B can be relabelled so that A and B are cospectral and every vertex has the same number of length-k incoming walks in both, for every k ≤ n-1** | PROVED; the two partitions coincide on all 34,071 survivors (s9) |
| J.3 | necessary conditions: in-degree sequence, #edges, #2-cycles, #3-cycles, `chi_A`, `chi_(A+sJ)` for all s, `chi_(A|perp)`, total walk counts. Out-degree sequence is NOT invariant | PROVED / CHECKED |
| J.4 | **bidirected in-twin switching preserves every fixed point of every support (same vector, same Jacobian charpoly), for every eps, delta > 0 and theta; hence FP(G), the stable supports and the target-free cliques. G1 and G2 are related this way** | PROVED; CHECKED (s10, s10b, s11, s17) |
| J.4' | two non-adjacent in-twins inside a support give the exact Jacobian eigenvalue +delta, so that support is never stable | PROVED |
| J.5 | general D-equivalence (even non-adjacent in-twin switching) does NOT preserve FP(G) or the stable supports | CHECKED with exact counterexamples (s11, s14-s16) |
| gen | generating moves on the 28,350 strict survivors: in-twin switching leaves 741 of the 17,307 D-classes split; rank-one moves leave 95 split (in-regular cospectral classes etc.); no finite list of local moves is claimed | CHECKED / OBSERVED |
| anom | `certificates/n6_survivors.idx` mixes two positivity filters (switch at idx ≈ 318,000 in every shard) | CHECKED (s5 + shard reasons) |

---

## 1. The G1 / G2 relation (CHECKED: s1_relabel.py, s2_g1g2.py, s17_bidir_n6_and_G1.py)

Relabel G2 by p = (4,3,5,2,1,0) (`graphs.permute`: new node i = old node p[i]):
```
G1  = 011011 / 101011 / 110101 / 111000 / 011000 / 000110
G2p = 011011 / 101001 / 110111 / 111000 / 011000 / 000110
```
They differ in exactly two entries, which is the minimum over all 720 relabellings. **The edge 1->4 of G1
is replaced by 2->4.** In G1 nodes 1 and 2 are bidirected (1<->2) and have the same in-neighbours among
the others (`in(1) \ {2} = in(2) \ {1} = {0,3,4}`): they are *bidirected in-twins* (§3). Hence
`N(G2p) = N(G1) + e_4 (e_1 - e_2)^T` with `(e_1 - e_2)^T N(G1) = 0` and `(e_1 - e_2)^T 1 = 0`.
`beta_data(G1)` and `beta_data(G2p)` agree entry for entry: `a` agrees as a VECTOR in this labelling.
Only the out-degrees of the twins change (1: 4 -> 3, 2: 4 -> 5). No in-degree changes: node 4 trades
in-neighbour 1 for in-neighbour 2. The class {G1, G2} is exactly the in-twin switching orbit of G1:
the out-neighbourhoods of 1 and 2 differ only at w ∈ {3, 4}; switching at 4 gives G2p, switching at 3 gives
its image under 1 <-> 2, and switching at both is the relabelling 1 <-> 2.

This is the same relation that ANOMALIES.md (entry "Candidate idx 1424629 verified; the equivalence operation
is a 'twin move'") describes in another labelling: twins 3 and 5 there, 1 and 2 here. That entry left open whether
twin moves generate all the classes. They do not (§7). The exact characterisation is J.2.

Other candidate operations from the brief, checked:
- **Transposition.** It is not this, and it does not preserve D here: `D(G1^T) != D(G1)`, because the
  in-degree multiset of G1^T is G1's out-degree multiset [2,2,3,4,4,4], not [2,3,3,3,4,4] (J.3).
- **`P N Q` with permutations.** Also not this: `N(G2p) - N(G1)` is a single row-local swap.
- **Similarity fixing 1.** Yes, here: `N(G2p) = S^{-1} N(G1) S` with
  `S = I + z (e_1 - e_2)^T`, `z = (-1,1,1,0,0,0)^T` (so `N(G1) z = e_4`), and `S 1 = 1`.
  S is an integer matrix with det S = 1. The general mechanism is J.0 (a similarity is
  not always available; see the remark after J.0).

---

## 2. The block lemma (PROVED)

**Lemma J.0.** Let `N, N'` be real `n x n` matrices and `V` a subspace with `1 in V`, `N V ⊆ V`, and
`N v = N' v` for all `v in V` (so `N' V ⊆ V` too). Let `Z, Z'` be the maps induced on `R^n / V`.
If `chi_Z = chi_Z'` then:
  (a) `det(cI + N + sJ) = det(cI + N' + sJ)` as polynomials in `(c, s)`;
  (b) `d' = d`, `phi' = phi`, `omega' = omega` as polynomials in `(beta, mu)`;
  (c) `adj(I + beta N' + tJ) 1 = adj(I + beta N + tJ) 1` as polynomials in `(beta, t)`. In particular
      `a' = a` as vectors (t = 0). Also `det(mu I + I + beta N' + tJ) = det(mu I + I + beta N + tJ)`.

*Proof.* In a basis whose first `r = dim V` vectors span `V`:
`N = [[X, Y], [0, Z]]` and `N' = [[X, Y'], [0, Z']]`. The two matrices have the same `X` because they agree
on `V`. Because `1 in V`, `J = 1 1^T` maps `R^n` into `V`, so in this basis
`J = [[J_1, J_2], [0, 0]]`. Hence `cI + N + sJ = [[cI + X + sJ_1, *], [0, cI + Z]]`, and the same holds with
primes, so `det(cI + N + sJ) = det(cI + X + sJ_1) det(cI + Z)`. Here `det(cI + Z) = (-1)^(n-r) chi_Z(-c)`, which is
the same for `Z'` because `chi_Z = chi_Z'`; that proves (a).
(b) `d = det(I + beta N)`, `phi = det((mu+1) I + beta N)`, and `phi + omega = det((mu+1) I + beta N + J)`.
For `beta != 0` each equals `beta^n det((c/beta) I + N + sJ)` with `s ∈ {0, 1/beta}`, so (a) makes them equal
for `N` and `N'`. Polynomial identities that hold for all `beta != 0` hold identically.
(c) Put `M := I + beta N + tJ` and `M' := I + beta N' + tJ`. By (a) `det M = det M'`, and as a polynomial
in `(beta, t)` it takes the value 1 at (0,0), so it is nonzero on a dense open set. Fix `(beta, t)` there.
`V` is `M`-invariant and `M|_V = M'|_V`. `M` is injective, so `M|_V` is a bijection of `V`, and
`u := (M|_V)^{-1} 1 ∈ V` is the unique solution of `M u = 1`. It also solves `M' u = 1`. So
`adj(M) 1 = det(M) u = det(M') u = adj(M') 1`. Both sides are polynomials that agree on a dense set, so they
are equal. The last claim is (a) with `c = mu + 1`, rescaled as in (b). ∎

**Consequence (dynamics).** `I - W = eps K`, so for every `eps, delta, theta` the full-support solution
`x* = (theta/eps) K^{-1} 1` is the same vector for `N` and `N'`, and `-I + W = -eps K` has the same
characteristic polynomial. The full-support verdict therefore coincides everywhere: positivity, stability,
the whole 2-D region and every critical value. This goes beyond the eps -> 0 limit.

**Remark (when it is a similarity).** Suppose moreover `Z' = T^{-1} Z T` for some invertible T, and
`spec X ∩ spec Z = ∅`. Then the Sylvester equation `X R - R Z' = Y' - Y T` has a solution R, and
`S := [[I, R], [0, T]]` satisfies `N S = S N'`, i.e. `N' = S^{-1} N S`, with `S|_V = id` (so `S 1 = 1`). In general (shared eigenvalues, or quotients that are cospectral but not similar)
J.0 does not give a similarity. Only the cospectral statement is claimed.

---

## 3. Rank-one moves and in-twin switching (PROVED)

**Lemma J.1.** Let `x, y ∈ R^n` with `x ⊥ N^k 1` and `x ⊥ N^k y` for all `k >= 0`, i.e.
`x ⊥ Kr(N; 1, y)`. Then `N' := N + y x^T` satisfies J.0 with `V = Kr(N; 1, y)` and `Z' = Z`. Hence
`D(N') = D(N)` with `a' = a` as vectors, and everything in the Consequence above holds. The same `x, y`
satisfy the hypothesis for `N'` (`Kr(N'; 1, y) = V`), so the move is reversible.

*Proof.* For `v ∈ V`, `N' v = N v + y (x^T v) = N v`. So `N, N'` agree on `V`, and `V` is invariant under
both. `N' - N = y x^T` has range in `span(y) ⊆ V`, so it induces 0 on `R^n / V`, giving `Z' = Z`. ∎

**Converse for rank one (PROVED).** Suppose `D(N + y x^T) = D(N)` with `a` equal as vectors, and `y != 0`.
Then `x ⊥ Kr(N; 1, y)`. Proof: `d' = d` and `a' = a` give `C'^{-1} 1 = C^{-1} 1`. Sherman-Morrison gives
`C'^{-1} 1 - C^{-1} 1 = -beta (x^T C^{-1} 1) C^{-1} y / (1 + beta x^T C^{-1} y)`. Since `C^{-1} y != 0`, this
forces `x^T C^{-1} 1 = sum_k (-beta)^k x^T N^k 1 ≡ 0`. Also `phi' = phi` gives `chi_(N') = chi_N`, and
`chi_(N + y x^T)(z) = chi_N(z) - x^T adj(zI - N) y`, so `x^T (zI - N)^{-1} y ≡ 0`. Expanding at
`z = ∞` gives `x^T N^k y = 0` for all k.

**In-twins.** Call `u != v` *in-twins* if `A[k][u] = A[k][v]` for every `k ∉ {u, v}` and `A[u][v] = A[v][u]`.
"Bidirected" means `A[u][v] = 1`; "non-adjacent" means `A[u][v] = 0`. Rows `u, v` of `N` then agree off
columns `{u, v}`, and `N[u][v] = N[v][u] =: c` (`c = 0` bidirected, `c = 1` non-adjacent). So `x := e_u - e_v`
is a left eigenvector, `x^T N = -c x^T`, with `x^T 1 = 0`. Then `x ⊥ N^k 1` always, and `x ⊥ N^k y` whenever
`y_u = y_v`, because `x^T N^k y = (-c)^k x^T y`.

**In-twin switching.** Take `w ∉ {u, v}` with exactly one of `u -> w`, `v -> w`, and replace that edge by
the other. In `N` this is `N' = N ± e_w (e_u - e_v)^T`, a J.1 move (`y = ±e_w`). It preserves D with `a`
equal as a vector. `u, v` stay in-twins, since their in-neighbourhoods are untouched, so switches at
different `w` commute and combine. In words: **for in-twins u, v, the out-edges of u and v to third
vertices may be redistributed arbitrarily, provided each third vertex w keeps the same count
`A[u][w] + A[v][w]`.** Every in-degree is preserved; only the out-degrees of u, v change.
Switching at every w at once is the relabelling `u <-> v`.

---

## 4. The characterisation (PROVED; partition identity CHECKED)

**Krylov agreement lemma.** If `X^k v = Y^k v` for all `k`, then `X` and `Y` agree on `Kr(X; v) = Kr(Y; v)`.
Proof: `X (X^k v) = X^{k+1} v = Y^{k+1} v = Y (Y^k v) = Y (X^k v)`. ∎

**Theorem J.2.** For digraphs A, B on n vertices, `D(A) = D(B)` iff there is a relabelling of B such that
  (i)  `(A^T)^k 1 = (B^T)^k 1` for `k = 0, ..., n-1`. Equivalently for all k ≥ 0, given (ii). Entry i of
       `(A^T)^k 1` is the number of directed walks of length k in A that END at i.
  (ii) `chi_A = chi_B` (A and B are cospectral).
Equivalently: `chi_A = chi_B`, and the multisets over vertices of the in-walk sequences
`(((A^T)^k 1)_i)_{k=0..n-1}` agree. When this holds, the relabelling can be chosen so that `a` agrees as a vector.

*Proof.* Write `N_A = J - I - A^T`. If `V ∋ 1` then `J V ⊆ V`, so `V` is `N_A`-invariant iff it is
`A^T`-invariant. On such a `V`, `N_A` and `N_B` agree iff `A^T` and `B^T` do.
(⇒) Relabel B so that `a_B = a_A` as vectors (the multisets agree; D is relabelling-invariant).
`d_A = d_B` and `d(0) = 1`, so as formal power series `a/d = C^{-1} 1 = sum_k (-beta)^k N^k 1`. Hence
`N_A^k 1 = N_B^k 1` for all k. By the Krylov lemma, `N_A` and `N_B` agree on `V := Kr(N_A; 1) ∋ 1`. So do
`A^T` and `B^T`, and `V = Kr(A^T; 1)`, which gives (i). For (ii), work in a basis adapted to `V`:
`N = [[X, *], [0, Z]]` and `A^T = J - I - N = [[J_1 - I - X, *], [0, -I - Z]]`, with `X` common to A and B.
`phi_A = phi_B` gives `chi_(N_A) = chi_(N_B)`, so `chi_(Z_A) = chi_(Z_B)`, so `chi_(A^T) = chi_(B^T)`.
(⇐) By (i) and the Krylov lemma, `A^T` and `B^T` agree on `V := Kr(A^T; 1)`, hence so do `N_A` and `N_B`.
In the adapted basis `A^T = [[P, *], [0, Q]]` with `P` common, so (ii) gives `chi_(Q_A) = chi_(Q_B)`. The quotient
of `N` is `-I - Q`, so the quotients are cospectral and J.0 applies.
(k ≤ n-1 suffices) Cayley-Hamilton with the common `chi` gives both sequences `v_k = (A^T)^k 1` and
`w_k = (B^T)^k 1` the same order-n linear recurrence, so they agree for all k once they agree for k < n. ∎

**CHECKED (s9_invariant_check.py).** Partitioning the 34,071 survivors by
`(chi_A, sorted per-vertex in-walk sequences)` gives exactly the 21,800 classes of
`certificates/n6_survivor_equivalence_classes.json` (identical partitions).
**CHECKED (s18_weaker_invariants.py).** The walk sequences cannot be shortened on this data.
Truncating at length 1, 2, 3, 4 merges 3,136 / 1,271 / 385 / 54 D-classes. Length 5 = n-1 merges none.
`(chi_A, chi_(A+J), in-degree multiset)`, i.e. generalised spectrum plus in-degrees, merges 1,229.

**Reading.** The full-support CTLN data sees the graph only through its spectrum and its in-walk
counts. The equilibrium `x* ∝ K^{-1} 1` is a generating function of the in-walk vectors of the
complement graph: `C^{-1} 1 = sum (-beta)^k N^k 1`, and `(N^k 1)_i` counts length-k walks of the
complement digraph ending at i. Out-walk data (the out-degrees, for one) never enters.

**In-regular graphs.** If `A^T 1 = r 1` (every in-degree r), then `V = span(1)`, (i) holds automatically,
and D-equivalence is plain cospectrality among r-in-regular graphs; cospectrality also forces equal r,
since r is the spectral radius. The largest class (259 graphs) is exactly of this kind: every member has
in-degree 3, and they show 28 different out-degree sequences (s13).

---

## 5. Necessary conditions and non-invariants (J.3)

PROVED from J.2. If `D(A) = D(B)` then A and B have equal:
- in-degree multiset (k = 1), and more strongly the multiset of per-vertex in-walk sequences;
- number of edges `1^T A 1` (sum of in-degrees) and total walk counts `1^T A^k 1` for all k;
- characteristic polynomial; hence `tr A^k` for all k, so equal numbers of 2-cycles (bidirected pairs,
  `tr A^2 / 2`) and of directed 3-cycles (`tr A^3 / 3`);
- generalised characteristic polynomial `det(zI - A - sJ)` for all s (J.0(a), rewritten via
  `N + sJ = (1+s)J - I - A^T` and transposition); hence `chi_(A|perp)` (lemma_I_spectral.md I.1b:
  `m chi_(X|perp) = chi_X - chi_(X+J)`), the limit-stability data.

CHECKED (s14). Each of these is constant on all 4,653 non-singleton classes of the 34,071 survivors.
**Not invariant:** the out-degree multiset is constant on only 214 of the 4,653 (G1/G2 is an example).
OBSERVED, not proved: domination-freeness was constant on all 4,653 classes. This is not explained.
"Admits a positive full-support solution for some parameters" is a D-invariant, and domination rules
that out (Lemma D), but a domination-free graph can also lack such a solution, so this is not a proof.

---

## 6. Is the equivalence CTLN-natural? (Task 3)

### 6.1 Bidirected in-twin switching: yes, completely (PROVED)

**Theorem J.4.** Let u, v be bidirected in-twins of G (`u <-> v`, `A[k][u] = A[k][v]` for all `k ∉ {u,v}`), let
`w ∉ {u, v}` with `u -> w`, `v -/-> w`, and let G' be G with `u -> w` replaced by `v -> w`. Fix any
`eps > 0`, `delta > 0`, `theta > 0`; legality is not needed. Then for EVERY support `sigma ⊆ [n]`:
sigma satisfies conditions (1) and (2) of CLAUDE.md §2 in G iff it does in G'. When it does, the fixed point is
the same vector and `-I + W_sigma` has the same characteristic polynomial in G and G'. Also
`det(I - W_sigma)` agrees whenever sigma does not contain exactly one of u, v. Hence FP(G) = FP(G'),
the stable supports coincide, and so does every three-valued verdict of ctln.analyze_support on every
fixed-point support. The target-free cliques of G and G' coincide too.

*Proof.* W and W' differ only in row w: `W[w][u] = -1+eps`, `W[w][v] = -1-delta`, and in W' these
are swapped. In-twins give `W[u][k] = W[v][k]` for `k ∉ {u,v}`. Bidirected gives `W[u][v] = W[v][u] = -1+eps`.
Both facts also hold in G'.
- *Exactly one of u, v in sigma*, say `u ∈ sigma`, `v ∉ sigma`. If x is a fixed point with support sigma,
  then `x_u = theta + sum_(j ∈ sigma, j != u) W[u][j] x_j > 0`. The drive at v is
  `theta + sum_(j ∈ sigma) W[v][j] x_j = theta + sum_(j != u) W[u][j] x_j + W[v][u] x_u = x_u + (-1+eps) x_u = eps x_u > 0`,
  so condition (2) fails at v. The same holds in G'. So sigma is in neither FP(G) nor FP(G').
- *Neither in sigma.* `W_sigma = W'_sigma`, and every off-support drive uses only columns in sigma, which
  are unchanged.
- *Both in sigma, w ∈ sigma.* Let `M := I - W_sigma` and `xi := e_u - e_v`. Then
  `M' = M + (eps+delta) e_w xi^T`, and `xi^T M = eps xi^T` (twin rows; `W[u][u] = 0`, `W[v][u] = -1+eps`).
  Also `xi^T 1 = 0` and `xi^T e_w = 0`. Lemma J.1 in its M-form (same proof: `V = Kr(M; 1, e_w)` and
  `xi ⊥ V`, since `xi^T M^k = eps^k xi^T`) gives `det(mu I + M') = det(mu I + M)` for all mu, and
  `M'^{-1} 1 = M^{-1} 1` when M is invertible. So x, positivity and the Jacobian charpoly are identical.
  Off-support drives at `k ∉ sigma` use row k, and `k != w`, so they are unchanged.
- *Both in sigma, w ∉ sigma.* `W_sigma = W'_sigma`, so x is the same, and drives at `k != w` are unchanged.
  At w the drive changes by `(eps+delta)(x_v - x_u)`. From `M x = theta 1` and `xi^T M = eps xi^T` we get
  `eps (x_u - x_v) = theta xi^T 1 = 0`, so `x_u = x_v` and the drive is unchanged.
For target-free cliques: a clique containing u but not v has v as a target (every member sends to u,
hence to v by twinship, and u -> v), so target-free cliques contain both twins or neither. A clique
containing both has `A[u][w] + A[v][w] = 2` for every member w, so the switch touches no edge inside it,
and the switched w (count 1) is not a target before or after. ∎

**Remark J.4' (exact twin eigenvalue; PROVED).** Let u, v in sigma be in-twins *relative to sigma*
(`A[k][u] = A[k][v]` for `k ∈ sigma \ {u,v}`, and `A[u][v] = A[v][u]`). Then `xi = e_u - e_v` is a left eigenvector
of `-I + W_sigma`. The eigenvalue is `-eps` if u, v are bidirected and `+delta` if they are non-adjacent:
`xi^T W_sigma = (1 - eps) xi^T` resp. `(1 + delta) xi^T`, by the computation in the proof above with
`W[v][u] = -1 - delta`. Hence **a support that contains two non-adjacent in-twins (relative to itself) is never
stable, at any eps, delta > 0.** This is a cheap exact prune for Stage 4. Whether it is already in the literature
is unchecked. The bidirected case explains the exact factor `(101 lam + 1)` seen in the slices of G1/G2
(eigenvalue `-eps = -1/101`).

**Corollary.** G1 and G2 (hit 3, idx 1424629) are related by a bidirected in-twin switch (§1). So they have
the same fixed points with the same stability for every support at every (eps, delta, theta). Hit 3 is not
an independent counterexample: it is Geneson's counterexample with one inert edge moved. For Stage 5
counting, counterexamples should be reported modulo bidirected in-twin switching, and better still
modulo D-equivalence for the full-support question.

CHECKED, exact:
- s10_fp_g1g2.py / s10b: G1 vs G2p, all 63 supports, at (eps, delta) = (1/2,3), (1/3,1), (1/4,1/2), (1/10,1/8),
  (1/101, 27/2626), (1/5,1), (1/20,1/10). Identical FP sets, stable sets (full support stable only at
  the Geneson point eps = 1/101), x vectors and charpolys on every fixed-point support. Every differing
  record is a support with exactly one of {1, 2}, and it is a non-fixed point in both graphs. No undecided.
- s11 (n=5, 150 bidirected moves, 2 parameters) and s17 (n=6, 120 strict survivors, 2 parameters):
  0 violations.

### 6.2 Everything else: no (CHECKED, exact counterexamples)

- **Non-adjacent in-twins (n=5, s11).** A = 00000/00001/00010/10000/10000 (edges 1->4, 2->3, 3->0, 4->0).
  Nodes 1 and 2 have no in-neighbours and are not adjacent, so they are non-adjacent in-twins.
  B = A with 1->4 replaced by 2->4. D(A) = D(B). At eps = 1/4, delta = 1/2 the stable supports are {0} for A
  and {0}, {1} for B: in B, {1} is a sink singleton, i.e. a target-free clique. The proof of J.4 breaks
  exactly here: for a non-adjacent twin pair the drive at v is `-delta x_u < 0`, so supports containing one
  twin can be fixed points.
- **Other J.1 moves (n=5, s11).** Non-twin pairs `x = e_u - e_v` changed the FP set in 59 of 150 tested
  moves; support-4 x in 116 of 150. Witnesses are printed by s11.
- **Inside the n=6 survivor classes.** All examples are strict survivors with identical full-support data.
  - idx 353 vs 302942 (s15): |FP| = 5 vs 7 at eps = 1/4, delta = 1/2.
  - idx 260864 vs 340100 (s16), the same D-class (with 355 and 538128): stable supports {0}, {4,5} vs
    {2,4}, at both (1/4, 1/2) and (1/10, 1/8).
  - Combinatorially, 1,492 of the 4,029 non-singleton strict D-classes have members with different
    target-free-clique size profiles.

**Conclusion.** D-equivalence is a statement about the full support only. Bidirected in-twin switching is
the sub-equivalence that is fully CTLN-natural (PROVED). Every other move tested changes FP(G), and
general D-classes mix graphs with different stable supports (CHECKED).

**For Stage 4 (remark, follows from Lemma A of notes/lemmas.md).** For a support sigma, the verdict at
every `(beta, t)` depends only on `D(G|sigma)` (with `a` as a vector) and the multiset of off-support
polynomials `o_k = d - beta sum_j N[k][j] a_j`. Any equivalence used to save work at Stage 4 must preserve
that pair. J.4 does, for all supports at once.

---

## 7. Which moves generate the classes? (CHECKED / OBSERVED)

The right universe is the 28,350 **strict** survivors (all `a_i > 0` at some legal beta in the limit). That
set is closed under every D-preserving move, because `a` is preserved. Union-find over canonical forms
(exact; s6, s8, s13):

| moves used | components | D-classes split |
|---|---|---|
| none | 28,350 | -- |
| bidirected in-twin switching only (CTLN-natural, J.4) | 23,055 | -- |
| all in-twin switching (J.1, `x = e_u - e_v` twin) | 19,194 | 741 |
| all rank-one J.1 moves with `x, y ∈ {0,±1}^6` | 17,459 | 95 |
| D-classes (= J.2 classes) | 17,307 | 0 |

- 13,278 of the 17,307 classes are singletons. Rank-one move types used, counted as directed moves by
  (|supp x|, |supp y|): (2,1) 74,986; (2,2) 49,475; (2,3) 16,305; (2,4) 2,382; (4,1) 10,882; (4,2) 842.
  No move ever left the strict set or changed D (0 of 50,855 twin moves, 0 of 154,872 rank-one moves).
  This is a consistency check of J.1 against `beta_data`.
- In each of the 95 residual classes the members can be relabelled to match the in-walk vectors (consistent
  with J.2), but the smallest rank of `N_B - N_A` over all such relabellings is 2 (91 classes) or 3
  (4 classes) (s12). The Krylov dimension `dim Kr(A^T; 1)` is 1 in 55 of them (in-regular: plain
  cospectrality), 2 in 25 and 3 in 15. These are "global" cospectral coincidences. I do not claim any
  finite list of local moves generates D-equivalence. The 259-member class is connected by rank-one moves
  (7 in-twin components).
- On the mixed 34,071-graph file the same computation gives 23,941 in-twin components, 878 split classes,
  and 3,050 moves that land outside the file. All 3,050 are explained by the anomaly below.

---

## 8. Anomaly for the orchestrator (not logged in ANOMALIES.md -- I was told to edit no other file)

`certificates/n6_survivors.idx` (34,071 graphs) mixes two positivity filters. In each limit shard
(`region_n6_limit_shard{0,1,2}of3.jsonl`), every record before about idx 318,000 that has no cert carries
the reason "no legal beta with a common sign for all a_i". After that point the reason is
"no legal beta with all a_i > 0 and all o_k < 0". The switch is at idx 318,984/318,987,
318,085/318,091 and 318,236/318,239 in the three shards. So the first part of the run used the
pre-b68c5af filter. Current `fastregion.region(..., positivity_only=True)` rejects 5,721 of the
34,071 graphs, e.g. idx 0, 18, 65, 144, 161; graph 65 has every `a_i < 0` on `beta > 2`.
Only 28,350 pass the strict filter (s5_survivor_filter.py).
Soundness is NOT affected: both filters are sound by Lemma A.2, and the loose one only keeps extra graphs.
But the count 34,071 cannot be reproduced by the current code. The set is not closed under
D-equivalence: a loose-only graph below idx 318k can have D-equivalent partners above 318k that were
filtered out, which is why 3,050 in-twin moves "leave" the file. And the 21,800-class count is a count
on a mixed set. On the strict set it is 17,307 classes (13,278 singletons, max 259).

---

## 9. Scripts (exact; run with python3 from anywhere; they chdir to the repo)

common.py (helpers: G1, G2, N, decoding of indices) ·
s1_relabel.py (min-distance relabelling) · s2_g1g2.py (one-edge move, vector equality of a) ·
s3_twin_components.py (twin union-find on the 34,071 file) · s4_debug_leave.py · s5_survivor_filter.py (anomaly) ·
s6_strict_components.py (strict set, D-classes, twin components) · s7_split_examples.py ·
s8_rank1_moves.py (rank-one J.1 moves) · s9_invariant_check.py (Theorem J.2 partition identity) ·
s10_fp_g1g2.py, s10b_fp_g1g2_diffs.py (Task 3 on G1/G2) · s11_fp_moves.py (FP preservation per move type, n=5) ·
s12_split95.py · s13_bidir_components.py · s14_invariants_and_fp6.py · s15_fp6_strict.py · s16_tfc_profile.py ·
s17_bidir_n6_and_G1.py (J.4 at n=6; transpose; the similarity S) · s18_weaker_invariants.py ·
s19_twin_eigen.py (Remark J.4': 120,146 twin pairs in supports, 0 failures; 59,410 non-adjacent-twin supports unstable by Hurwitz, none stable).
