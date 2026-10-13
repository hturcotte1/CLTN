# Lemma J: when do two graphs have the same full-support data D(A)?

Prover notes, 2026-10-13. WORK IN PROGRESS -- saved incrementally.
Conventions exactly CLAUDE.md §2 and notes/lemma_I_spectral.md §0:
`A[i][j] = 1` iff `i -> j`; `N := J - I - A^T` (`N[i][j] = 1` iff `i != j` and `j -/-> i`);
`C := I + beta N`, `K := C + tJ`, `a := adj(C) 1`, `d := det C`, `phi(mu) := det(mu I + C)`,
`omega(mu) := 1^T adj(mu I + C) 1 = det(mu I + C + J) - det(mu I + C)`.
`D(A) := (d, phi, omega, multiset{a_i})`, all polynomials in `beta` (and `mu`), as computed by
`fastregion.beta_data(A, range(n))`. `chi_X(z) := det(zI - X)`. `Kr(X; v_1..v_r)` := the smallest
`X`-invariant subspace containing `v_1..v_r` (= span of all `X^k v_j`). Vectors are columns.
Scratch scripts (all exact; re-runnable): /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/ .

Labels: **PROVED** (proof in this file, linear algebra only), **CHECKED** (exact computation, script
named), **OBSERVED** (data only, no proof), **OPEN**.

---

## 0. Summary

| item | status |
|---|---|
| J.0 block lemma: common invariant subspace containing 1 + cospectral quotients => equal D, equal `a` as vectors | PROVED |
| J.1 rank-one moves `N' = N + y x^T`, `x ⊥ Kr(N; 1, y)`; in-twin switching is the case `x = e_u - e_v` | PROVED |
| G2 = G1 with one edge moved between the bidirected in-twins 1, 2 | CHECKED (s1, s2) |
| J.2 characterisation: `D(A) = D(B)` iff some relabelling of B is cospectral with A and has the same in-walk vectors `(A^T)^k 1` for all k | PROVED; partition identity CHECKED on all 34,071 survivors (s9) |
| J.3 necessary conditions (in-degree sequence, #edges, #2-cycles, #3-cycles, chi_A, chi_(A|perp), total walk counts) and non-invariants (out-degree sequence, FP structure in general) | PROVED / CHECKED |
| J.4 bidirected in-twin switching preserves FP(G), every fixed point and its Jacobian spectrum, for every support and every (eps, delta, theta) | PROVED |
| (see later sections for FP results on other moves, the generating-set question, and an anomaly in the survivor set) | |

---

## 1. The G1 / G2 relation (CHECKED: s1_relabel.py, s2_g1g2.py)

Relabel G2 by p = (4,3,5,2,1,0) (`graphs.permute`: new node i = old node p[i]). Then
```
G1  = 011011 / 101011 / 110101 / 111000 / 011000 / 000110
G2p = 011011 / 101001 / 110111 / 111000 / 011000 / 000110
```
differ in exactly two entries (the minimum over all 720 relabellings; the only other relabelling at
distance 2 is the image under swapping the roles of 3 and 4): **the edge 1->4 of G1 is replaced by
2->4.** In G1 nodes 1 and 2 are bidirected (1<->2) and have the same in-neighbours among the others
(`in(1) \ {2} = in(2) \ {1} = {0,3,4}`). So `N(G2p) = N(G1) + e_4 (e_1 - e_2)^T`, and `x = e_1 - e_2`
satisfies `x^T N(G1) = 0` and `x^T 1 = 0`. Exact: `beta_data(G1)` and `beta_data(G2p)` agree entry for
entry -- `a` agrees as a VECTOR in this labelling, not only as a multiset.
Out-degrees change only at the two twins (1: 4 -> 3, 2: 4 -> 5), in-degrees change nowhere
(node 4 trades in-neighbour 1 for in-neighbour 2). This is the whole difference between Geneson's
graph and hit 3.

---

## 2. The block lemma (PROVED)

**Lemma J.0.** Let `N, N'` be real `n x n` matrices and `V` a subspace with `1 in V`, `N V ⊆ V`, and
`N v = N' v` for all `v in V` (so also `N' V ⊆ V`). Let `N̄, N̄'` be the maps induced on `R^n / V`.
If `chi_N̄ = chi_N̄'` then, for all real `beta, mu, s`:
  (a) `det(zI - N - sJ) = det(zI - N' - sJ)` as polynomials in `(z, s)`;
  (b) `d' = d`, `phi' = phi`, `omega' = omega` (as polynomials in `beta, mu`);
  (c) `adj(I + beta N') 1 = adj(I + beta N) 1` (as vectors of polynomials in `beta`); more generally
      `adj(I + beta N' + tJ) 1 = adj(I + beta N + tJ) 1` and `det(mu I + I + beta N' + tJ) = det(mu I + I + beta N + tJ)`.

*Proof.* Choose a basis whose first `r = dim V` vectors span `V`. In it
`N = [[X, Y], [0, Z]]`, `N' = [[X, Y'], [0, Z']]` with the same `X` (agreement on `V`) and
`chi_Z = chi_N̄ = chi_N̄' = chi_Z'`. Because `1 in V`, `J = 1 1^T` maps `R^n` into `V`, so in this basis
`J = [[J_1, J_2], [0, 0]]`. Hence for every scalar `c, s`, `cI + N + sJ = [[cI + X + sJ_1, *], [0, cI + Z]]`
and the same with primes, with identical diagonal blocks up to the cospectral `Z, Z'`; so
`det(cI + N + sJ) = det(cI + X + sJ_1) det(cI + Z) = det(cI + N' + sJ)`. Taking the polynomial
identity in `c, s` gives (a). For (b): `d = det(I + beta N)`, `phi = det((mu+1) I + beta N)`,
`phi + omega = det((mu+1) I + beta N + J)`; each is `det(cI + beta(N + sJ))` with `s = 0` or `s = 1/beta`,
i.e. `beta^n det((c/beta) I + N + sJ)` for `beta != 0`, equal for `N` and `N'` by (a); polynomial
identities in `beta` that hold for all `beta != 0` hold identically. For (c): fix `beta, t` with
`M := I + beta N + tJ` invertible (all but finitely many `beta` for fixed `t`, since `det M` is a nonzero
polynomial in `beta`: it is `1 + nt` at `beta = 0`... for `t > -1/n`; for the identity it suffices
to vary `beta` and `t` over an open set). `V` is `M`-invariant and `M|_V = M'|_V` (both `beta N + tJ`
and `I` agree on `V`). `M|_V` is injective, hence invertible on `V`, so the unique solution of
`M u = 1` lies in `V` (`u = (M|_V)^{-1} 1`) and then also solves `M' u = 1`. With `det M = det M'`
(by (a)) this gives `adj(M) 1 = det(M) u = adj(M') 1`. Both sides are polynomial in `(beta, t)` and agree
on a nonempty open set, hence identically. ∎

Note what this says dynamically (full support, any `eps, delta, theta`): `I - W = eps K`, so the full-
support equilibrium `x* = theta K^(-1) 1 / eps` is **the same vector** for `N` and `N'`, and the
Jacobian `-I + W = -eps K` has **the same characteristic polynomial**. So the full-support verdict
(positivity, stability, the whole 2-D region, every critical value) coincides, not just the limit.

---

## 3. Rank-one moves and in-twin switching (PROVED)

**Lemma J.1.** Let `x, y in R^n` with `x ⊥ N^k 1` and `x ⊥ N^k y` for all `k >= 0`
(i.e. `x ⊥ Kr(N; 1, y)`). Then `N' := N + y x^T` satisfies the hypotheses of J.0 with
`V = Kr(N; 1, y)` and `Z' = Z`; hence `D(N') = D(N)` with `a' = a` as vectors.
Moreover `x ⊥ Kr(N'; 1, y)` (the same space), so the move is reversible.

*Proof.* For `v in V`, `N' v = N v + y (x^T v) = N v`. So `V` is `N'`-invariant and `N, N'` agree on it.
`N' - N = y x^T` has range in `span(y) ⊆ V`, so it induces 0 on `R^n / V`: `Z' = Z`. Apply J.0. ∎

Conversely, for a rank-one change the condition is necessary: if `D(N + y x^T) = D(N)` with `a` equal as
vectors then `x ⊥ Kr(N; 1, y)`. (From `a' = a`, `d' = d`: `C'^(-1) 1 = C^(-1) 1`; Sherman-Morrison gives
`C'^(-1) 1 - C^(-1) 1 = -beta C^(-1) y (x^T C^(-1) 1) / (1 + beta x^T C^(-1) y)`, so
`x^T C^(-1) 1 = sum_k (-beta)^k x^T N^k 1 ≡ 0` unless `C^(-1) y ≡ 0`, impossible; and
`chi_(N') = chi_N - x^T adj(zI - N) y`, so `x^T (zI - N)^(-1) y ≡ 0`, i.e. `x ⊥ N^k y` for all k.)

**In-twins.** Call `u != v` *in-twins* if `A[k][u] = A[k][v]` for every `k ∉ {u, v}` and
`A[u][v] = A[v][u]` (the pair is bidirected or non-adjacent). Then rows `u, v` of `N` agree off columns
`{u, v}` and `N[u][v] = N[v][u] =: c`, so `x := e_u - e_v` satisfies `x^T N = -c x^T` and `x^T 1 = 0`.
Hence `x ⊥ N^k 1` and `x ⊥ N^k y` for every `y` with `y_u = y_v` (`x^T N^k y = (-c)^k x^T y`).

**In-twin switching.** For a vertex `w ∉ {u, v}` with exactly one of `u -> w`, `v -> w`, replace that
edge by the other one. In `N` this is `N' = N ± e_w (e_u - e_v)^T`, a move of Lemma J.1 (`y = ±e_w`,
`y_u = y_v = 0`). So it preserves `D(A)`, with `a` equal as a vector. `u, v` remain in-twins in the new
graph (their in-neighbourhoods are untouched), so moves at different `w` commute and can be combined:
**for in-twins `u, v`, the out-edges of `u` and `v` to third vertices can be redistributed arbitrarily
subject to keeping, for every third vertex `w`, the number `A[u][w] + A[v][w]` fixed.** Each in-degree
is preserved; only the out-degrees of `u, v` change.

The move from G1 to G2p is exactly this with `(u, v, w) = (1, 2, 4)` and `c = 0`.
