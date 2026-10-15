
---

## 14. Every legal point carries an infinite family (PROVED modulo three cited lemmas; instance CERTIFIED at (1/2, 2))

Added at the end of the continuation session. The project summary (`notes/FINAL_SUMMARY.md`) records Theorem U as refereed.
It was proved by the parallel workflow in
`/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/notes_prover_q.md` and states
`q_n* -> infinity`. It is built from "arms-forest" cores with one sink class. Combined with Theorem FP(i), its ingredients
settle O1. Script: `extend_lambda.py` -> `extend_lambda.json`.

**Cited (from notes_prover_q.md §§2-4; not re-proved here).**
- **(T) Lemma T.** The arms-forest core built from `c = [1, 1, c_2, ..., c_d]` has `chi_{D - e_c2 1^T} = z^(m-d) p`, with
  `p = sum c_j z^(d-j)`. It is proved there by the matrix determinant lemma, and **re-verified here exactly** on the 65-node core
  `c = [1,1,11,7,10,3,1]` (E1).
- **(S) Lemma S (Rouché).** Fix `N >= 3` and `s = 1/(2N)`. For every integer `W > W0(N) = 2N 4^N (N+1)^(2N-2)`, the
  coefficientwise rounding `p~_W` of `prod_{k=1..N} ((z+s)^2 + W^2 k^2)` has the following properties:
  - integer coefficients with `c~_0 = c~_1 = 1` and `c~_j >= 1`;
  - it satisfies the forest constraint;
  - all its roots lie in `-3/(4N) < Re z < -1/(4N)`.
- **(P') closed-form positivity (their Lemma P).** At `r = 1/beta`, the core is positive and the class is activated iff
  `lam > 0`, `lam fmax < 1`, `u_c0, u_c1 > 0` and `0 < u_c2 < 1`. Here:
  - `lam = r^3 (1 - r + r^2)/Den` and `Den = F0 - r F1 + r^2 F2 + r^2 + r^5`;
  - `F2 = c_2 - 1`, `F1 = sum_{odd j >= 3} c_j r^(3-j)` and `F0 = sum_{even j >= 4} c_j r^(4-j)`;
  - `u_c2 = 1 - r + r lam`, `u_c1 = 1 - r(1 - e2) - lam F2/r` and `u_c0 = 1 - r(1 - e1) - lam F1/r`.

  This matches §1's formulation `u > 0`, `U > 1`, `b^T u < 1` with `b = e_c2`. **Re-verified here by direct exact solves**
  (`ctln.solve_exact`) on the 65-node core at beta = 5, 11/2, 28/5 (E3).

**Theorem Λ.** For every legal `(eps, delta)` there is a non-clique arms-forest core C and the in-set S = all nodes but c2
such that the full support of `C[S, all; k]` is a stable positive fixed point at `(eps, delta)` for every `k >= k0`.
Hence the set Lambda of §8 (O1) is the **entire legal region**. At no legal point is the size of a stable non-clique full
support bounded.

*Proof.* Put `beta = 1 + delta/eps > 2`, `r = 1/beta`, `t = (1-eps)/eps`.
1. Choose an integer `N >= 3` with `N > 3 beta / 4`, so that `3/(4N) < r`.
2. For every `W > W0(N)`, (S) and (T) put every root of `chi_L = z^(m-d) p~_W` either at 0 or in `Re z > -3/(4N) > -r`.
3. **Positivity at this r for large W.** This step extends their Step 3, which is written for `r = 1/N`; the argument needs
   only a fixed `r in (0, 1/2)`.
   - Let `W -> inf` with N and r fixed. `c~_2N >= W^(2N) (N!)^2 - 1/2`. Every other coefficient is `O(W^(2N-2))`: for
     `j < 2N` the coefficient of `z^(2N-j)` uses at most N - 1 of the factors `W^2 k^2`. So `F1/F0`, `F2/F0`, `1/F0 -> 0`.
   - Hence `Den/F0 -> 1`, `lam -> 0+` and `lam F0 -> r^3(1 - r + r^2) < 1/8`, while `lam F1, lam F2 -> 0`.
   - `fmax <= max(F0, F1, F2, 1)`: a depth-1 arm's f is at most its exit's F, and `f_parent > f_child`. So
     `lam fmax -> r^3(1-r+r^2) < 1`.
   - `u_c2 -> 1 - r`, which lies in `(1/2, 1)`; `u_c1 -> 1 - r + r^2 > 0`; and `u_c0 -> 1 - r + r^2 - r^3 > 0`.
   - By (P'), the core is positive and the class activated at r for all large W. Fix one such W.
4. Theorem FP(i), at the fixed legal `(beta, t)`, gives k0 such that `C[S, all; k]` is an SFP for every `k >= k0`. Positivity
   is Lemma X2'; the "if" direction needs no rationality of beta (§1 remark).
5. `D != 0`, so no member is a clique. ∎

*Remarks.*
- k0 comes from Rouché and is not explicit. The cores are enormous for large beta (`m = 2 p~(1) - 1`, `W > W0(N)`). Theorem Λ
  is an existence statement and says nothing about the smallest size.
- It is consistent with Geneson's theorem, which only excludes `n <= (q + 2 - 3eps/2)/(1 - eps/2)`.
- Status: **PROVED conditional on (T), (S), (P')**. (T) and (P') are re-verified here on an instance, and (S) is cited as
  refereed. The extension of their Step 3 to arbitrary fixed r is new here and short. It should be refereed before it is
  quoted.

**Exact instance beyond `beta_Y` (CERTIFIED, reduced; `extend_lambda.py`).** Take the 65-node core `c = [1,1,11,7,10,3,1]`,
`p = z^6 + z^5 + 11z^4 + 7z^3 + 10z^2 + 3z + 1`.
- (E2, E3) At beta = 5, 11/2 and 28/5, every root of p has `Re > -1/beta` (exact Hurwitz test of `p(-1/beta - x)`), and the core
  is positive with the class activated (exact solve). So by FP(i), every legal point with one of these betas carries the
  family. All three exceed `beta_Y = 4.9607`.
- At beta = 6 the shifted polynomial has a zero coefficient. That is **UNDECIDED** (boundary), not unstable. It fits the trace
  constraint: the six roots of p sum to -1.
- (E4) At `(eps, delta) = (1/2, 2)` (beta = 5, t = 1, q = 2), Lemma P gives `Psi_k = z^59 g_k` with `deg g_k = 7`. The 7 Hurwitz
  minors of `-g_k(-1/5 - x)` are interpolated exactly in k. **Every member with `k >= 1581` (n >= 1646) is a stable positive
  non-clique full support**, and the verdict at k0 is "all Hurwitz minors positive". The minors are stored in
  `extend_lambda.json`.
- (E5) This validates the reduction on an explicit graph. On the 66-node member (k = 1) at `(1/2, 2)`, `ctln.analyze_support`
  gives positive = True (as E3 + X2' predict) and stable = False (Hurwitz minor 7 negative). This matches the Lemma P
  verdict (k = 1 < k0).
- An instance check by `ctln.analyze_support` at n >= 1646 is computationally out of reach (n = 107 takes 13 min), so this
  instance is certified through Lemma P (proved, §1) rather than by ctln.
