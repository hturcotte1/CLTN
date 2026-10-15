
---

## 13. Continuation session: re-run audit, independent re-proof, and remarks on (a) and (c)

**13.1 Re-run of the previous attempt (CHECKED).** Before anything was re-run, every JSON and log was copied to
`PF/orig_outputs/`. Then every proof and certification script was re-run from scratch, sequentially on one core
(`rerun_chain2.sh`; logs in `PF/rerun_logs/`). `diff_outputs.py` compares each regenerated JSON with the original, ignoring
only `seconds` fields.

| script | exit | regenerated output vs original |
|---|---|---|
| `theorem_A.py` | 0 | `theorem_A_certificate.json` identical |
| `conjP_refutation.py` | 0 | identical |
| `negative.py` | 0 | `negative_certificate.json` identical |
| `ceilings.py` | 0 | identical |
| `positivity_ranges.py` | 0 | identical |
| `pointwise_thresholds.py` | 0 | identical |
| `check_lemmaP.py` | 0 | identical (ALL CHECKS PASSED) |
| `check_X2.py` | 0 | identical (61 positive / 19 non-positive agree with ctln) |
| `theorem_B.py R`, `Q Q2`, `Y` | 0 | all three certificates identical |
| `Tk_table.py` | 0 | identical |
| `crosscheck_region.py` | CROSSCHECK_STATUS | CROSSCHECK_DIFF |
| `verify_all.py --full` | VERIFY_STATUS | VERIFY_DETAIL |

Not re-run: `scan_k0_catalogue.py` (about 18 min, float preselection; §9 is supplementary) and the HEURISTIC scripts.

**13.2 Independent re-proof of Theorem A step 3 and Theorem SJ (PROVED; `independent_A.py` -> `independent_A.json`).** This
uses a separate code path, with no sympy and only `fractions.Fraction` plus the repository's
`ctln.charpoly / det_exact / hurwitz_stable`.
- R is entered from its rows and compared with the digraph6 decoding. `chi_{D+aJ}` and `chi_{D - e_2 1^T}` are computed by
  Faddeev-LeVerrier and cross-checked against ctln's independent interpolation routine.
- `Psi_k = z chi_{D+aJ} - a k chi_L`, and `h_k(x) = Psi_k(-r - x)`.
- Each coefficient of `h_k` is affine in k, so the j-th Hurwitz minor is a polynomial of degree at most j in k. It is
  recovered exactly by Lagrange interpolation at `k = 0..j`, and the recovery is confirmed at 6 further k.
- All 8 minors have **only positive coefficients**, and they agree with the polynomials stored in
  `theorem_A_certificate.json` at k = 0..11.
- As a redundant sweep, `ctln.hurwitz_stable(h_k)` is True for every k = 0..2000.
- The same interpolation in c re-proves Theorem SJ's locus statement: all 7 minors of the undivided `det(zI - D - cJ)` have
  positive coefficients in c.
- Positivity (`w_R > 0`, `sigma_S = 381654/295819 > 1`) is recomputed with `ctln.solve_exact`.

**13.3 Question (a) for the multi-class families (CERTIFIED; `Tk_extra.py` -> `Tk_extra.json`).** Stable t-sets are computed exactly at fixed beta. Endpoints are
real roots of the stated minors, isolated exactly; decimals are display only.

| family, beta (legal boundary) | k: T_k |
|---|---|
| Q2, beta = 3 (t > 1) | k = 0, 1, 2: no stable legal t. Then T_3 = 13.837, T_4 = 2.0932, T_5 = 1.1321, and **T_k = 1 for k = 6, 7, 10, 20** |
| Y, beta = 3 (t > 1) | k = 0, 1, 2, 5, 10: no stable legal t. Then T_15 = 3.1002, T_20 = 1.1461, T_21 = 1.0181, and **T_k = 1 for k = 22, 30** |
| Y, beta = 9/2 (t > 2/5) | k = 20, 50, 100: no stable legal t. Then T_150 = 3.4497, T_200 = 1.3961, T_300 = 0.63779, T_400 = 0.41340, and **T_502 = 2/5** |

Every stable set computed is a single ray `(T_k, inf)`, with no bounded windows, unlike F at beta = 21/10. In each row the
sequence decreases to the legal boundary. **Theorem B's k0 is sharp at beta1 = 3 for both families:**
- Q2: `T_5 > 1` while `T_6 = 1` (Theorem B: k0 = 6).
- Y: `T_21 > 1` while `T_22 = 1` (Theorem B: k0 = 22).
- At `beta = 9/2`, `T_400 > 2/5` shows that the region theorem's k0 = 502 cannot be lowered below 401 there.

This is the bounded-sequence answer to question (a): along each family at fixed beta below its ceiling, `T_k(beta)` is
eventually constant, equal to the legal boundary (Theorem FP(iv)).

**13.4 Why idea (c) cannot succeed where families exist, and what remains (PROVED remarks).**
- Families exist for every legal point with `beta < beta_Y`, so no bound `n <= f(eps, delta)` on stable non-clique full
  supports holds anywhere in that region.
- The known necessary conditions are scale-free, and the families satisfy them with growing slack:
  - The `e_2` bound (`B > Nn`: more bidirected than non-adjacent pairs) holds easily. In `R_k`, `B >= k(k-1)/2`, while
    `Nn` is the number of non-adjacent pairs of R, which is 0.
  - The trace identity `tr(D' + aJ) = an` is balanced by the single escaping Perron root (§12, structural summary).
- A size bound at fixed `(eps, delta)` could therefore only hold at large beta. There it would say
  `Lambda ⊆ {beta < beta*}`. That is the fixed-parameter form of the open question "is `sup_n q_n*` finite?": a point with
  `q > sup_n q_n*` carries no non-clique stable full support of any size. This remains **OPEN** here.
- Proven caps that apply to restricted constructions:
  - Joins and self-joins of an m-node core: `beta < m - 1` (§8, Lemma M2).
  - Sink-class families without a join class: `beta_inf <= m` (§8).
  - Blow-ups whose fastest classes are not complete fail at every point (Lemma U'').
- Mixtures of a join class with non-join classes over an m-node core are not covered by a proved cap.
  - HEURISTIC (`heur_mixed_cap.py` -> `heur_mixed_cap.json`, floats): exhaustive over all D for m <= 4, plus 3000 / 1500
    random D for m = 5 / 6, with one non-join class of weight w and the rest join, and positivity checked. The largest
    `beta_inf / m` found is 0.67, 0.50, 0.45, 0.33 for m = 3, 4, 5, 6. Y reaches 0.62 at m = 8.
  - Nothing suggests mixtures beat the cap `beta_inf <= m`.
  - **OPEN (conjecture):** for every sink-class family over an m-node core, `beta_inf <= m`. If true, together with Lemmas
    U, U'' and M2 it would mean that every infinite family of the forms studied here, at large beta, needs a core of
    at least beta nodes.
