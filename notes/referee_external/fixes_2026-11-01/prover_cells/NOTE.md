# Issue 1 of the external report of 2026-11-01: SFP points on lower-dimensional cells

Report: `notes/referee_external/geneson_report_2026-11-01.txt`, point 1 (paper Section 3.1, `paper.tex` l. 427-428).
Script: `notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py`
(output `geneson_section_check.json` in the same folder; about 70 s single-threaded; re-run with
`cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py`).
All decisions in the script are exact (sympy integers and rationals, `fractions.Fraction`). The decimals below are
display approximations only.

**Verdict: the report is correct.** The sentence "lower-dimensional cells are boundary points and are never SFP" is
false. No count, verdict or threshold in the paper depends on it. The implementation never used it. The project's
own proof note (`notes/lemma_Bprime_2d.md`, Corollary B''(4)-(5)) already states the correct description, so the
error is confined to the paper's one-line summary.

## 1. Exact reproduction of the example

The script rebuilds everything from the definitions in paper Section 2: `N = J - I - A^T`, `C = I + beta N`,
`K = C + tJ`, `p(mu) = det(mu I + K)`, `a = adj(C) 1`, `Delta_k = det[c_{2i-j}]`, `L = t(beta-2) - 1`.
It does not use the project engines, except for the cross-checks in (e) and (f).

- **(a) Graph.** `G_1` has the Table 1 rows `011000 000011 000111 110011 110101 101110`. Its `p(mu; beta, t)` is
  identical to that of Geneson's labelling `A_ours` in `certificates/stage0_geneson_2607.21396.json`.
- **(b) `Delta_5` on the legality boundary.**
  - `Delta_5` has degree 15 in `beta` and degree 5 in `t`.
  - `P(beta) := (beta-2)^5 Delta_5(beta, 1/(beta-2))` has degree 20, and `Res_t(Delta_5, L) = -P` exactly.
  - `P` factors as (degree 6)(degree 14). Both sympy `count_roots` and an independent Sturm sequence over
    `Fraction` find **exactly one** root in `(17/8, 15/7)`. As a check, `P(17/8) > 0 > P(15/7)`.
  - The root `beta_0` belongs to the degree-14 factor
    `43b^14 - 498b^13 + 2780b^12 - 9966b^11 + 25641b^10 - 50226b^9 + 77654b^8 - 96910b^7 + 98882b^6 - 82792b^5 + 56784b^4 - 31408b^3 + 13760b^2 - 4352b + 1024`.
  - This polynomial is the one recorded as the `beta_critical` entry in `[17/8, 15/7]` of the `G_1` certificate
    (`certificates/region2d_n6_shard2of3_HITS.jsonl`, idx 1000255). So `beta_0` is the paper's `beta_1` (l. 600).
  - Isolating interval after 60 exact bisections:
    `beta_0 in (68668602371678778193/32281802128991715328, 137337204743357556387/64563604257983430656)`
    (width about 1.5e-20). Display values: `beta_0 ~ 2.1271613678` and `1/(beta_0 - 2) ~ 7.8640`.
- **(c) Positivity on the whole strip.** Substitute `beta = (17/8 + (15/7)u)/(1+u)`, multiply by `(1+u)^deg`
  and clear denominators. At `t = 100`, every coefficient is strictly positive for all six adjugate coordinates
  `a_1..a_6` (`beta`-degrees 4, 3, 4, 4, 2, 2; they do not depend on `t`) and for all six Hurwitz determinants
  `Delta_1..Delta_6` (`beta`-degrees 0, 3, 6, 10, 15, 20). The same holds at `t = 1000`.
  - Two further checks, beyond what the report asks:
    - The bivariate version, with `t = 100 + s`, has all coefficients positive for `Delta_1..Delta_6`. So
      **every point of `[17/8, 15/7] x [100, inf)` is an SFP**.
    - The coefficient test is a sufficient condition only. It is not conclusive at `t = 10`, so nothing is
      claimed there from this test.
- **(d) Legality.** `100(beta_0 - 2) > 100(17/8 - 2) = 25/2 > 1`. Conclusion: by Lemma A.2, with condition (2)
  vacuous, **the full support of `G_1` is an SFP at `(beta_0, 100)`**. The point lies on no boundary curve:
  `a_i`, `D` and `Delta_5` are nonzero there, and `L > 0`.
- **(e) Cross-check with the project verifier.** `ctln.analyze_support`, which builds `W` from the project conventions (§2) by
  an independent code path, was run at `eps = 1/101` (`t = 100`) and `beta` in {17/8, both ends of the isolating
  interval, 15/7}. It returns positive and stable at all four. At both ends of the isolating interval and
  `t` in {8, 9, 10, 20, 50}, all legal, it also returns positive and stable.
- **(f) `beta_0` is a projection root of our decomposition.**
  - `region2d.projection` on `G_1` gives the curves `F = {D (deg_t 1), Delta5_f0 (deg_t 1), Delta5_f1 (deg_t 4), L}`.
  - It gives 11 projection polynomials: `lc_D, lc_Delta5_f0, lc_Delta5_f1, disc_Delta5_f1, lc_L, res_D_Delta5_f0,
    res_D_Delta5_f1, res_D_L, res_Delta5_f0_Delta5_f1, res_Delta5_f0_L, res_Delta5_f1_L`.
  - The minimal polynomial of `beta_0` divides **exactly one** of them, `res_Delta5_f1_L`.
  - The repository's `Delta_5` equals the one computed here (ratio 1).
  - So the line `beta = beta_0` is a lower-dimensional cell (a vertical section) of the paper's own CAD, and it
    passes through the open SFP region.

**The whole legal section is SFP (proof, from the computed facts in (f)).**
1. The interval `I = (2, beta_0)` contains no projection root, since the first joint isolating interval is
   `[17/8, 15/7]`, so Lemma B''' applies on `I`. The recorded cell at the sample 33/16 has a single `t`-cell
   above `t_min = 16`. Hence no curve of `F` has a real `t`-root above `1/(beta-2)` for any `beta` in `I`.
2. At `beta_0`, every `lc` and the discriminant are nonzero, and every resultant except `res_Delta5_f1_L` is
   nonzero. So each real `t`-root of each curve at `beta_0` is simple and lies on a real continuous branch that is
   defined on both sides (implicit function theorem). Its values for `beta < beta_0` are below `1/(beta-2)`.
   Letting `beta` increase to `beta_0`, the root is at most `1/(beta_0-2)`.
3. So the open ray `{beta_0} x (1/(beta_0-2), inf)` meets no curve. It does not meet `a_i = 0` either, because
   `beta_0` lies inside the positivity cell.
4. The ray is therefore a connected subset of `R \ Z_0`. By Lemma B' its verdict is constant, and it is SFP at
   `t = 100`. So it is SFP at every legal `t`.

This agrees with the paper's own description of the `G_1` region (l. 599-601), which includes `beta = beta_1`. It
contradicts l. 428, so the paper was also internally inconsistent here.

## 2. How the implementation and the certificates treat lower-dimensional cells

`region2d.py` **does not decide, does not mark, and does not record** cells of dimension 0 or 1. They are
neither "undecided" nor "boundary"; they simply do not appear.

- **`beta`-cells.**
  - `beta`-cells are the open gaps between consecutive *isolating intervals*.
  - The sample is the midpoint of the gap, a rational point inside the true open cell `(r_k, r_{k+1})`.
  - The vertical lines over the roots `r_k` are never sampled.
- **`t`-cells.** At the sample, the `t`-cells are likewise the gaps between `t`-isolating intervals, so the
  sections (curve graphs) are not sampled either.
- **What a certificate records.**
  - `beta_critical`: the minimal polynomial and isolating interval of every projection root.
  - `cells`: the 2-cells with their verdicts.
  - `stable_cells`: the SFP 2-cells, given by rational hulls `(lo, hi, sample, t-range at the sample)`.
  - `verdict`: `stable_somewhere` / `never` / `escalate`, a function of the 2-cells only.
- **Example, `G_1`.** The certificate's recorded cells are `(2, 17/8)` and `(15/7, 1029/475)`.
  `any_recorded_cell_contains_beta0 = false`: no recorded hull mentions any `beta` in `[17/8, 15/7]`.
- **What a stable entry certifies.** It certifies the whole true 2-cell, whose `beta`-ends are the algebraic
  roots listed in `beta_critical`. It does not certify only the rational hull. `n7/qstar.py` reads entries this
  way (docstring S0-S2).
- **Reading the certificates.** The recorded cells certify the open SFP set **up to lower-dimensional
  sections**. Reading a certificate as "the SFP set is the union of the listed hulls" would be an inner
  approximation, and that reading is not used anywhere. A point on a section is decided pointwise:
  - by the three-valued exact test at that point (for an algebraic `beta`, by an exact sign certificate such as
    (c) above, or by Lemma B along the section);
  - equivalently, by Corollary B''(5): it is SFP iff it lies off `Z_0` and in the closure of the SFP 2-cells.

  It is never decided by cell membership.

## 3. Nothing in the paper depends on the false sentence

Notation: `S = S_H` is the SFP set of one (graph, support). It is open (Lemma `lem:open`), and the union of the
2-cells is dense (lemma_Bprime_2d.md §2; every lower-dimensional CAD cell is nowhere dense).

- **"never" verdicts are sound.** These cover Theorem `thm:n5`, Theorem `thm:n6` ("every other six-node graph is
  certified never stable"), Theorem `thm:n7` and stage 6 of Table `tab:chain`, the 389 / 705 five-node records
  and the 34,067 six-node records.
  - If `S` is nonempty then, being open, it meets some 2-cell `U`.
  - The verdict is constant on `U` (Corollary B''(2)), so `U` is an SFP cell.
  - Contrapositive: no SFP 2-cell implies `S` is empty. Sections are never needed.
- **"stable somewhere" and the counts.** An SFP 2-cell is a nonempty subset of `S`. These are the exactly 3
  non-complete six-node graphs, the 124 seven-node graphs and `K_k`.
- **Suprema `q_6^*`, `q_7^* = rho`, `q_8^*(inf)` and the per-graph suprema.** The thresholds are Theorems
  `thm:q6`, `thm:q7` and `thm:q8`; there is no label `thm:thr6`.
  - `S` lies inside the closure of `S_c`, the union of the SFP 2-cells, and `q` is continuous. So
    `sup_S q = sup_{S_c} q`. This is `qstar.py` S0.
  - The rest of the computation uses only 2-cells and their algebraic endpoints.
  - Sections can neither raise nor lower the supremum.
- **"not attained".** This follows from openness alone. `dq/dbeta = t/(t+1) > 0`, so `q` has no local maximum
  on an open set, and `Bad(n)` is a finite union of open `S_H`. Sections play no role.
- **Pointwise region descriptions.**
  - **`G_1`/`G_2`** (l. 599-601): the description includes the section `beta = beta_1 = beta_0`, and §1 proves that
    section is SFP. It is consistent with the corrected statement.
  - **`H_3`** (l. 601-604: "SFP iff `2 < beta < beta_1'`, every legal `t`"): its positivity cell contains no projection
    root (1 `beta`-cell and 1 `t`-cell in the certificate), so there is no section inside. At `beta_1'`, `a_5 = 0`.
  - **`E_0`** (l. 713, "bounded in `t` at every `beta`"): a single `beta`-cell `(2, beta_1')` with no interior
    projection root (`certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json`).
  - **Bulge statements and cell counts** (Proposition `prop:bulge`; 258 / 40 / 298 stable cells) are statements
    about 2-cells by definition.
- **One-parameter engine** (Lemma `lem:B`, the `eps`-slices, the `t = inf` limit maps at n = 6, 7, 8, and the
  Lemma-R curve test). Here the critical set consists of members of the boundary family itself: `a_i`, `o_k`,
  `D_t`, `Delta_{m-1}`, or their `t = inf` and curve analogues. So "false at those roots" in Lemma B is
  **correct**. The defect is specific to the 2-D projection: leading coefficients, discriminants and resultants
  are not boundary polynomials.

## 4. Corrected statement

Let `S_c` be the union of the SFP full-dimensional cells and `Z` the union of the zero sets of the `a_i`, `o_k`,
`D` and `Delta_{m-1}`. Then

    S = ( cl(S_c) ∩ {L > 0} ) \ Z          (lemma_Bprime_2d.md, Corollary B''(5); proof: S open + Lemma B')

Lower-dimensional cells lying on `Z` or `{L = 0}` are never SFP. A vertical section above a root of a leading
coefficient, discriminant or resultant need not lie on `Z` and can be SFP, as in §1. In particular
`S != ∅` iff some full-dimensional cell is SFP.

## 5. Proposed text edits (old -> new). The document was not edited.

### 5.1 `notes/paper/paper.tex`, l. 427-428 (required)

OLD
```
of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact
rational sample; lower-dimensional cells are boundary points and are never SFP.
```
NEW
```
of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact
rational sample; lower-dimensional cells are not sampled. Points on a boundary curve or a positivity root are never SFP, but
a vertical section above a root of a projection polynomial need not lie on any boundary curve and can pass through the open
SFP region: as pointed out by J.~Geneson (personal communication), for $G_1$ of Table~\ref{tab:n6} the full support is an
SFP at $(\beta_0,100)$, where $\beta_0$ is the unique root in $(17/8,15/7)$ of $\Delta_5(\beta,1/(\beta-2))$, a root of the
projection polynomial $\operatorname{Res}_t(\Delta_5,L)$.\footnote{\cf{notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py}.}
Since the SFP set $S$ is open (Lemma~\ref{lem:open}) and the full-dimensional cells are dense, they determine it exactly:
$S=\bigl(\overline{S_c}\cap\{L>0\}\bigr)\setminus Z$, where $S_c$ is the union of the SFP full-dimensional cells and $Z$ the
union of the zero sets of the $a_i$, $o_k$, $D$ and $\Delta_{m-1}$; in particular $S\neq\emptyset$ iff some full-dimensional
cell is SFP, and every ``never'' verdict, count and supremum of $q$ below is a statement about full-dimensional
cells.\footnote{\cf{notes/lemma_Bprime_2d.md}, Corollary~B$''$(4)--(5).}
```

### 5.2 `notes/paper/paper.tex`, l. 397 (recommended; the decomposition's walls also include vertical lines over projection roots)

OLD
```
are decided instead on open cells of a decomposition whose walls are the zero sets of a proved boundary set.
```
NEW
```
are decided instead on open cells of a decomposition whose walls contain the zero sets of a proved boundary set.
```

### 5.3 `notes/paper/paper.tex`, l. 600 (optional; records the section in the `G_1` region)

OLD
```
$\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$: these graphs are stable at \emph{every} legal $t$ when $\delta/\eps$
```
NEW
```
$\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$, the root there of $\operatorname{Res}_t(\Delta_5,L)$ (the
region also contains the whole legal vertical section $\beta=\beta_1$; Section~\ref{sec:machinery}): these graphs are stable at \emph{every} legal $t$ when $\delta/\eps$
```

**Build check.** All three edits were applied to a scratch copy of `notes/paper` and built with
`pdflatex`, `bibtex`, `pdflatex`, `pdflatex`. The result has 0 errors, 0 undefined references and 0 overfull
boxes, and is still 33 pages. The macros used (`\operatorname`, `\overline`, `\cf`) are already available in the
preamble.

### 5.4 Other repository documents

Search: `grep -rn -i 'never SFP\|never a stable\|lower-dimensional\|lower dimensional\|boundary point\|critical values themselves'`
over `notes/*.md`, `notes/paper/*.tex`, `README.md` and `n7/*.md`.

- **No other current document repeats the false claim.** `REPORT_n6.md`, `REPORT_n7.md`, `REPORT_n8.md`,
  `FINAL_SUMMARY.md`, `n6_map.md` and `README.md` contain no such sentence. `REPORT_n6.md` l. 191-193 correctly
  says the note proves "correctness of deciding only full-dimensional cells via openness".
- `notes/lemma_Bprime_2d.md` is **correct as written** and is the source of the fix. See Corollary B''(4)-(5),
  l. 293-306 ("A point of R on a lower-dimensional piece is stable iff it lies off Z_0 and in the closure of the
  stable cells"), §3 l. 576 and §4 l. 698. No edit is needed.
- `notes/lemmas.md` l. 257 and l. 721-722 ("critical values themselves are never SFP (B(d))") concern the
  one-parameter engine. Its critical values are roots of boundary polynomials, so these lines are correct and
  need no edit.
- `notes/paper/paper_draft_before_check.tex` l. 389 contains the same false sentence. It is the frozen
  pre-check snapshot of the draft. I recommend leaving it unchanged as a historical record. If it is to be kept
  consistent, apply 5.1 verbatim (its l. 388-389 have the same text).
