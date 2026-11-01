# Applied fixes for the external correctness report of 2026-11-01

Report: `notes/referee_external/geneson_report_2026-11-01.txt` (three points: lower-dimensional CAD cells can be SFP;
Theorem 5.4 invoked at `nu_inf = 1/2` in Lemma 7.4; the primitive proportions of Section 6.6 decrease).

Applied on 2026-11-01 by exact single-occurrence substitution with
`notes/referee_external/fixes_2026-11-01/fixer/apply_edits.py` (the script refuses to run unless every `old` text occurs
exactly once in its file; the machine-readable log is `fixer/edits_applied.json` next to it). Analysis and exact checks:
`prover_cells/` (point 1, `geneson_section_check.py`), `prover_AS/` (point 2, `check_AS_endpoint.py`),
`doc_search/` (point 3, `primitive_recount.py`).

## Build

`cd notes/paper && pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper`: exit 0, **0 errors, 0 undefined
references or citations, 0 overfull boxes**. `paper.pdf` is **34 pages** (33 before; the Section 3.1 and 5.5 additions
push the last page over). Cross-references after the rebuild: Section 3.1 = `sec:boundary`, Theorem 5.4 = `thm:AS`,
Section 6.6 = `sec:thrpicture`, Section 7.3 = `sec:finiteeps`, Lemma 7.4 = `lem:Bmech`, Theorem 7.5 = `thm:famE`
(numbering unchanged). Auxiliary files (.aux .bbl .blg .log .out .toc) removed after the build.

Not edited: `notes/paper/paper_draft_before_check.tex` (frozen pre-check snapshot; its l. 388-389 keep the old
sentence as a historical record), `certificates/`, `ANOMALIES.md`, the project instructions file, `n8/`, `lean/`, `region2d.py` (the
implementation never samples, marks or records lower-dimensional cells, so no code change is needed). `REPORT_n6.md`,
`REPORT_n7.md`, `lemma_Bprime_2d.md`, `lemmas.md` and the outreach drafts contain none of the three defects and were left
unchanged.

## The 17 edits

### `notes/paper/paper.tex`

#### Edit 1 (prover_cells edit 1 + doc_search edit 5 (merged))

Reason: Report point 1. The old sentence is false (Geneson's G_1 example, reproduced exactly by prover_cells). Merged wording of prover_cells 5.1 and doc_search: correct description of the lower-dimensional cells, the example with credit, and the statement that every verdict, count and threshold rests on full-dimensional cells via openness.

Old:
```
of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact
rational sample; lower-dimensional cells are boundary points and are never SFP.
```

New:
```
of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact
rational sample; lower-dimensional cells (vertical sections above roots of the projection polynomials, and points on
them) are not sampled. Points on a boundary curve or a positivity root are never SFP, but a section above a root of a
leading coefficient, discriminant or resultant need not lie on any boundary curve and can pass through the open SFP
region: for $G_1$ of Table~\ref{tab:n6} the projection polynomial $\operatorname{Res}_t(\Delta_5,L)$ has a unique root
$\beta_0\in(17/8,15/7)$, the root there of $\Delta_5(\beta,1/(\beta-2))$, and the full support is an SFP at
$(\beta_0,100)$; we thank J.~Geneson for this example.\footnote{\cf{notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py}.}
Since the SFP set $S$ is open (Lemma~\ref{lem:open}) and the full-dimensional cells are dense, they determine it exactly:
$S=\bigl(\overline{S_c}\cap\{L>0\}\bigr)\setminus Z$, where $S_c$ is the union of the SFP full-dimensional cells and $Z$ the
union of the zero sets of the $a_i$, $o_k$, $D$ and $\Delta_{m-1}$; in particular $S\neq\emptyset$ iff some full-dimensional
cell is SFP, and every ``never'' verdict, count and supremum of $q$ below is a statement about full-dimensional
cells.\footnote{\cf{notes/lemma_Bprime_2d.md}, Corollary~B$''$(4)--(5).}
```

#### Edit 2 (prover_cells edit 2)

Reason: Report point 1. The CAD walls also include vertical lines over projection roots, which are not boundary zero sets.

Old:
```
are decided instead on open cells of a decomposition whose walls are the zero sets of a proved boundary set.
```

New:
```
are decided instead on open cells of a decomposition whose walls contain the zero sets of a proved boundary set.
```

#### Edit 3 (prover_cells edit 3)

Reason: Report point 1. beta_1 has the same minimal polynomial as Geneson's beta_0 (unique root in the interval), and the G_1 region includes that vertical section (proved from computed facts in prover_cells/NOTE.md section 1).

Old:
```
$\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$: these graphs are stable at \emph{every} legal $t$ when $\delta/\eps$
```

New:
```
$\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$, the root there of $\operatorname{Res}_t(\Delta_5,L)$ (the
region also contains the whole legal vertical section $\beta=\beta_1$; Section~\ref{sec:boundary}): these graphs are stable at \emph{every} legal $t$ when $\delta/\eps$
```

#### Edit 4 (fixer (supporting doc_search edits 7-8))

Reason: Label added so the corrected sections can be cross-referenced without hard-coded numbers.

Old:
```
\subsection{Boundary sets and three-valued decisions}

```

New:
```
\subsection{Boundary sets and three-valued decisions}\label{sec:boundary}

```

#### Edit 5 (fixer (supporting doc_search edits 7-8))

Reason: Label added so the corrected sections can be cross-referenced without hard-coded numbers.

Old:
```
\subsection{The threshold picture for \texorpdfstring{$n=6,7,8$}{n = 6, 7, 8}}

```

New:
```
\subsection{The threshold picture for \texorpdfstring{$n=6,7,8$}{n = 6, 7, 8}}\label{sec:thrpicture}

```

#### Edit 6 (fixer (supporting doc_search edits 7-8))

Reason: Label added so the corrected sections can be cross-referenced without hard-coded numbers.

Old:
```
\subsection{Finite-\texorpdfstring{$\eps$}{eps}-only counterexamples on every \texorpdfstring{$n\ge7$}{n >= 7}}

```

New:
```
\subsection{Finite-\texorpdfstring{$\eps$}{eps}-only counterexamples on every \texorpdfstring{$n\ge7$}{n >= 7}}\label{sec:finiteeps}

```

#### Edit 7 (prover_AS edit 1)

Reason: Report point 2. The expansion and the boundary curve hold at nu_inf = 1/2; part (a) uses nu_inf > 1/2 and is false at the endpoint when 0 < R* < 1/8, so the endpoint gets its own part (c) with threshold 1/8.

Old:
```
Assume $\nu_\infty:=-\alpha\in(1/2,1)$ and that every root $z_0$ of $\psi$ with $\Rea z_0=\nu_\infty$ is simple and not a root of
$\chi$. Let $R^*=\min\Rea\bigl(1/\hat s\,'(z_0)\bigr)$ over these roots. Then $\nu_A(c)=\nu_\infty+R^*/c+O(c^{-2})$, and the
boundary of the stable region near $(\beta_c,\infty)$ is $\beta=\beta_c+\beta_c^3R^*/t+O(t^{-2})$.
\begin{enumerate}
\item[(a)] If $R^*>0$, then for every $\beta$ slightly above $\beta_c$ the stable $t$-set is nonempty and bounded above; if
moreover $a(\beta_c)>0$, these points are SFPs.
\item[(b)] If $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\beta_c$.
\end{enumerate}
```

New:
```
Assume $\nu_\infty:=-\alpha\in[1/2,1)$, so $\beta_c=1/(1-\nu_\infty)\in[2,\infty)$, and that every root $z_0$ of $\psi$ with
$\Rea z_0=\nu_\infty$ is simple and not a root of $\chi$. Let $R^*=\min\Rea\bigl(1/\hat s\,'(z_0)\bigr)$ over these roots. Then
$\nu_A(c)=\nu_\infty+R^*/c+O(c^{-2})$, and near $(\beta_c,\infty)$ the curve $\nu_A(c)=\gamma$, on which $K$ loses positive
stability, is $\beta=\beta_c+\beta_c^3R^*/t+O(t^{-2})$; for $\nu_\infty>1/2$ it is the boundary of the stable region there.
\begin{enumerate}
\item[(a)] If $\nu_\infty>1/2$ and $R^*>0$, then for every $\beta$ slightly above $\beta_c$ the stable $t$-set is nonempty and
bounded above; if moreover $a(\beta_c)>0$, these points are SFPs.
\item[(b)] If $\nu_\infty>1/2$ and $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\beta_c$.
\item[(c)] If $\nu_\infty=1/2$, then $\beta_c=2$ and the legality boundary $\beta=2+1/t$ is of the same order as the curve
above. If $R^*>1/8$, there are legal stable points with $t\to\infty$ and $\beta\to2^+$, and they are SFPs if $a(2)>0$
(Lemma~\ref{lem:Bmech}); if $R^*<1/8$, no legal point with $\beta$ bounded and $t$ large is stable.
\end{enumerate}
```

#### Edit 8 (prover_AS edit 2)

Reason: Report point 2. The proof sketch now says where the range restriction is used and gives the argument for part (c).

Old:
```
theorem gives an analytic branch $z(e)=z_0+e\,\chi_1(z_0)/\psi_1'(z_0)+O(e^2)$, and $\chi_1(z_0)/\psi_1'(z_0)=1/\hat s\,'(z_0)$.
\cf{notes/lemma_M_refuted.md}, \S4.
```

New:
```
theorem gives an analytic branch $z(e)=z_0+e\,\chi_1(z_0)/\psi_1'(z_0)+O(e^2)$, and $\chi_1(z_0)/\psi_1'(z_0)=1/\hat s\,'(z_0)$.
The other bounded eigenvalues stay near roots of $\psi$ with real part $>\nu_\infty$ and one eigenvalue tends to $+\infty$, so
$\nu_A(c)$ is the minimum of $\Rea z(1/c)$ over the critical branches. None of this uses $\nu_\infty>1/2$; that bound enters
only (a) and (b), where it makes every large $c$ legal for all $\gamma$ near $\nu_\infty$. For (c), legal stability means
$c-\sqrt{c^2-c}<\gamma<\nu_A(c)$, and $c-\sqrt{c^2-c}=1/(1+\sqrt{1-1/c})=\tfrac12+\tfrac1{8c}+O(c^{-2})$ exceeds
$\tfrac12+\tfrac1{8c}$ for every $c>1$. \cf{notes/lemma_M_refuted.md}, \S4.
```

#### Edit 9 (prover_AS edit 3)

Reason: Report point 2. lem:Bmech now cites exactly what it uses (the endpoint expansion) and shows its hypotheses (i), (ii) give nu_inf = 1/2 and R* = R. No circularity: part (c) of thm:AS only points to lem:Bmech.

Old:
```
By Theorem~\ref{thm:AS}, $\nu_A(c)=\tfrac12+R/c+O(c^{-2})$, while the legality boundary is
```

New:
```
By (i) and (ii), $\nu_\infty=1/2$ and the roots of $\psi$ on $\Rea z=1/2$ are exactly $z_0,\bar z_0$, simple and not roots
of $\chi$; at a root of $\psi$, $1/\hat s\,'(z_0)=\chi(z_0)/\psi'(z_0)$, and conjugate roots give the same real part, so
$R^*=R$. The expansion of Theorem~\ref{thm:AS}, which holds at $\nu_\infty=1/2$ (parts (a) and (b) are not used), gives
$\nu_A(c)=\tfrac12+R/c+O(c^{-2})$, while the legality boundary is
```

#### Edit 10 (doc_search edit 1)

Reason: Report point 3. 1 > 34/41 > 6720/9133: the share decreases; only the counts grow. Recounted exactly in notes/referee_external/fixes_2026-11-01/doc_search/primitive_recount.json.

Old:
```
the limit counterexamples for $n=6,7,8$ on one axis. The counts grow by factors of about $41$ and $74$, the share of
primitive counterexamples (containing no smaller one) grows from $3/3$ to $102/123$ to $6{,}720/9{,}133$, and the
threshold grows from $1.1664$ to $1.3247$ to at least $1.7107$, against Geneson's $n-2$.
```

New:
```
the limit counterexamples for $n=6,7,8$ on one axis. The counts grow by factors of about $41$ and $74$; the number of
primitive counterexamples (containing no six- or seven-node counterexample as an induced subgraph) grows from $3$ to $102$ to
$6{,}720$, while their share falls from $3/3=1$ to $102/123=34/41\approx0.83$ to $6{,}720/9{,}133\approx0.74$; and the
threshold grows from $1.1664$ to $1.3247$ to at least $1.7107$, against Geneson's $n-2$.
```

#### Edit 11 (doc_search edit 7 (section numbers replaced by \ref))

Reason: Acknowledge the external correctness report without overstating it (a check of the draft, not of the repository).

Old:
```
The $n=5$ reproduction uses the published package of Curto and collaborators; nauty (geng, directg, labelg) was used for graph enumeration and canonical labelling.
```

New:
```
The $n=5$ reproduction uses the published package of Curto and collaborators; nauty (geng, directg, labelg) was used for graph enumeration and canonical labelling. A correctness report by J.~Geneson on the October 2026 draft led to the corrections in Sections~\ref{sec:boundary}, \ref{sec:thrpicture} and~\ref{sec:finiteeps}.
```

#### Edit 12 (doc_search edit 8 (section numbers replaced by \ref))

Reason: Section 8.2 is where the paper logs defects found by reviewers; records the three external findings there.

Old:
```
by the referees (e.g.\ a misprinted node count, $2p(1)-1$ for $2p(1)-3$, in the notes on Theorem~\ref{thm:U}) are logged and
corrected here.
```

New:
```
by the referees (e.g.\ a misprinted node count, $2p(1)-1$ for $2p(1)-3$, in the notes on Theorem~\ref{thm:U}) are logged and
corrected here. An external correctness report on the October 2026 draft (J.~Geneson) found that Section~\ref{sec:boundary} misdescribed the
lower-dimensional cells of the decomposition, that Lemma~\ref{lem:Bmech} invoked Theorem~\ref{thm:AS} at the endpoint
$\nu_\infty=1/2$ outside its stated range, and that the primitive proportions of Section~\ref{sec:thrpicture} decrease rather than grow; all
three are corrected in this version and none affects a count or a threshold.
```

### `notes/n8_structure.md`

#### Edit 13 (doc_search edit 2)

Reason: Report point 3; this is the source document the paper sentence was copied from.

Old:
```
* The share of primitive counterexamples grows: 3/3, 102/123, 6,720/9,133. Most counterexamples at size n+1
```

New:
```
* The number of primitive counterexamples grows (3, 102, 6,720) while their share falls: 3/3 = 1, 102/123 = 34/41 (about 0.83), 6,720/9,133 (about 0.74). Most counterexamples at size n+1
```

### `notes/FINAL_SUMMARY.md`

#### Edit 14 (doc_search edit 3)

Reason: Report point 3; clarification (the line was not wrong, but reads as if the share were rising).

Old:
```
  (primitive share 3/3, 102/123, 6,720/9,133).
```

New:
```
  (primitive counts 3, 102, 6,720; share 3/3, 102/123, 6,720/9,133, decreasing).
```

### `notes/REPORT_n8.md`

#### Edit 15 (doc_search edit 4)

Reason: Report point 3; clarification.

Old:
```
   whether the growth 3, 123, 9,133 and the primitive share 3/3, 102/123, 6,720/9,133 have a structural
```

New:
```
   whether the growth 3, 123, 9,133 and the decreasing primitive share 3/3, 102/123, 6,720/9,133 have a structural
```

### `notes/lemma_M_refuted.md`

#### Edit 16 (prover_AS edit 4)

Reason: Report point 2; gives the paper's restated thm:AS part (c) a source. Theorem AS in the notes stays correct as stated.

Old:
```
5. The boundary formula follows from `nu_A(c) = gamma`, using `gamma_c - gamma ≈ (beta - beta_c)/beta_c^2` and
   `c ≈ t/beta_c`. ∎
```

New:
```
5. The boundary formula follows from `nu_A(c) = gamma`, using `gamma_c - gamma ≈ (beta - beta_c)/beta_c^2` and
   `c ≈ t/beta_c`. ∎

**Remark (endpoint nu_inf = 1/2; added 2026-11-01).** Step 1 and the boundary formula use only Lemma AS and the
root hypothesis, not nu_inf > 1/2, so `nu_A(c) = 1/2 + R*/c + O(c^-2)` and `beta = 2 + 8 R*/t + O(t^-2)` also hold at
nu_inf = 1/2 (beta_c = 2). The bound nu_inf > 1/2 is used only in step 2. At nu_inf = 1/2, legality
`gamma > nu_P2(c) = 1/(1 + sqrt(1 - 1/c)) = 1/2 + 1/(8c) + O(c^-2)`, with `nu_P2(c) > 1/2 + 1/(8c)` for every c > 1, is of
the same order, so (a) holds with `R* > 1/8` in place of `R* > 0` (this is Lemma B of
notes/research_wf3/theory/notes_theory.md §4.2). If `R* < 1/8`, no legal point with c large is stable; (b) is vacuous.
Exact check on `&FNFL~[udo?` (R* = 1/4): notes/referee_external/fixes_2026-11-01/prover_AS/check_AS_endpoint.py.
```

### `notes/research_wf3/theory/notes_theory.md`

#### Edit 17 (prover_AS edit 5)

Reason: Report point 2; optional clarification so the same conflation does not recur. The note was already correct.

Old:
```
by Lemma AS and (i) those branches are
```

New:
```
by Lemma AS (the local lemma of lemma_M_refuted.md §4, which has no hypothesis on Re z0; the range nu_inf in (1/2, 1) of Theorem AS is not needed) and (i) those branches are
```

## Second pass (referee findings), 2026-11-01

The internal referee (`REFEREE.json`; exact checks in `referee/referee_check.py`, 38/38 pass) confirmed the fixes
above and raised three MINOR wording points and one optional NOTE. Applied by exact single-occurrence substitution with
`notes/referee_external/fixes_2026-11-01/fixer2/apply_edits.py` (aborts unless every `old` text occurs exactly once;
machine-readable log `fixer2/edits_applied.json`). No new numbers are introduced; the primitive-count wording was checked
against the definition used by `doc_search/primitive_recount.py` (n = 7: no six-node counterexample induced, 102; n = 8:
no six- or seven-node limit counterexample induced, 6,720, unchanged when the FNFL graph is included).

### `notes/paper/paper.tex`

#### Edit 18 (referee proposed_edits[0]; MINOR, Section 3.1)

Reason: the parenthetical definition of the lower-dimensional cells omitted the 1-cells that are pieces of the curves
over each open beta-cell, and "on a positivity root" should be "on the vertical line above a positivity root" (a_i
depends only on beta, Lemma 2.3 = Lemma A.2 of the notes).

Old:
```
rational sample; lower-dimensional cells (vertical sections above roots of the projection polynomials, and points on
them) are not sampled. Points on a boundary curve or a positivity root are never SFP, but a section above a root of a
```

New:
```
rational sample; lower-dimensional cells (the pieces of the curves over each open $\beta$-cell, the vertical lines above roots
of the projection polynomials, and the points where these meet) are not sampled. Points on a boundary curve or on the vertical
line above a positivity root are never SFP, but a vertical line above a root of a
```

#### Edit 19 (referee proposed_edits[1]; MINOR, Section 6.6)

Reason: read literally, "no six- or seven-node counterexample as an induced subgraph" makes every seven-node
counterexample non-primitive, since each graph is its own induced subgraph.

Old:
```
(containing no six- or seven-node counterexample as an induced subgraph)
```

New:
```
(containing no counterexample on fewer nodes as an induced subgraph)
```

#### Edit 20 (supporting label for Edit 21)

Old:
```
\subsection{\texorpdfstring{$n=6$}{n = 6}}
```

New:
```
\subsection{\texorpdfstring{$n=6$}{n = 6}}\label{sec:n6map}
```

#### Edit 21 (referee proposed_edits[2]; MINOR, Acknowledgments)

Reason: the list of corrected places omitted Theorem 5.4 (Section 5.5) and the G_1 region paragraph (Section 4.3).
The referee's text cited `sec:n6` (renders "Section 4"); the new label `sec:n6map` renders "Section 4.3", which is
where the paragraph is. Lemma 7.4 is the corrected part of Section 7.3, so `sec:finiteeps` is cited through it.

Old:
```
led to the corrections in Sections~\ref{sec:boundary}, \ref{sec:thrpicture} and~\ref{sec:finiteeps}.
```

New:
```
led to the corrections in Section~\ref{sec:boundary} and the description of the region of $G_1$ (Section~\ref{sec:n6map}), in Theorem~\ref{thm:AS} and Lemma~\ref{lem:Bmech}, and in Section~\ref{sec:thrpicture}.
```

#### Edit 22 (referee NOTE "thm:AS(c) vs lem:Bmech(iv)", optional; applied)

Reason: part (c) of Theorem 5.4 states its SFP condition as a(2) > 0 but cited only Lemma 7.4, whose hypothesis is
w(2) > 0. By Lemma 2.3 (positivity is t-free: a full support is an SFP iff all a_i(beta) > 0 and the characteristic
polynomial is Hurwitz), a(2) > 0 gives a(beta) > 0 for beta near 2 by continuity, so the legal stable points of (c)
are SFPs. Lemma 2.3 is now cited, with Lemma 7.4 kept as a cross-reference.

Old:
```
(Lemma~\ref{lem:Bmech}); if $R^*<1/8$
```

New:
```
(Lemma~\ref{lem:A2}; cf.\ Lemma~\ref{lem:Bmech}); if $R^*<1/8$
```

### Not applied

* Referee NOTE on process ("the fixer ran one read-only git status"): no document change; this pass ran no git.
* Referee remark that J. Geneson is now named three times (Section 3.1, Section 8.2, Acknowledgments) and that his
  permission should be confirmed before the paper circulates: an action for the author, not an edit.
* `referee_build/` (the referee's scratch build copy, including its aux files and `paper.pdf`) was left as the
  referee's record; it is not part of `notes/paper`.

### Build

`cd notes/paper && nice -n 10 pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper`: all four steps exit 0,
**0 errors** (`grep -c '^!' paper.log` = 0), 0 undefined references or citations, 0 overfull boxes, `paper.out` stable
(rerunfilecheck). `paper.pdf` is **34 pages**. Labels after the rebuild: `lem:A2` = Lemma 2.3, `sec:boundary` = 3.1,
`sec:n6map` = 4.3, `thm:AS` = 5.4, `sec:thrpicture` = 6.6, `lem:Bmech` = 7.4. The rendered text contains the three new
phrases and the corrected acknowledgment. Auxiliary files (.aux .bbl .blg .log .out .toc) removed after the build;
`paper.tex` and the `fixer2/` files contain no model or AI-product names.
