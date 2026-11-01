#!/usr/bin/env python3
"""Apply the 2026-11-01 external-report fixes by exact single-occurrence substitution.

Run from the repository root:  python3 notes/referee_external/fixes_2026-11-01/fixer/apply_edits.py [--check]

Each edit is (file, old, new, reason).  The script refuses to run unless every `old` occurs exactly once in
its file (so it is idempotent in the sense that a second run fails loudly instead of double-applying).
With --check it only verifies the counts and exits.  It also writes the edit log used for APPLIED.md to
notes/referee_external/fixes_2026-11-01/fixer/edits_applied.json.
"""
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))

EDITS = []


def edit(file, old, new, reason, source):
    EDITS.append({"file": file, "old": old, "new": new, "reason": reason, "source": source})


# ----------------------------------------------------------------------------------------------------------
# paper.tex, point 1 (Section 3.1): merged wording of prover_cells (5.1) and doc_search.
# ----------------------------------------------------------------------------------------------------------
edit(
    "notes/paper/paper.tex",
    "of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact\n"
    "rational sample; lower-dimensional cells are boundary points and are never SFP.",
    "of a cylindrical algebraic decomposition are sign-invariant cells. Each full-dimensional cell is decided at one exact\n"
    "rational sample; lower-dimensional cells (vertical sections above roots of the projection polynomials, and points on\n"
    "them) are not sampled. Points on a boundary curve or a positivity root are never SFP, but a section above a root of a\n"
    "leading coefficient, discriminant or resultant need not lie on any boundary curve and can pass through the open SFP\n"
    "region: for $G_1$ of Table~\\ref{tab:n6} the projection polynomial $\\operatorname{Res}_t(\\Delta_5,L)$ has a unique root\n"
    "$\\beta_0\\in(17/8,15/7)$, the root there of $\\Delta_5(\\beta,1/(\\beta-2))$, and the full support is an SFP at\n"
    "$(\\beta_0,100)$; we thank J.~Geneson for this example.\\footnote{\\cf{notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py}.}\n"
    "Since the SFP set $S$ is open (Lemma~\\ref{lem:open}) and the full-dimensional cells are dense, they determine it exactly:\n"
    "$S=\\bigl(\\overline{S_c}\\cap\\{L>0\\}\\bigr)\\setminus Z$, where $S_c$ is the union of the SFP full-dimensional cells and $Z$ the\n"
    "union of the zero sets of the $a_i$, $o_k$, $D$ and $\\Delta_{m-1}$; in particular $S\\neq\\emptyset$ iff some full-dimensional\n"
    "cell is SFP, and every ``never'' verdict, count and supremum of $q$ below is a statement about full-dimensional\n"
    "cells.\\footnote{\\cf{notes/lemma_Bprime_2d.md}, Corollary~B$''$(4)--(5).}",
    "Report point 1. The old sentence is false (Geneson's G_1 example, reproduced exactly by prover_cells). "
    "Merged wording of prover_cells 5.1 and doc_search: correct description of the lower-dimensional cells, the example "
    "with credit, and the statement that every verdict, count and threshold rests on full-dimensional cells via openness.",
    "prover_cells edit 1 + doc_search edit 5 (merged)",
)

edit(
    "notes/paper/paper.tex",
    "are decided instead on open cells of a decomposition whose walls are the zero sets of a proved boundary set.",
    "are decided instead on open cells of a decomposition whose walls contain the zero sets of a proved boundary set.",
    "Report point 1. The CAD walls also include vertical lines over projection roots, which are not boundary zero sets.",
    "prover_cells edit 2",
)

edit(
    "notes/paper/paper.tex",
    "$\\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$: these graphs are stable at \\emph{every} legal $t$ when $\\delta/\\eps$",
    "$\\beta_1$ a degree-$14$ algebraic number in $(17/8,15/7)$, the root there of $\\operatorname{Res}_t(\\Delta_5,L)$ (the\n"
    "region also contains the whole legal vertical section $\\beta=\\beta_1$; Section~\\ref{sec:boundary}): these graphs are stable at \\emph{every} legal $t$ when $\\delta/\\eps$",
    "Report point 1. beta_1 has the same minimal polynomial as Geneson's beta_0 (unique root in the interval), and the G_1 "
    "region includes that vertical section (proved from computed facts in prover_cells/NOTE.md section 1).",
    "prover_cells edit 3",
)

# subsection labels so that the acknowledgment and audit sentences use \ref instead of hard-coded numbers
edit(
    "notes/paper/paper.tex",
    "\\subsection{Boundary sets and three-valued decisions}\n",
    "\\subsection{Boundary sets and three-valued decisions}\\label{sec:boundary}\n",
    "Label added so the corrected sections can be cross-referenced without hard-coded numbers.",
    "fixer (supporting doc_search edits 7-8)",
)
edit(
    "notes/paper/paper.tex",
    "\\subsection{The threshold picture for \\texorpdfstring{$n=6,7,8$}{n = 6, 7, 8}}\n",
    "\\subsection{The threshold picture for \\texorpdfstring{$n=6,7,8$}{n = 6, 7, 8}}\\label{sec:thrpicture}\n",
    "Label added so the corrected sections can be cross-referenced without hard-coded numbers.",
    "fixer (supporting doc_search edits 7-8)",
)
edit(
    "notes/paper/paper.tex",
    "\\subsection{Finite-\\texorpdfstring{$\\eps$}{eps}-only counterexamples on every \\texorpdfstring{$n\\ge7$}{n >= 7}}\n",
    "\\subsection{Finite-\\texorpdfstring{$\\eps$}{eps}-only counterexamples on every \\texorpdfstring{$n\\ge7$}{n >= 7}}\\label{sec:finiteeps}\n",
    "Label added so the corrected sections can be cross-referenced without hard-coded numbers.",
    "fixer (supporting doc_search edits 7-8)",
)

# ----------------------------------------------------------------------------------------------------------
# paper.tex, point 2 (thm:AS / lem:Bmech), per prover_AS.
# ----------------------------------------------------------------------------------------------------------
edit(
    "notes/paper/paper.tex",
    "Assume $\\nu_\\infty:=-\\alpha\\in(1/2,1)$ and that every root $z_0$ of $\\psi$ with $\\Rea z_0=\\nu_\\infty$ is simple and not a root of\n"
    "$\\chi$. Let $R^*=\\min\\Rea\\bigl(1/\\hat s\\,'(z_0)\\bigr)$ over these roots. Then $\\nu_A(c)=\\nu_\\infty+R^*/c+O(c^{-2})$, and the\n"
    "boundary of the stable region near $(\\beta_c,\\infty)$ is $\\beta=\\beta_c+\\beta_c^3R^*/t+O(t^{-2})$.\n"
    "\\begin{enumerate}\n"
    "\\item[(a)] If $R^*>0$, then for every $\\beta$ slightly above $\\beta_c$ the stable $t$-set is nonempty and bounded above; if\n"
    "moreover $a(\\beta_c)>0$, these points are SFPs.\n"
    "\\item[(b)] If $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\\beta_c$.\n"
    "\\end{enumerate}",
    "Assume $\\nu_\\infty:=-\\alpha\\in[1/2,1)$, so $\\beta_c=1/(1-\\nu_\\infty)\\in[2,\\infty)$, and that every root $z_0$ of $\\psi$ with\n"
    "$\\Rea z_0=\\nu_\\infty$ is simple and not a root of $\\chi$. Let $R^*=\\min\\Rea\\bigl(1/\\hat s\\,'(z_0)\\bigr)$ over these roots. Then\n"
    "$\\nu_A(c)=\\nu_\\infty+R^*/c+O(c^{-2})$, and near $(\\beta_c,\\infty)$ the curve $\\nu_A(c)=\\gamma$, on which $K$ loses positive\n"
    "stability, is $\\beta=\\beta_c+\\beta_c^3R^*/t+O(t^{-2})$; for $\\nu_\\infty>1/2$ it is the boundary of the stable region there.\n"
    "\\begin{enumerate}\n"
    "\\item[(a)] If $\\nu_\\infty>1/2$ and $R^*>0$, then for every $\\beta$ slightly above $\\beta_c$ the stable $t$-set is nonempty and\n"
    "bounded above; if moreover $a(\\beta_c)>0$, these points are SFPs.\n"
    "\\item[(b)] If $\\nu_\\infty>1/2$ and $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\\beta_c$.\n"
    "\\item[(c)] If $\\nu_\\infty=1/2$, then $\\beta_c=2$ and the legality boundary $\\beta=2+1/t$ is of the same order as the curve\n"
    "above. If $R^*>1/8$, there are legal stable points with $t\\to\\infty$ and $\\beta\\to2^+$, and they are SFPs if $a(2)>0$\n"
    "(Lemma~\\ref{lem:Bmech}); if $R^*<1/8$, no legal point with $\\beta$ bounded and $t$ large is stable.\n"
    "\\end{enumerate}",
    "Report point 2. The expansion and the boundary curve hold at nu_inf = 1/2; part (a) uses nu_inf > 1/2 and is false at the "
    "endpoint when 0 < R* < 1/8, so the endpoint gets its own part (c) with threshold 1/8.",
    "prover_AS edit 1",
)

edit(
    "notes/paper/paper.tex",
    "theorem gives an analytic branch $z(e)=z_0+e\\,\\chi_1(z_0)/\\psi_1'(z_0)+O(e^2)$, and $\\chi_1(z_0)/\\psi_1'(z_0)=1/\\hat s\\,'(z_0)$.\n"
    "\\cf{notes/lemma_M_refuted.md}, \\S4.",
    "theorem gives an analytic branch $z(e)=z_0+e\\,\\chi_1(z_0)/\\psi_1'(z_0)+O(e^2)$, and $\\chi_1(z_0)/\\psi_1'(z_0)=1/\\hat s\\,'(z_0)$.\n"
    "The other bounded eigenvalues stay near roots of $\\psi$ with real part $>\\nu_\\infty$ and one eigenvalue tends to $+\\infty$, so\n"
    "$\\nu_A(c)$ is the minimum of $\\Rea z(1/c)$ over the critical branches. None of this uses $\\nu_\\infty>1/2$; that bound enters\n"
    "only (a) and (b), where it makes every large $c$ legal for all $\\gamma$ near $\\nu_\\infty$. For (c), legal stability means\n"
    "$c-\\sqrt{c^2-c}<\\gamma<\\nu_A(c)$, and $c-\\sqrt{c^2-c}=1/(1+\\sqrt{1-1/c})=\\tfrac12+\\tfrac1{8c}+O(c^{-2})$ exceeds\n"
    "$\\tfrac12+\\tfrac1{8c}$ for every $c>1$. \\cf{notes/lemma_M_refuted.md}, \\S4.",
    "Report point 2. The proof sketch now says where the range restriction is used and gives the argument for part (c).",
    "prover_AS edit 2",
)

edit(
    "notes/paper/paper.tex",
    "By Theorem~\\ref{thm:AS}, $\\nu_A(c)=\\tfrac12+R/c+O(c^{-2})$, while the legality boundary is",
    "By (i) and (ii), $\\nu_\\infty=1/2$ and the roots of $\\psi$ on $\\Rea z=1/2$ are exactly $z_0,\\bar z_0$, simple and not roots\n"
    "of $\\chi$; at a root of $\\psi$, $1/\\hat s\\,'(z_0)=\\chi(z_0)/\\psi'(z_0)$, and conjugate roots give the same real part, so\n"
    "$R^*=R$. The expansion of Theorem~\\ref{thm:AS}, which holds at $\\nu_\\infty=1/2$ (parts (a) and (b) are not used), gives\n"
    "$\\nu_A(c)=\\tfrac12+R/c+O(c^{-2})$, while the legality boundary is",
    "Report point 2. lem:Bmech now cites exactly what it uses (the endpoint expansion) and shows its hypotheses (i), (ii) give "
    "nu_inf = 1/2 and R* = R. No circularity: part (c) of thm:AS only points to lem:Bmech.",
    "prover_AS edit 3",
)

# ----------------------------------------------------------------------------------------------------------
# paper.tex, point 3 (Section 6.6), per doc_search.
# ----------------------------------------------------------------------------------------------------------
edit(
    "notes/paper/paper.tex",
    "the limit counterexamples for $n=6,7,8$ on one axis. The counts grow by factors of about $41$ and $74$, the share of\n"
    "primitive counterexamples (containing no smaller one) grows from $3/3$ to $102/123$ to $6{,}720/9{,}133$, and the\n"
    "threshold grows from $1.1664$ to $1.3247$ to at least $1.7107$, against Geneson's $n-2$.",
    "the limit counterexamples for $n=6,7,8$ on one axis. The counts grow by factors of about $41$ and $74$; the number of\n"
    "primitive counterexamples (containing no six- or seven-node counterexample as an induced subgraph) grows from $3$ to $102$ to\n"
    "$6{,}720$, while their share falls from $3/3=1$ to $102/123=34/41\\approx0.83$ to $6{,}720/9{,}133\\approx0.74$; and the\n"
    "threshold grows from $1.1664$ to $1.3247$ to at least $1.7107$, against Geneson's $n-2$.",
    "Report point 3. 1 > 34/41 > 6720/9133: the share decreases; only the counts grow. Recounted exactly in "
    "notes/referee_external/fixes_2026-11-01/doc_search/primitive_recount.json.",
    "doc_search edit 1",
)

# ----------------------------------------------------------------------------------------------------------
# paper.tex, acknowledgment of the external report (Acknowledgments + Section 8.2), per doc_search.
# ----------------------------------------------------------------------------------------------------------
edit(
    "notes/paper/paper.tex",
    "The $n=5$ reproduction uses the published package of Curto and collaborators; nauty (geng, directg, labelg) was used for graph enumeration and canonical labelling.",
    "The $n=5$ reproduction uses the published package of Curto and collaborators; nauty (geng, directg, labelg) was used for graph enumeration and canonical labelling. A correctness report by J.~Geneson on the October 2026 draft led to the corrections in Sections~\\ref{sec:boundary}, \\ref{sec:thrpicture} and~\\ref{sec:finiteeps}.",
    "Acknowledge the external correctness report without overstating it (a check of the draft, not of the repository).",
    "doc_search edit 7 (section numbers replaced by \\ref)",
)

edit(
    "notes/paper/paper.tex",
    "by the referees (e.g.\\ a misprinted node count, $2p(1)-1$ for $2p(1)-3$, in the notes on Theorem~\\ref{thm:U}) are logged and\n"
    "corrected here.",
    "by the referees (e.g.\\ a misprinted node count, $2p(1)-1$ for $2p(1)-3$, in the notes on Theorem~\\ref{thm:U}) are logged and\n"
    "corrected here. An external correctness report on the October 2026 draft (J.~Geneson) found that Section~\\ref{sec:boundary} misdescribed the\n"
    "lower-dimensional cells of the decomposition, that Lemma~\\ref{lem:Bmech} invoked Theorem~\\ref{thm:AS} at the endpoint\n"
    "$\\nu_\\infty=1/2$ outside its stated range, and that the primitive proportions of Section~\\ref{sec:thrpicture} decrease rather than grow; all\n"
    "three are corrected in this version and none affects a count or a threshold.",
    "Section 8.2 is where the paper logs defects found by reviewers; records the three external findings there.",
    "doc_search edit 8 (section numbers replaced by \\ref)",
)

# ----------------------------------------------------------------------------------------------------------
# Notes files.
# ----------------------------------------------------------------------------------------------------------
edit(
    "notes/n8_structure.md",
    "* The share of primitive counterexamples grows: 3/3, 102/123, 6,720/9,133. Most counterexamples at size n+1",
    "* The number of primitive counterexamples grows (3, 102, 6,720) while their share falls: 3/3 = 1, 102/123 = 34/41 (about 0.83), 6,720/9,133 (about 0.74). Most counterexamples at size n+1",
    "Report point 3; this is the source document the paper sentence was copied from.",
    "doc_search edit 2",
)
edit(
    "notes/FINAL_SUMMARY.md",
    "  (primitive share 3/3, 102/123, 6,720/9,133).",
    "  (primitive counts 3, 102, 6,720; share 3/3, 102/123, 6,720/9,133, decreasing).",
    "Report point 3; clarification (the line was not wrong, but reads as if the share were rising).",
    "doc_search edit 3",
)
edit(
    "notes/REPORT_n8.md",
    "   whether the growth 3, 123, 9,133 and the primitive share 3/3, 102/123, 6,720/9,133 have a structural",
    "   whether the growth 3, 123, 9,133 and the decreasing primitive share 3/3, 102/123, 6,720/9,133 have a structural",
    "Report point 3; clarification.",
    "doc_search edit 4",
)
edit(
    "notes/lemma_M_refuted.md",
    "5. The boundary formula follows from `nu_A(c) = gamma`, using `gamma_c - gamma ≈ (beta - beta_c)/beta_c^2` and\n"
    "   `c ≈ t/beta_c`. ∎",
    "5. The boundary formula follows from `nu_A(c) = gamma`, using `gamma_c - gamma ≈ (beta - beta_c)/beta_c^2` and\n"
    "   `c ≈ t/beta_c`. ∎\n"
    "\n"
    "**Remark (endpoint nu_inf = 1/2; added 2026-11-01).** Step 1 and the boundary formula use only Lemma AS and the\n"
    "root hypothesis, not nu_inf > 1/2, so `nu_A(c) = 1/2 + R*/c + O(c^-2)` and `beta = 2 + 8 R*/t + O(t^-2)` also hold at\n"
    "nu_inf = 1/2 (beta_c = 2). The bound nu_inf > 1/2 is used only in step 2. At nu_inf = 1/2, legality\n"
    "`gamma > nu_P2(c) = 1/(1 + sqrt(1 - 1/c)) = 1/2 + 1/(8c) + O(c^-2)`, with `nu_P2(c) > 1/2 + 1/(8c)` for every c > 1, is of\n"
    "the same order, so (a) holds with `R* > 1/8` in place of `R* > 0` (this is Lemma B of\n"
    "notes/research_wf3/theory/notes_theory.md §4.2). If `R* < 1/8`, no legal point with c large is stable; (b) is vacuous.\n"
    "Exact check on `&FNFL~[udo?` (R* = 1/4): notes/referee_external/fixes_2026-11-01/prover_AS/check_AS_endpoint.py.",
    "Report point 2; gives the paper's restated thm:AS part (c) a source. Theorem AS in the notes stays correct as stated.",
    "prover_AS edit 4",
)
edit(
    "notes/research_wf3/theory/notes_theory.md",
    "by Lemma AS and (i) those branches are",
    "by Lemma AS (the local lemma of lemma_M_refuted.md §4, which has no hypothesis on Re z0; the range nu_inf in (1/2, 1) of Theorem AS is not needed) and (i) those branches are",
    "Report point 2; optional clarification so the same conflation does not recur. The note was already correct.",
    "prover_AS edit 5",
)


def main():
    check_only = "--check" in sys.argv
    problems = []
    contents = {}
    for e in EDITS:
        path = os.path.join(ROOT, e["file"])
        if path not in contents:
            with open(path, encoding="utf-8") as fh:
                contents[path] = fh.read()
        n = contents[path].count(e["old"])
        if n != 1:
            problems.append(f"{e['file']}: old text occurs {n} times (need exactly 1): {e['old'][:70]!r}")
        else:
            contents[path] = contents[path].replace(e["old"], e["new"], 1)
    if problems:
        print("REFUSING TO APPLY:")
        for p in problems:
            print("  " + p)
        return 1
    print(f"all {len(EDITS)} edits match exactly once")
    if check_only:
        return 0
    for path, text in contents.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        print("wrote", os.path.relpath(path, ROOT))
    with open(os.path.join(HERE, "edits_applied.json"), "w", encoding="utf-8") as fh:
        json.dump(EDITS, fh, indent=1, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
