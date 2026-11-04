# Edits applied for the 4 November 2026 correctness report (J. Geneson)

Point 1 (very minor): in Section 3.1, the paragraph "The two-parameter version" listed the
boundary curves a_i = 0, D = 0, Delta_{m-1} = 0, L = 0 but omitted the off-support curves
o_k = 0 (k not in sigma) for proper supports. Confirmed exactly by the prover
(notes/referee_external/fixes_2026-11-04/prover/NOTE.md); the implementation (fastregion,
region2d) and the certificates already include the o_k, so no count or threshold changes.

Script: notes/referee_external/fixes_2026-11-04/fixer/apply_edits.py (exact single-occurrence
substitutions on notes/paper/paper.tex). Edits:

1. Section 3.1 boundary list, lines ~422-426:
   old: "the complement of the curves $\{a_i=0\}$, $\{D=0\}$, $\{\Delta_{m-1}=0\}$ and the legality boundary
         $L(\beta,t)=t(\beta-2)-1=0$ in the open legal region, and that"
   new: "the complement of the curves $\{a_i=0\}$, $\{o_k=0\}$ ($k\notin\sigma$), $\{D=0\}$, $\{\Delta_{m-1}=0\}$ and the
         legality boundary $L(\beta,t)=t(\beta-2)-1=0$ in the open legal region (for $\sigma=[n]$ there are no $o_k$ and the
         list reduces to the $a_i$, $D$, $\Delta_{m-1}$ and $L$), and that"

2. Section 3.1, line ~431 (optional consistency edit; "positivity root" was undefined):
   old: "line above a positivity root are never SFP, but a vertical line above a root of a"
   new: "line above a root of some $a_i$ or $o_k$ are never SFP, but a vertical line above a root of a"

3. Section 8.2 defect log, after "...none affects a count or a threshold.":
   added: "A second report (4 November) found that the boundary list for proper supports in
   Section~\ref{sec:boundary} omitted the off-support curves $o_k=0$; this is corrected here and
   affects no count or threshold."

4. Acknowledgments, end of the Geneson sentence:
   added: "; a second report (4 November) led to the correction of the boundary list for proper
   supports in Section~\ref{sec:boundary}."

Not edited: notes/paper/paper_draft_before_check.tex (archived draft, same omission, left alone).

Rebuild: pdflatex + bibtex + pdflatex x2 in notes/paper, nice -n 10. Zero LaTeX errors, zero
undefined references, no rerun warnings. Output: paper.pdf, 34 pages (unchanged page count).
Removed .aux .bbl .blg .log .out .toc.
