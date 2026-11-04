"""Apply the 2026-11-04 edits to notes/paper/paper.tex by exact single-occurrence substitution."""
import sys
P = "/home/user/CLTN/notes/paper/paper.tex"
EDITS = [
 ("the complement of the curves $\\{a_i=0\\}$, $\\{D=0\\}$, $\\{\\Delta_{m-1}=0\\}$ and the legality boundary\n$L(\\beta,t)=t(\\beta-2)-1=0$ in the open legal region, and that",
  "the complement of the curves $\\{a_i=0\\}$, $\\{o_k=0\\}$ ($k\\notin\\sigma$), $\\{D=0\\}$, $\\{\\Delta_{m-1}=0\\}$ and the\nlegality boundary $L(\\beta,t)=t(\\beta-2)-1=0$ in the open legal region (for $\\sigma=[n]$ there are no $o_k$ and the\nlist reduces to the $a_i$, $D$, $\\Delta_{m-1}$ and $L$), and that"),
 ("line above a positivity root are never SFP, but a vertical line above a root of a",
  "line above a root of some $a_i$ or $o_k$ are never SFP, but a vertical line above a root of a"),
 ("three are corrected in this version and none affects a count or a threshold.",
  "three are corrected in this version and none affects a count or a threshold. A second report (4 November) found that the boundary list\nfor proper supports in Section~\\ref{sec:boundary} omitted the off-support curves $o_k=0$; this is corrected here and affects no count or threshold."),
 ("in Theorem~\\ref{thm:AS} and Lemma~\\ref{lem:Bmech}, and in Section~\\ref{sec:thrpicture}.",
  "in Theorem~\\ref{thm:AS} and Lemma~\\ref{lem:Bmech}, and in Section~\\ref{sec:thrpicture}; a second report (4 November) led to the correction of the boundary list for proper supports in Section~\\ref{sec:boundary}."),
]
s = open(P).read()
for old, new in EDITS:
    c = s.count(old)
    if c != 1:
        print("FAIL count", c, "for", old[:60]); sys.exit(1)
    s = s.replace(old, new)
open(P, "w").write(s)
print("applied", len(EDITS), "edits")
