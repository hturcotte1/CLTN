"""Second fixer, 2026-11-01: apply the referee's three MINOR findings (and one optional NOTE) to notes/paper/paper.tex.

Each edit is an exact substitution that must match exactly once in the file; otherwise the script aborts before
writing anything.  A machine-readable log is written next to this script (edits_applied.json).

Run from the repository root:  python3 notes/referee_external/fixes_2026-11-01/fixer2/apply_edits.py
"""
import json, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
PAPER = os.path.join(ROOT, 'notes', 'paper', 'paper.tex')
LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'edits_applied.json')

EDITS = [
    {
        'id': 1, 'finding': 'MINOR, Section 3.1 l. 428-429 (referee proposed_edits[0])',
        'reason': 'The parenthetical definition of the lower-dimensional cells omitted the 1-cells that are pieces of '
                  'the curves over each open beta-cell; and "on a positivity root" should be "on the vertical line '
                  'above a positivity root" (a_i depends only on beta, Lemma A.2).',
        'old': "rational sample; lower-dimensional cells (vertical sections above roots of the projection polynomials, and points on\n"
               "them) are not sampled. Points on a boundary curve or a positivity root are never SFP, but a section above a root of a",
        'new': "rational sample; lower-dimensional cells (the pieces of the curves over each open $\\beta$-cell, the vertical lines above roots\n"
               "of the projection polynomials, and the points where these meet) are not sampled. Points on a boundary curve or on the vertical\n"
               "line above a positivity root are never SFP, but a vertical line above a root of a",
    },
    {
        'id': 2, 'finding': 'MINOR, Section 6.6 l. 1196 (referee proposed_edits[1])',
        'reason': 'Read literally, "no six- or seven-node counterexample as an induced subgraph" makes every seven-node '
                  'counterexample non-primitive (each is its own induced subgraph). The new wording matches what '
                  'doc_search/primitive_recount.py counts: at n=7 no six-node counterexample, at n=8 no six- or '
                  'seven-node limit counterexample, as an induced subgraph (102 and 6,720).',
        'old': "(containing no six- or seven-node counterexample as an induced subgraph)",
        'new': "(containing no counterexample on fewer nodes as an induced subgraph)",
    },
    {
        'id': 3, 'finding': 'supporting edit for 4 (label)',
        'reason': 'Label on the n=6 subsection (4.3), which holds the corrected G_1 region paragraph, so the '
                  'acknowledgment can cite it without a hard-coded number.',
        'old': "\\subsection{\\texorpdfstring{$n=6$}{n = 6}}\n",
        'new': "\\subsection{\\texorpdfstring{$n=6$}{n = 6}}\\label{sec:n6map}\n",
    },
    {
        'id': 4, 'finding': 'MINOR, Acknowledgments l. 1618 (referee proposed_edits[2], with the G_1 reference pointed '
                            'at the subsection 4.3 rather than Section 4)',
        'reason': 'The list of corrected places omitted Theorem 5.4 (thm:AS, Section 5.5) and the G_1 region '
                  'paragraph (Section 4.3). Lemma 7.4 (lem:Bmech) is the corrected part of Section 7.3.',
        'old': "led to the corrections in Sections~\\ref{sec:boundary}, \\ref{sec:thrpicture} and~\\ref{sec:finiteeps}.",
        'new': "led to the corrections in Section~\\ref{sec:boundary} and the description of the region of $G_1$ "
               "(Section~\\ref{sec:n6map}), in Theorem~\\ref{thm:AS} and Lemma~\\ref{lem:Bmech}, and in Section~\\ref{sec:thrpicture}.",
    },
    {
        'id': 5, 'finding': 'NOTE (optional), thm:AS(c) vs lem:Bmech(iv)',
        'reason': 'Part (c) states its SFP condition as a(2) > 0 but cites lem:Bmech, whose hypothesis is w(2) > 0. '
                  'The condition a(2) > 0 is sufficient by Lemma A.2 (full-support SFP iff all a_i(beta) > 0 and '
                  'Hurwitz; a(2) > 0 gives a(beta) > 0 for beta near 2 by continuity), so cite Lemma A.2 and keep '
                  'lem:Bmech as a cross-reference.',
        'old': "(Lemma~\\ref{lem:Bmech}); if $R^*<1/8$",
        'new': "(Lemma~\\ref{lem:A2}; cf.\\ Lemma~\\ref{lem:Bmech}); if $R^*<1/8$",
    },
]

def main():
    with open(PAPER, encoding='utf-8') as f:
        text = f.read()
    for e in EDITS:
        c = text.count(e['old'])
        if c != 1:
            sys.exit(f"edit {e['id']}: old text occurs {c} times (need exactly 1); nothing written")
        if e['new'] in text:
            sys.exit(f"edit {e['id']}: new text already present; nothing written")
    for e in EDITS:
        text = text.replace(e['old'], e['new'], 1)
    for e in EDITS:
        assert text.count(e['new']) == 1 and e['old'] not in text
    with open(PAPER, 'w', encoding='utf-8') as f:
        f.write(text)
    with open(LOG, 'w', encoding='utf-8') as f:
        json.dump({'file': 'notes/paper/paper.tex', 'date': '2026-11-01', 'edits': EDITS}, f, indent=1)
    print(f"applied {len(EDITS)} edits to {PAPER}; log at {LOG}")

if __name__ == '__main__':
    main()
