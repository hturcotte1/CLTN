# Finalized draft (2026-10-29, short form) -- to be sent by the project owner from his own account; not sent from this repository

To: Jesse Geneson <jesse.geneson@sjsu.edu>  (SJSU Department of Mathematics and Statistics directory; confirm before sending)
Subject: The gap below your target-free clique bound: exact map for n <= 7, limit map at n = 8
Attachment: paper.pdf

Dear Dr. Geneson,

Following your counterexample (arXiv:2607.21396), I have mapped the gap between q = 1 and your bound
q >= n - 2 - (n - 3) eps / 2, with exact arithmetic and re-runnable certificates (draft attached):

- n = 6: exactly 3 counterexamples among all 1,540,944 digraphs (yours, its in-twin partner, one more);
  critical q about 1.1664, an explicit sextic root.
- n = 7: exactly 124 among 882,033,440, at every legal (eps, delta); critical q = the plastic number.
- n = 8, eps -> 0 limit only: exactly 9,133 among the 801,558,826,511 digraphs with >= 29 arcs; largest
  limit critical value about 1.7107 (degree 15). Machine-verified twice, not hand-checked.
- Counterexamples exist for every n >= 6, with critical q unbounded in n.

The pipelines were built and run by automated software agents under my direction; every claim rests on the
certificates in https://github.com/hturcotte1/CLTN (verifier, certificates, reports, Lean proofs of the lemmas).

Are you mapping this gap yourself? If so I would rather coordinate than duplicate, and I would welcome any
comment on the draft.

Best regards,
Henry Turcotte
