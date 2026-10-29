# Finalized draft (2026-10-29) -- to be sent by the project owner from his own account; not sent from this repository

To: Jesse Geneson <jesse.geneson@sjsu.edu>  (address as listed in the SJSU Department of Mathematics and Statistics faculty directory; confirm before sending)
Subject: Mapping the gap below q = n - 2 (target-free clique conjecture: exact for n <= 7, eps -> 0 limit at n = 8)
Attachment: paper.pdf (notes/paper/paper.pdf, 33 pages)

Dear Dr. Geneson,

Your counterexample (arXiv:2607.21396) leaves open the gap between q = 1 and your bound
q >= n - 2 - (n - 3) eps / 2. I have mapped it exactly for n <= 7 and in the eps -> 0 limit at n = 8,
using only exact arithmetic (rationals, exact root isolation, Hurwitz determinants) with re-runnable
certificates for every claim. The attached draft has the details; the headline numbers are:

- n <= 5: no counterexamples. n = 6: exactly 3 non-clique graphs among the 1,540,944 digraphs (yours,
  its bidirected in-twin partner, and one more). The critical q is the root of
  797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6), about 1.1664.
- n = 7: exactly 124 among the 882,033,440 digraphs, at every legal (eps, delta); 123 are visible in the
  eps -> 0 limit and one only at finite eps. q_7* is the plastic number (root of q^3 - q - 1, about 1.3247),
  not attained. The n = 6 and n = 7 lists were checked by hand and confirmed by four independent exact
  re-implementations.
- n = 8, eps -> 0 limit only: exactly 9,133 non-complete limit-stable full supports among the
  801,558,826,511 digraphs with at least 29 arcs (the sparser ones are excluded by a proved edge bound);
  the largest limit critical value is about 1.7107, a degree-15 algebraic number. This list was verified by
  two independent machine campaigns but not by hand, and the finite-eps map at n = 8 is not done.
- Every n >= 6 has counterexamples, and the critical q is unbounded in n (proofs in the draft; internally
  reviewed, not yet peer-reviewed).

A word on how this was done: the enumeration, certification and verification pipelines were built and run
by automated software agents under my direction, and every statement rests on the exact certificates in the
repository, which also holds the verifier, the region code, the reports and Lean 4 proofs of the pruning
lemmas: https://github.com/hturcotte1/CLTN

My question: are you mapping this gap yourself, or planning to? If so I would rather coordinate than
duplicate, and I would welcome any correction or comment on the draft. The n = 5 reproduction (zero
mismatches against the published package) and the n = 6 certificates are available as evidence that the
pipeline is sound.

Best regards,
Henry Turcotte
