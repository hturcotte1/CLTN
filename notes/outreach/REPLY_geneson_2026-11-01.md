# Reply draft (2026-11-01) -- to be sent by the author; not sent from this repository

To: Jesse Geneson <jesse.geneson@sjsu.edu>
Subject: Re: The gap below your target-free clique bound
Attachment: paper.pdf (revised, 34 pages)

Dear Jesse,

Thank you for the report; all three points are right, and I have corrected the draft (attached).

1. Section 3.1: you are right that a vertical section above a projection root can cross the open stable region.
   I reproduced your example exactly: Res_t(Delta_5, L) for G_1 has exactly one root beta_0 in (17/8, 15/7), on
   the same degree-14 factor as the certificate's critical value, and at (beta_0, 100) all a_i and all six Hurwitz
   determinants are positive. In fact the whole legal line beta = beta_0 is stable. The engine only ever samples
   full-dimensional cells and records neither sections nor curves, so, as you say, the verdicts, counts and
   thresholds are unaffected; the text now describes the regions as S = (closure of the stable open cells,
   intersected with the legal region) minus the zero set of the decision polynomials, with your example cited.

2. Theorem 5.4 / Lemma 7.4: the expansion only needs a simple root off the fixed spectrum, so the theorem is now
   stated for nu_inf in [1/2, 1) for the expansion and the boundary curve, with parts (a) and (b) kept for
   nu_inf > 1/2 and a new part (c) at the endpoint with the threshold R* > 1/8 (part (a) is false there for
   0 < R* < 1/8, because the legality boundary is 1/2 + 1/(8c) + O(c^-2)). Lemma 7.4 now cites only that.

3. Section 6.6: the counts grow (3, 102, 6,720) and the share falls (1, 34/41, 6720/9133); fixed.

I have credited the report in Section 3.1, Section 8.2 and the acknowledgments. Please tell me if you would rather
not be named, or prefer different wording. I would also be glad of any further comments before I post the preprint.

Best regards,
Henry Turcotte
