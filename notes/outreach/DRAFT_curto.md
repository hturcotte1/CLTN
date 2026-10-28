# DRAFT -- NOT SENT -- for the project owner to review

To: Carina Curto
Subject: Your n = 5 package reproduced; exact maps of Open Question 10 for n <= 7, eps -> 0 limit at n = 8

Dear Professor Curto,

We have been working on Open Question 10 (the target-free clique conjecture) since Geneson's
counterexample (arXiv:2607.21396). Your n5-graphs-package was our first validation gate.

- All 9,608 digraphs from n5_digraphs.mat (read as the transpose of the package's sA) give the
  package's fixed-point supports with zero mismatches at eps = 1/4, delta = 1/2 (24,442
  supports, 14,488 stable).
- n <= 5: no counterexample. n = 6: exactly 3 among 1,540,944 digraphs, critical q about
  1.1664. n = 7: exactly 124 among 882,033,440, critical q the plastic number, about 1.3247.
  All exact (Hurwitz determinants), checked by hand (n = 7 not with full confidence) and confirmed by
  independent exact re-implementations.
- Proper supports at n = 7: every stable non-clique proper support is a size-6 support on one of
  the three n = 6 counterexamples; 7,744 of 12,288 (graph, support) pairs are stable somewhere,
  4,544 never; a seventh node that is a target of the support is never stable.
- n = 8, eps -> 0 limit only: 9,133 non-complete limit-stable full supports among the
  801,558,826,511 digraphs with at least 29 arcs. Verified by two independent machine campaigns,
  no human hand check; the finite-eps map at n = 8 is not done.
- Counterexamples exist for every n >= 6; the critical q is unbounded in n (our proofs, internally
  reviewed, not yet peer-reviewed).

The repository holds the verifier, certificates, reports and Lean 4 proofs of the pruning
lemmas. Two questions: is the reproduction above the right reading of the package (parameters and
convention), and would you like the diff file? And are the proper-support results consistent with what
is known about stable supports beyond cliques?

Best regards,
[project owner]

Sources: see README.md here.
