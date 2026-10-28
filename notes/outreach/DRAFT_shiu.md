# DRAFT -- NOT SENT -- for the project owner to review

To: Anne Shiu
Subject: Exact maps of the target-free clique conjecture, n <= 7 (and the eps -> 0 limit at n = 8)

Dear Professor Shiu,

After Geneson's counterexample to the target-free clique conjecture (arXiv:2607.21396) we
computed exactly where it holds and fails for small n (exact arithmetic, Hurwitz determinants,
re-runnable certificates).

- The published n = 5 classification (ccurto/n5-graphs-package) is reproduced with zero
  mismatches on all 9,608 digraphs (24,442 fixed-point supports, 14,488 stable) at eps = 1/4,
  delta = 1/2.
- n <= 5: no counterexample. n = 6: exactly 3 non-clique stable full supports, critical q
  about 1.1664. n = 7: exactly 124, critical q the plastic number, about 1.3247. Both
  checked by hand (n = 7 not with full confidence) and confirmed by independent exact re-implementations.
- Proper supports at n = 7: every stable non-clique proper support is a size-6 support on one of
  the three n = 6 counterexamples plus a seventh node; 7,744 of the 12,288 (graph, support) pairs
  are stable somewhere, 4,544 never. A seventh node that is a target of the support is never
  stable; one receiving no edge from it gives a support that is stable at some legal parameter in all 192
  such cases (certified, not hand-verified).
- n = 8, eps -> 0 limit only: 9,133 non-complete limit-stable full supports among the
  801,558,826,511 digraphs with at least 29 arcs. Verified by two independent machine campaigns,
  no human hand check; the finite-eps map at n = 8 is not done.
- Counterexamples exist for every n >= 6; the critical q is unbounded in n (our proofs, internally
  reviewed, not yet peer-reviewed).

The repository holds the verifier, certificates, reports and Lean 4 proofs of the pruning
lemmas. Two questions: does the n = 5 reproduction at eps = 1/4, delta = 1/2 match the
package's own parameters? And is the proper-support classification of interest, or already
known in some form?

Best regards,
[project owner]

Sources: see README.md here.
