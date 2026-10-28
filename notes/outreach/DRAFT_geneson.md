# DRAFT -- NOT SENT -- for the project owner to review

To: Jesse Geneson
Subject: Mapping the gap below q = n - 2 (target-free clique conjecture: exact for n <= 7, eps -> 0 limit at n = 8)

Dear Dr. Geneson,

Your counterexample (arXiv:2607.21396) leaves open the gap between q = 1 and your bound
q >= n - 2 - (n - 3) eps / 2. We have mapped it exactly for n <= 7 and in the eps -> 0 limit at n = 8
(rationals, exact root isolation,
Hurwitz determinants, re-runnable certificates).

- n <= 5: none. n = 6: exactly 3 non-clique graphs among 1,540,944 digraphs (yours, its
  bidirected in-twin partner, one more). Critical q: the root of
  797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6), about 1.1664.
- n = 7: exactly 124 among 882,033,440 digraphs at every legal (eps, delta); 123 visible in the
  eps -> 0 limit, one only at finite eps. q_7* = the plastic number (root of q^3 - q - 1, about
  1.3247), not attained. n = 6 and 7 were checked by hand (n = 7 not with full confidence) and confirmed
  by four independent exact re-implementations.
- n = 8, eps -> 0 limit only: exactly 9,133 non-complete limit-stable full supports among the
  801,558,826,511 digraphs with at least 29 arcs; largest limit critical value about 1.7107.
  Verified by two independent machine campaigns, no human hand check; the finite-eps map at n = 8 is not done.
- Every n >= 6 has counterexamples, and the critical q is unbounded in n (our proofs, internally
  reviewed, not yet peer-reviewed).

The repository holds the verifier, region code, certificates, reports and Lean 4 proofs of
the pruning lemmas.

My question: are you mapping this gap yourself, or planning to? If so we would rather coordinate
than duplicate; the n = 5 reproduction and n = 6 certificates are available as evidence.

Best regards,
[project owner]

Sources: see README.md here.
