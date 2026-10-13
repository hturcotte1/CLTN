# Literature check on Lemmas J.2, J.4, J.4' (2026-10-13, literature agent)

| Result | Verdict |
|---|---|
| J.4' (in-twins inside a support give eigenvalue -eps / +delta; non-adjacent twins never stable) | KNOWN: special case of Curto-Geneson-Morrison, arXiv:1804.00794v2 Lemma 3.13 (simply-added split with uniform-in-degree tau inherits eig(I-W_tau) minus the Perron root) and arXiv:1909.02947v3 Lemma 2.6 / Corollary 2.7 (if tau is not a stable motif, sigma is not). The explicit values -eps/+delta are not displayed there. |
| J.4 (out-edge switching between bidirected in-twins preserves every fixed point of every support, vector, Jacobian charpoly, target-free cliques) | PARTLY KNOWN: ingredients exist (CM 2301.12638 Def 5.5 / Lemma 5.6 simply-embedded; CGM19 Thm 3.4 simply-added gives x_u = x_v; CGM19 Example 6.16 states the switch informally for domination w.r.t. one sigma; Parmelee et al. SIADS 2022 Thm 3.3 / Cor 3.4 is a sibling for removable nodes). The full all-supports statement with Jacobian charpoly and stability was not found. Novelty modest. |
| J.2 (same full-support data iff same charpoly and in-walk counts) | NOT FOUND in any of 15 papers read. Nearest: CM Open Question 4 (hyperplane arrangements), robust motifs (1902.10270), PLOS ONE master graphs with optional edges. |

Stronger prune for Stage 4 / n=7 (published, sound): CGM24 Corollary 2.7 -- any simply-added uniform-in-degree tau inside sigma that is not a stable motif (independent sets, 3-cycles, any tau with in-degree d < |tau|/2) makes sigma unstable.

Paper-side anomaly: CGM24 (arXiv:1909.02947v3) Theorem 2.5 prints eig(I - W_sigma) = {eps, |sigma|(1-eps) - eps} for a clique; the trace forces |sigma|(1-eps) + eps (as in CGM19a p.43). Do not copy the 2024 formula into a test.

Scoop re-check 2026-10-13: no new arXiv items (index reaches ~9 Oct); Geneson's recent papers are extremal combinatorics; Curto's only revision is 2606.00373 (E-I chaos); ipitchford and ChickenNungets repos unchanged since Oct 1 / Aug 8.
Not verified: Neural Computation numbering of CGM19; Parmelee SIADS 2022 beyond an HTML rendering; arXiv 10-13 Oct not indexed yet.
