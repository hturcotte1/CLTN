# Certificate index

Every result in this project exists only as a re-runnable certificate file. This index says what
each file certifies and which script regenerates it. Conventions: CLAUDE.md section 2. All values
exact (integers, rationals, algebraic numbers as minimal polynomial + isolating interval); decimals
are labelled approximations only.

## Validation gate
| file | certifies | regenerate |
|---|---|---|
| `stage0_geneson.json` | Geneson's 6-node counterexample (arXiv:2607.21396) reproduced: stable positive non-clique full support at eps = 1/100, delta = 29/2500 | `python3 test_ctln.py` (block 8) |
| `stage1_n5_diff.json`, `../data/n5_package_FP.csv` | the published n = 5 classification (ccurto/n5-graphs-package) reproduced with diff zero | `stage1_n5_diff.py` |
| `validate_fastregion_n4.json`, `_n5.json` | the fast exact region engine agrees with `symbolic.q_region` | `validate_fastregion.py` |
| `region2d_n2to4_all.json`, `region2d_v2_validation.json` | the two-parameter engine on all digraphs with n <= 4; v1/v2 agreement | `region2d.py` |

## n = 5 and n = 6 (final, audited; n = 6 hits hand-verified 2026-10-14)
| file | certifies |
|---|---|
| `digraphs_n5.keys`, `digraphs_n6.keys` | the enumerations (A000273: 9,608 and 1,540,944) |
| `region_n5_limit_*`, `region_n6_limit_*` (+ `.jsonl.gz`, `_HITS`, `_tally`) | eps -> 0 limit classification of every graph |
| `region_n6_eps*_only*` | exact beta-slices at eps = 1/2, 1/4, 1/10, 1/100, 1/1000 on the survivors |
| `region2d_n5_*`, `region2d_n6_shard{0,1,2}of3*` | complete two-parameter (beta, t) certification of every positivity survivor: n = 5: K5 only; n = 6: K6 + exactly three non-clique graphs |
| `COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json`, `_idx1424629.json` | the two new n = 6 counterexamples (signed off), with limit regions and pointwise witnesses |
| `n6_summary.json`, `n6_survivors.idx`, `n6_survivor_equivalence_classes.json`, `n6_never_spotcheck.json` | aggregation, survivor list, D-classes, random exact spot-check |
| `qstar/qstar_n6.json` | q_6* exact (sextic), the threshold theorem at n = 6 |

## n = 7 (complete; all 124 hits VERIFIED 2026-10-15: human hand-verification + four independent machine verifications)
| file | certifies | regenerate |
|---|---|---|
| `../n7/n7_survivors_all.txt.gz`, `../n7/n7_filter_stats.txt` | limit prefilter over all 882,033,440 digraphs: 638,042 spectral survivors (1,889 strict) | `geng 7 \| directg \| filter7` |
| `COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl`, `_summary.json`, `n7_limit_results.jsonl.gz`, `n7_limit_tally.json` | limit classification of the strict survivors: 123 non-clique limit-stable graphs + K7, with limit regions and pointwise witnesses | `n7/run_n7_limit.py` |
| `../n7/n7_pos_all_PU.txt.gz`, `../n7/n7_pos_filter_stats.txt`, `../n7/n7_pos_part?_resolved.stats.json` | exact positivity prefilter (all parameters): 9,632,133 kept, 8,909,123 survivors after exact resolution | `filter7p`, `n7/resolve_U.py` |
| `n7_curve_TOTAL.json`, `n7_curve_hits.jsonl`, `n7_curve_stats_shards.json` | Lemma R curve stage: 165 graphs stable somewhere in the legal plane, 0 escalations | `n7/run_n7_curve.py` |
| `n7_region2d/` | two-parameter certification of the 165: K7 + 124 non-clique stable positive, 40 never, 0 escalations | `n7/run_n7_region2d.py` |
| `n7_finite_summary.json`, `COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl` | the complete finite-eps map: accounting chain over all 882,033,440 digraphs; the 124 with full region certificates, limit status, D-classes | `n7/aggregate_n7.py` |
| `COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json` | the finite-eps-only graph `&FNFL~[udo?` (alpha = -1/2 exactly), region + witnesses | `n7/write_finite_certificates.py` |
| `n7_limit_hits_finite_eps_bulge_witnesses.json` | two limit candidates stable beyond beta_c at finite t | same |
| `nondegeneracy_witnesses.json` | nondegeneracy over all supports at a witness point of every n = 6 and n = 7 hit | `n7/nondegeneracy_check.py` |
| `qstar/qstar_fnfl.json`, `qstar/qstar_n7limit.json`, `qstar/qstar_n7_finite.json` | exact suprema of q (FNFL: degree-25 number; limit candidates: plastic number; finite map: see file) | `n7/qstar.py` |
| `n7_candidates_legal_witnesses.json` | one LEGAL exact witness per n = 7 counterexample (first legal stored witness + an interior cell point), as re-decided by four independent verifiers | `notes/verify_124/` |
| `n7_bulge_Rstar_exact.json` | exact sign of R* (Lemma AS) for all 124; bounded-in-t iff R* > 0 and positive at beta_c | `notes/bulge_exact/` |
| `../notes/verify_124_workflow.json`, `../notes/verify_124/` | the independent machine verification of all 124 (four implementations, 742 points, 0 undecided); statuses in the summary files read VERIFIED | |
| `theorem_U_explicit_instances.json`, `../notes/theorem_U/` | Theorem U (q_n* -> infinity): explicit n = 28 and n = 41 instances; arms-forest construction, thresholds; referee record `../notes/theorem_U_workflow.json` | |
| `theorem_A_fixed_parameter_instances.json`, `../notes/theorem_A_fixed/` | Theorem A (unbounded support size at the fixed point eps = 1/21, delta = 11/210): R_k for k = 0..15; referee record `../notes/theorem_A_fixed_workflow.json` | |
| `n8_instances_region2d.json` | two-parameter certificates of three 8-node family instances (candidates) | `n7/n8_instances_region2d.py` |
| `qstar/qstar_n7_finite.json` | q_7* = plastic number, exact, sharpened over all 298 stable cells | `n7/qstar.py --sharpen` |
| `SHA256SUMS_jsonl.txt` | checksums of the bulk files | |

## Theory (notes/)
`notes/lemmas.md` (A-E), `lemma_F_reduction.md`, `lemma_I_spectral.md`, `lemma_J_equivalence.md`,
`lemma_Bprime_2d.md`, `lemma_M_refuted.md` (Lemma M false; identity (E); Lemma R; Lemma AS; e_2 bound),
`research_wf3/theory/notes_theory.md` (Lemmas W/X, Theorems F, E, Q, Q2, Y; catalogue; refuted
conjectures), with scripts and family certificates alongside; referee and audit records in
`audit_n6_workflow_*.json`, `audit_n7_finite_workflow.json`, `research_n7_parallel_workflow.json`,
`referee_n7_theory_workflow.json`.
