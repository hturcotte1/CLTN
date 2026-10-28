# Referee fixes applied to the n = 8 finish documents

Date: 2026-10-28. Input: `notes/referee_finish_numbers.json` (191 entries: 183 OK, 7 WRONG, 1 UNSOURCED)
and `notes/referee_finish_overclaims.json` (73 findings, numbered 0-72 in file order). Files edited:
`notes/REPORT_n8.md`, `notes/FINAL_SUMMARY.md`, `notes/paper/paper.tex` (and its rebuilt `paper.pdf`),
`notes/n8_structure.md`, `notes/outreach/DRAFT_geneson.md`, `notes/outreach/DRAFT_shiu.md`,
`notes/outreach/DRAFT_curto.md`. Nothing else was touched (no certificates, no ANOMALIES.md, no
`notes/verify_n8*`, no git). Every replacement was applied by exact single-occurrence string substitution;
every new number comes from the file named next to it. Decimals remain display approximations.

Paper rebuild after the edits (`nice -n 10`, pdflatex, bibtex, pdflatex x2, in `notes/paper/`): exit 0 on
every pass, **0 LaTeX errors** (`^!` lines in `paper.log`), **0 undefined references or citations**, 0 overfull
boxes (one introduced by the new table cell was removed by moving the explanation into the caption),
50 underfull boxes (unchanged from the build checker's count), 30 pages, 621,593 bytes.

Forbidden-name scan of the seven edited files after the edits (case-insensitive, model and product names;
digraph6 strings containing `gpT` excluded): no hit. The three references to the project instructions file by
its upper-case file name were replaced by "the project instructions file".

## 1. Numbers referee: WRONG and UNSOURCED entries (8)

| # | file, place | verdict | what was changed |
|---|---|---|---|
| N1 | REPORT_n8.md 62-63, filter8 "byte-identical on all n <= 7 data" | WRONG | now "byte-identical to that filter on all 1,540,944 n = 6 digraphs and on the 638,042 n = 7 spectral survivors; ANOMALIES.md 2026-10-15" |
| N2 | REPORT_n8.md 65-67, non-printed lines "alpha > -1/2" | UNSOURCED | now "alpha >= -1/2 is what is proved ...; in the 300 sampled non-printed lines alpha > -1/2 strictly, `notes/verify_n8_2/never_side/check_ab_summary.json`" |
| N3 | REPORT_n8.md 73-74, "230 of 300 ... (ANOMALIES.md 2026-10-28)" | WRONG citation | 230/300 now cited to `notes/verify_n8_2/never_side/check_ab_summary.json` (`sample_B`) and `notes/verify_n8_2/SYNTHESIS.json`; the ANOMALIES entry's own figure (23 of 122 cross-checked B lines) is quoted separately with its source |
| N4 | REPORT_n8.md 123-124, "K8 - 2K2: 19,701" | WRONG | now "K8 - P3 (graph6 `GV~~~{`): 19,701; K8 - 2K2 (`G]~~~{`): 16,123", citing the referee's `labelg` recount and noting that the label in ANOMALIES.md 2026-10-27/28 is a mislabel (ANOMALIES.md 2026-10-22 identifies `GV~~~{` as K8 minus P3). **For the orchestrator to log in ANOMALIES.md.** |
| N5 | REPORT_n8.md 144-145, "seven container restarts" | WRONG | now "eleven container restarts (four on 15-16 October during the geng-shard phase, seven from 19 to 27 October during the per-job phase; ANOMALIES.md 2026-10-15, 2026-10-16, 2026-10-22, 2026-10-23, 2026-10-27)"; the four early restarts are the ones ANOMALIES.md 2026-10-16 lists as 19:55, 21:58, 00:1x and 01:42 |
| N6 | FINAL_SUMMARY.md 93-95, "seven container restarts" | WRONG | now "eleven container restarts: four before the per-job driver, seven during it" |
| N7 | paper.tex 846-847, "byte-identical ... on all $n\le7$ data" | WRONG | now "on all $1{,}540{,}944$ six-node digraphs and on the $638{,}042$ seven-node spectral survivors" |
| N8 | paper.tex 1409-1410 (trust base), same phrase | WRONG | now "on all six-node digraphs and on the $638{,}042$ seven-node spectral survivors" |

Numbers referee OK-with-comment entries, also acted on:

| | file, place | change |
|---|---|---|
| N9 | REPORT_n8.md Route A (a), "roughly 530-580 CPU-hours" | now "525-595 CPU-hours (131-148 wall-hours times four cores)", the referee's recomputation 524-594 rounded outward |
| N10 | REPORT_n8.md Route A (b), 15.6 ms/line | both denominators now stated: 15.6 ms per evaluated line (3,773,348) and REPORT_n7's 6.6 ms over all 8,909,123 positivity survivors, with the reason the evaluated-line figure is used |
| N11 | FINAL_SUMMARY.md "about 2,400 CPU-hours" | now "about 2,300-2,400 CPU-hours (regeneration plus curve test)", the sum of the two REPORT_n8 addends |
| N12 | REPORT_n8.md n = 9 count A000273(9) | not changed: the report already names `notes/n9_feasibility_count.json` as the only source and its validation at n <= 8 |
| N13 | paper.tex title | not changed (editorial, not a number; the referee did not ask for a change) |
| N14 | paper.tex 919-929, whether E_1 is one of the 15 | settled by the overclaim referee's exact check (E_1 isomorphic to `&GNQf^zr\IsNw`); wording changed under O16 below |

## 2. Overclaim referee: all 73 findings

Index = position in `notes/referee_finish_overclaims.json`. "Applied" means the suggested wording or a close
variant; "tightened" means a stricter wording of my own with the same content.

### Theorem Q comparison (HIGH; the record's -0.3689 lies to the LEFT of Q's limit -0.3516, so Q's limit q is larger)

| idx | file | change |
|---|---|---|
| 0 | REPORT_n8.md row 4 | applied: "lies to the left of (closer to -1/2 than) the k -> infinity limit of Theorem Q (about -0.3516, i.e. q -> 1.844); the sink-class families overtake the n = 8 record only for large k (q close to 9/5 needs k >= 286, n >= 293), and at n = 8 itself the record is attained with no sink class at all", citing `notes/referee_finish_overclaims_checks.json` block C |
| 1 | REPORT_n8.md section 5, "better core" | applied: both limits lie to the right of -0.3689, the record does not beat them; Q and Q2 exceed it only in the limit of large class size (Q needs k >= 286) |
| 2 | FINAL_SUMMARY.md section 2 | applied |
| 3 | paper.tex section n8, "lies to the right of" | applied: "lies to the \emph{left} of ... tends to $q\approx1.844$; family (Q) overtakes the eight-node record only for very large class size ($k\ge286$)" |
| 4 | paper.tex threshold picture, "further from -1/2 than any sink-class limit" | applied: "sink-free core, without any sink class; the sink-class families ... have larger limits (-0.3516 for Q, -1/4 for Q2) but reach them only for large class sizes" |
| 5 | n8_structure.md 1.8 | applied ("further LEFT (closer to -1/2) ..."), and the following sentence changed from "better cores" to "cores alone ..., although the cores themselves do not beat the asymptotic sink-class limits" |
| 6 | n8_structure.md section 2 | applied |
| 7 | n8_structure.md section 4 item 5 | applied ("does NOT beat ... -0.3689 lies to its left") |
| - | n8_structure.md section 0 item 5 (not in the findings, same topic) | tightened: added that -0.3516 "lies further from -1/2 (larger q)" so the summary cannot be read as a win for the record |

### Thresholds stated with decimals or without a recorded argument

| idx | file | change |
|---|---|---|
| 8 | FINAL_SUMMARY.md "holds ... with q >= 1.3247" (false as written: rho > 1.3247) | applied: q >= rho (the plastic number, real root of q^3 - q - 1, about 1.3247), fails for q in (1, rho) |
| 9 | FINAL_SUMMARY.md "closure contains (1, 1.7107)" | applied: fails at legal points with q arbitrarily close to every value in (1, q_8*(limit)); for q < rho by the n <= 7 map, above that by the interval of `&GMx|\Aw^qz]w`, whose lower q end is below 5/4 |
| 10 | REPORT_n8.md section 4, closure claim | applied with the full argument: (1, rho) from the n <= 7 map and the monotonicity of the failure set (Theorem F); [rho, q_8*(limit)) from the champion's lower q end in [1, 5/4] (`n8/n8_limit_HITS.jsonl`) with 5/4 < 49/37 < rho (checks block D); and the statement that no file records whether the union of all 9,133 intervals covers (1, q_8*(limit)) |

### Finite-eps-only counterexamples at n = 8 are known to exist (15 distinct, E_1 among them)

| idx | file | change |
|---|---|---|
| 11 | REPORT_n8.md comparison table | applied: "at least 15 known (the Lemma B graphs, E_1 among them; certified witnesses, not hand-verified); complete list not determined (B lines not stored)" |
| 12 | REPORT_n8.md Route B "decides whether ... exist" | applied: "lists the finite-eps-only counterexamples ... among the B and `never_full` lines (15 such graphs are already known, E_1 among them; which others exist is open)"; the next sentence (Route B is a complete map only with the alpha > -1/2 proof) was already there and stands |
| 13 | REPORT_n8.md "E_1 of Theorem E and the 15" | applied: "The 15 ... (E_1 of Theorem E is one of them, isomorphic to `&GNQf^zr\IsNw`; ...)" |
| 14 | FINAL_SUMMARY.md "settles whether finite-eps-only counterexamples exist" (HIGH) | applied: lists them among those lines (15 already known, E_1 among them); a complete finite-eps map additionally needs the alpha > -1/2 proof |
| 15 | paper.tex after Theorem n8, "may sit among the unstored B lines" | applied: they exist (15 distinct graphs of Lemma Bmech, E_1 among them, alpha = -1/2, among the unstored B lines; at n = 7 exactly one), further ones may sit among the B lines, the 7,078 never_full lines or, unless Question char is answered, the alpha > -1/2 lines |
| 16 | paper.tex "$E_1$ ... and the $15$ graphs" (lifts paragraph) | applied: "The $15$ graphs satisfying Lemma Bmech ($E_1$ of Theorem famE is one of them)" |
| 17 | paper.tex open question, same double count | applied: "$E_1$ of Theorem famE is one of the $15$ graphs satisfying Lemma Bmech" |
| 18 | paper.tex threshold table, "$\ge9{,}133$" | applied as "$\ge9{,}148$, not mapped" in the cell (the long form overflowed the column) with the explanation "counts the 9,133 limit counterexamples plus the 15 known finite-eps-only graphs of Lemma Bmech" moved into the caption |
| 19 | n8_structure.md section 2, "E_1 and the 15 'cube' candidates" | applied: the 15 finite-eps-only graphs (E_1 one of them; exact legal witnesses); they are n = 8 analogues of `&FNFL~[udo?` and live among the unstored B lines |
| - | REPORT_n8.md status box and FINAL_SUMMARY.md n = 8 row (not in the findings, same topic) | tightened: both scope statements now say that 15 such eight-node graphs are already known (Lemma B, E_1 among them, exact legal witnesses, no hand check) and that the complete finite-eps list is not determined |

### The paper's open question on the finite map (HIGH/MEDIUM)

| idx | file | change |
|---|---|---|
| 20 | paper.tex "would settle whether $q_8^*=q_8^*(\infty)$" | applied: the re-run plus curve test on the B lines and on the 7,078 never_full lines plus the bulge rule settles it for all graphs with alpha <= -1/2; the alpha > -1/2 lines need a proof answering Question char or the full finite pipeline |
| 21 | paper.tex "the bulge criterion decided the finite map from the limit data" | applied (variant): "At $n=7$ the finite map was computed by the full chain of Table chain, and the bulge criterion (Theorem AS) accounts for its difference from the limit map" |
| 22 | paper.tex "out of reach by enumeration alone" | applied: the full finite map needs the n = 7 pipeline over the 801,558,826,511 restricted digraphs, of the order of $4\times10^4$ CPU-hours extrapolated from measured rates |
| - | REPORT_n8.md Route B bullet "Either route also yields the finite-eps q_8*" (not in the findings, same defect) | tightened: Route A yields it; Route B only under the same proviso on the alpha > -1/2 lines |

### The 797 blow-ups not in the list (limit-only statement; sampled points)

| idx | file | change |
|---|---|---|
| 23 | FINAL_SUMMARY.md "are not counterexamples (797 of 858)" | applied: "not limit counterexamples ...; their finite-eps status is open for those neither dominated nor never-positive" |
| 24 | paper.tex "exactly not limit-stable at legal points inside their parent's interval" | applied: "not limit counterexamples (prefilter: 690 dropped, 47 B, 60 never in the limit stage), and none is stable and positive at 4,782 sampled legal points" (`notes/verify_n8_2/completeness_twins/part2_summary.json`: DROPPED 690, B 47, S 60 = 43 never_positivity + 17 never_full) |
| 25 | REPORT_n8.md section 3 | applied with the same prefilter breakdown and the note that the 4,782 points are a sample, not an exact decision of the whole interval; also REPORT_n8.md row 6 ("the other 797 are not limit counterexamples ...") and section 6 ("797 not limit counterexamples; 4,782 sampled legal points, 0 stable and positive"), and paper.tex section 8.3 ("not limit counterexamples and not stable and positive at 4,782 sampled legal points") |
| 26 | n8_structure.md "most blow-ups die" | applied: "most blow-ups are not limit counterexamples" |

### Sink deletions compared only with the 123 limit counterexamples

| idx | file | change |
|---|---|---|
| 27 | REPORT_n8.md section 5 | applied: 1,885 / 1,972 refer to n = 7 limit counterexamples; in 2 of the 1,972 (`&GNbR^x~[zR[?`, `&GNbR^z~[rR[_`) a sink deletion leaves the finite-eps-only `&FNFL~[udo?`; 1,970 have no sink deletion that is any n = 7 counterexample (checks block A) |
| 28 | paper.tex "leaves a graph that is not a counterexample" | applied: "not a seven-node limit counterexample ($1{,}970$ leave no seven-node counterexample at all; in two, a sink deletion gives &FNFL~[udo?)" |
| 29 | n8_structure.md section 0 item 4 | applied; and section 1.4 ("Deleting one sink node") tightened the same way, noting that `families.py` compares with the limit list only |

### "Consequences of the necessary conditions" (unproved causal claim)

| idx | file | change |
|---|---|---|
| 30 | FINAL_SUMMARY.md | applied: "do not distinguish counterexamples (whether they follow from the necessary conditions is open: S1+ is checked for n <= 9 only, and U1/U2 are not proved)" |
| 31 | REPORT_n8.md section 8 item 3 | applied: "shared by almost all strict survivors and so are not a signature; that U1-U3 follow from the necessary conditions is not proved" |
| 32 | n8_structure.md section 0 item 3 | applied: "shared by the graphs that pass the necessary conditions (whether they are implied by them is not proved for U1-U3)" |

### Conjecture S1 at n <= 9 is verified in the limit sense only

| idx | file | change |
|---|---|---|
| 33 | FINAL_SUMMARY.md section 2 | applied: "S1 (for limit counterexamples) is verified exhaustively ..." |
| 34 | FINAL_SUMMARY.md section 4 | applied: "S1 (limit sense) is already verified there" |
| 35 | paper.tex introduction | applied: "for limit counterexamples the last property is verified exhaustively up to $n=9$" |
| 36 | paper.tex after Conjecture S | applied: "(for limit counterexamples, as the conjecture is stated)" |
| 37 | REPORT_n8.md row 5 | applied: "verified, in the limit sense, for every n <= 9" |

### Heredity extrapolated beyond n = 8

| idx | file | change |
|---|---|---|
| 38 | paper.tex | applied: "At $n=7$ and $n=8$ most counterexamples are not one-node extensions of smaller ones" |
| 39 | FINAL_SUMMARY.md | applied: "the n = 7 and n = 8 data say most counterexamples at those sizes are new rather than inherited" |
| 40 | REPORT_n8.md section 8 | applied: "we expect, but do not know, that most n = 9 counterexamples are not extensions" |
| 41 | n8_structure.md section 3 | applied: "we expect them to be the majority; this is not known" |

### Catalogue membership is a recall check, not a completeness check

| idx | file | change |
|---|---|---|
| 42 | paper.tex | applied: "an independent recall check on the proved families (it cannot detect graphs outside them)" |
| 43 | n8_structure.md 1.6 | applied |

### Human verification caveats

| idx | file | change |
|---|---|---|
| 44 | FINAL_SUMMARY.md section 3 item 4 | applied: the user recorded the n = 7 hand check as not fully confident and conditioned sign-off on the machine re-verification (ANOMALIES.md 2026-10-15); also the n = 7 row of the table in section 1 |
| 45 | REPORT_n8.md comparison table, verification row | applied: "hand (2026-10-15, recorded by the user as not fully confident)"; also section 6 ("a check the user recorded as made without full confidence") |

### "Refereed" means internal review

| idx | file | change |
|---|---|---|
| 46 | paper.tex theorem labels (famF, famE, QQ2Y, U, A) and the audits subsection | applied: all five labels now "proved, internally reviewed"; the audits subsection says the theorem reviewers were "internal to the project (not external peer review)"; the sentence after Lemma limit now says "reviewed internally ... internal reviewers". FINAL_SUMMARY.md: "proved and refereed" -> "proved, internally reviewed (not peer-reviewed)", and the Theorem A / Theorem U / audit lines likewise |
| 47 | paper.tex second campaign synthesis | applied: "(one referee's MAJOR, an unfinished re-derivation of five never_full certificates, was covered by two other checks and reclassified)" |

### Literature claims

| idx | file | change |
|---|---|---|
| 48 | REPORT_n8.md row 8 | applied (limited arXiv search, index to 27 Oct 2026, no other preprint servers; status "searched, limits in section 7") |
| 49 | FINAL_SUMMARY.md section 3 item 6 | applied ("no overlapping work found in a limited arXiv search, index to 27 Oct 2026; nothing on n >= 7 or n = 8 was found") |

### Other wording

| idx | file | change |
|---|---|---|
| 50 | FINAL_SUMMARY.md preamble "every number ... certificate" | applied: "either backed by a re-runnable exact file (named in the report that owns the number) or explicitly an estimate or a display approximation" |
| 51 | FINAL_SUMMARY.md "grows with n without any known bound" | applied: unbounded in n (Theorem U); the only upper bound is Geneson's n - 2; the growth rate is open |
| 52 | FINAL_SUMMARY.md n = 9 cost | applied: "at least about 4 million CPU-hours at the n = 8 per-digraph rate (an optimistic extrapolation)" |
| 53 | REPORT_n8.md n = 9 cost | applied: "measured n = 8 rate ... (optimistic for n = 9, whose matrices are larger) ... at least about 4.1 million CPU-hours" |
| 54 | FINAL_SUMMARY.md n = 5 reproduction | applied: "at eps = 1/4, delta = 1/2 with zero mismatches" (`certificates/stage1_n5_diff.json`) |
| 55 | REPORT_n8.md chain table, "Lemma A.2" for the edge bound | applied: "`edge_bound_ctln` in `notes/lean/CtlnCert/CtlnCert/EdgeBound.lean`, Lean-checked; hypotheses: legal point, m >= 2, full support Hurwitz stable" (the correct Lemma A.2 references in sections 1 and 9 were left) |
| 56 | REPORT_n8.md row 1, "four pointwise witnesses" | applied: "four stored pointwise checks (all four stable and positive for 9,051 records; section 3, 'Witness hygiene'), and a certified legal interior witness" |
| 57 | paper.tex Theorem q8 proof sketch, Lemma open cited for non-attainment | applied: "at $\beta_c$ itself $\alpha=-\gamma$, so the graph is not limit-stable there (Lemma limit)" |
| 58 | paper.tex "about 1.8 million digraphs per second" | applied: "up to $1.8$ million digraphs per second per process ($1.5$--$1.7$ million per second aggregate in production)" |

### Outreach drafts

| idx | file | change |
|---|---|---|
| 59 | DRAFT_geneson.md subject | applied: "exact for n <= 7, eps -> 0 limit at n = 8" |
| 60 | DRAFT_geneson.md "mapped it exactly" | applied: "exactly for n <= 7 and in the eps -> 0 limit at n = 8" |
| 61 | DRAFT_geneson.md hand check | applied: "checked by hand (n = 7 not with full confidence) and confirmed by four independent exact re-implementations" |
| 62 | DRAFT_geneson.md families bullet | applied: "(our proofs, internally reviewed, not yet peer-reviewed)" |
| 63 | DRAFT_shiu.md subject | applied |
| 64 | DRAFT_shiu.md "always is" | applied: "stable at some legal parameter in all 192 such cases (certified, not hand-verified)" (`certificates/n7_proper_supports_summary.json`, `node6_no_in_verdicts.stable_somewhere` = 192) |
| 65 | DRAFT_shiu.md hand check | applied |
| 66 | DRAFT_shiu.md families bullet | applied |
| 67 | DRAFT_curto.md subject | applied |
| 68 | DRAFT_curto.md hand check | applied |
| 69 | DRAFT_curto.md families bullet | applied |

### Product-name rule

| idx | file | change |
|---|---|---|
| 70 | REPORT_n8.md preamble | applied: "Conventions are those of the project instructions file, section 2" |
| 71 | REPORT_n8.md section 6 | applied: "four verifiers written from the project instructions file's definitions" |
| 72 | n8_structure.md preamble | applied: "Conventions are those of the project instructions file, section 2" |

## 3. Not changed, and why

* `notes/verify_n8_2/SYNTHESIS.json`, field `synthesizer`, contains a model name (reported by the build
  checker and both referees). The file is under `notes/verify_n8*`, which this task may not modify. **For the
  orchestrator.**
* The `notes/outreach/README.md` was not in the referee findings and needed no change (it already states the
  limit-only scope and the no-hand-check status).
* ANOMALIES.md entries that the referees found wrong or incomplete were not edited (append-only, off limits).
  Items for the orchestrator to log there: (a) the "K8 - 2K2 (GV~~~{): 19,701" label in the 2026-10-27/28
  entries is K8 - P3; K8 - 2K2 (G]~~~{) has 16,123 strict lines; (b) E_1 of Theorem E (`&GNbR^z~[zP[_`) is
  isomorphic to the Lemma B graph `&GNQf^zr\IsNw`, so 15 distinct finite-eps-only eight-node graphs are known;
  (c) in five listed n = 8 graphs a one-vertex deletion is `&FNFL~[udo?`, always at a sink node, two of them
  among the 1,972; (d) the Theorem Q comparison in the finish documents was backwards (now corrected);
  (e) the run saw eleven container restarts, not seven; (f) the filter8 "byte-identical" validation covered
  all n = 6 digraphs and the 638,042 n = 7 spectral survivors, not all n = 7 data.
* Pre-existing n <= 7 passages of the paper that neither referee flagged were left as they were, except the
  "refereed" wording (finding 46), which was made consistent throughout.
