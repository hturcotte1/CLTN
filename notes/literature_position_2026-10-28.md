# Literature position, re-check of 2026-10-28

Retrieval date 2026-10-28. Baseline: notes/REPORT_n7.md section 9 (as of 2026-10-14), paper.tex "Related work",
notes/literature_twins_check.md (re-check of 2026-10-13), notes/overlap_report_github_repos.md.
Nothing was contacted, submitted or posted. No code from third parties was run.

## 1. Bottom line

Nothing new that bears on this project has appeared since 2026-10-14 (and nothing since the 2026-10-13 re-check).
No n = 8 result, no general-n classification, no improvement of Geneson's threshold, no follow-up by Geneson, Curto,
Morrison, Shiu, Parmelee, Santander or Londono Alvarez on target-free cliques or stable supports. Scoop check: negative.
Limits of this check are in section 5.

## 2. What was checked and what came back

Method: WebFetch of arXiv's own search pages (sorted newest first; the search index showed entries up to 27 Oct 2026),
arXiv abstract pages, the Zenodo record, `git ls-remote` on the two GitHub repos, and WebSearch.

| Query / page | Result |
|---|---|
| arXiv search "threshold-linear networks", all fields | 25 hits; newest relevant is Geneson 2607.21396 (23 Jul 2026), then Londono Alvarez 2606.07336, Zang-Curto 2606.00373. Nothing after 23 Jul. |
| arXiv search "combinatorial threshold-linear" | 11 hits; newest is again 2607.21396. |
| arXiv search "CTLN" | 12 hits; same list, nothing new. |
| arXiv search "target-free" | 76 hits; only 2607.21396 is on topic (others are unrelated: ML, astronomy). |
| arXiv author search Geneson (69 hits) | Since 2607.21396: 2608.12604, 2608.17294, 2609.25107, 2610.02707, all extremal combinatorics / sequences, none on CTLNs. No v2 of 2607.21396: abs page lists v1 only (23 Jul 2026, 12 KB). |
| arXiv abs 2511.05517 (Geneson, non-minimal fixed points) | v1 only (26 Oct 2025). |
| arXiv author search Curto (411 hits, many other Curtos) | Only CTLN item in 2026 is Zang-Curto 2606.00373 (29 May 2026). |
| arXiv author search Morrison, Shiu | No CTLN or TLN items in the first 25 hits (all other fields). |
| arXiv author page /a/curto_c_1 (a possibly cached page) | newest item 2606.00373. |
| Parmelee / Santander / Londono Alvarez author search | Query returned HTTP 400 (not retrieved). Covered indirectly by the topic searches above, which list Londono Alvarez 2606.07336 (v2 revised 28 Sep 2026) and Parmelee 2109.03198, 2107.10244 as the newest items of those authors. |
| export.arxiv.org API (several queries) | HTTP 429 (rate limited), not retrieved. The arxiv.org/search pages were used instead. |
| arxiv.org/a/geneson_j_1 | HTTP 404 (not retrieved). |
| Zenodo 10.5281/zenodo.23089682 ("Six Is Minimal") | Version 1.1.0-candidate, published and last modified 2026-10-01; no later version listed; no mention of n >= 7. |
| github.com/ipitchford/ctln-six-is-minimal (`git ls-remote`) | HEAD e28a309ed6de9da7d93dfebac29c35a92a65a312, unchanged from the 2026-10-13 report (tag v1.1.0-candidate). |
| github.com/ChickenNungets/CTLN (`git ls-remote`) | HEAD 52ad0003a34699ce47fe9a26cfa4b8c6919f088b, unchanged from the report (52ad000, 2026-08-08). |
| WebSearch (n = 8 classification, follow-ups, "Six Is Minimal", Curto lab 2026) | Nothing on n >= 7, n = 8 or the gap. |

## 3. Items that are in the record but not in our baseline text

These are not scoop risks; they are optional citations.

1. Londono Alvarez, "Fixed point compositionality via low-rank gluing rules in inhibition-dominated threshold-linear
   networks", arXiv:2606.07336, v1 5 Jun 2026, v2 28 Sep 2026, q-bio.NC, 41 pages, "Submitted to SIAM Journal on Life
   Sciences" (arXiv comment). Abstract (from the abs page): global fixed points of low-rank gluings are combinations of the
   local fixed points of the modules; complete characterization for rank-1 gluings; extends CTLN fixed-point decomposition
   rules to generalized CTLNs. It is about fixed-point decomposition, not target-free cliques, stable supports or small-n
   classification. Relation to us: nearest in spirit to our closed-twin blow-up and in-twin-switch results (Lemmas J, W);
   NOT read in full, so no claim of overlap or non-overlap is made beyond the abstract. Recommended: read before
   finalizing the Related work text.
2. Zang and Curto, arXiv:2606.00373 (29 May 2026), sequential chaotic oscillations in excitatory-inhibitory TLNs. Already
   noted in the 2026-10-13 re-check; unrelated.
3. Geneson, arXiv:2511.05517 (Oct 2025): stable non-minimal fixed points of TLNs; the Geneson 2607.21396 introduction cites
   it as a non-combinatorial disproof of the minimal-support conjecture (secondhand from a Pith page, see below). Already
   in notes/geneson_2607.21396.md as "not retrieved"; the abs page confirms v1 only.
4. Pith page for 2607.21396 (https://pith.science/paper/2607.21396): lists no citing paper other than 2511.05517, shows no
   human review; its reference list names Londono Alvarez, Morrison and Curto 2026 as reference [10] (title not given in
   the retrieved summary). Secondhand; the Pith page is an automated aggregator.
5. A search-engine snippet (not the paper) described 2511.05517 as stating the conjecture was proven for n <= 4. Unverified
   and probably a paraphrase error; our baseline already treats supports of size <= 4 as known (CGM). Do not cite.

## 4. Recommended wording changes

### 4.1 REPORT_n8, section "Literature position" (proposed text, as of 2026-10-28)

> As of 2026-10-28 we know of no published or posted work that classifies n = 8, treats n = 7, exhibits finite-eps-only
> counterexamples, or improves Geneson's threshold (arXiv:2607.21396, v1 only, 23 Jul 2026). A search of arXiv (topic
> searches "threshold-linear networks", "combinatorial threshold-linear", "CTLN", "target-free", and author searches for
> Geneson, Curto, Morrison, Shiu) on 2026-10-28 found no CTLN item dated after 23 Jul 2026 apart from a revision of
> Londono Alvarez, arXiv:2606.07336 (v2, 28 Sep 2026), which concerns fixed-point composition under gluings and not
> target-free cliques. The unrefereed preprint "Six Is Minimal" (Zenodo 10.5281/zenodo.23089682, v1.1.0-candidate,
> 2026-10-01) and its repository were unchanged at the time of the check; it covers n <= 5, Geneson's graph and a family on
> 5r + 6 nodes. The comparison with it in the n = 7 report (section 9) is unchanged. Our n = 8 result is a limit map
> only; it is verified by machine (two independent campaigns) and not by hand.

Changes relative to the n = 7 section 9: update the date; add the one-sentence arXiv search statement; add the
2606.07336 sentence; state that the Pitchford and Vazquez repositories were re-checked by `git ls-remote` and are
unchanged (commit ids in section 2). Keep "claims no priority either way" for Six Is Minimal.

### 4.2 paper.tex, Related work

- Keep the existing paragraph on Geneson, Curto-Morrison, Curto and "Six Is Minimal" unchanged.
- Add one sentence after the Geneson sentence: "Geneson's preprint has no later version at the time of writing (checked
  2026-10-28)." Delete if the paper is edited later than that date without a re-check.
- Optional, only after reading it: cite Londono Alvarez (arXiv:2606.07336) for fixed-point decomposition under gluings,
  and say that our closed-twin blow-up and in-twin-switch lemmas are special cases of local-to-global compositions of
  that kind only if the paper confirms it. Until read, do not write that sentence.
- Do not say "no one has done n = 8" as an absolute; say "we know of no".
- The paper's statement that it is based on reading their repository, not running their code, stays true; this check also
  did not run anything.

## 5. What could NOT be verified

- Coverage of arXiv search is only as fresh as its index (entries up to 27 Oct 2026 were visible); papers posted on
  27-28 Oct may not yet be searchable. Repeat the check before any posting.
- Author searches for Parmelee, Santander, Londono Alvarez returned HTTP 400, and the arXiv API returned HTTP 429, so those
  authors are covered only through topic searches. The Morrison and Shiu searches match many unrelated authors; only the
  first 25 hits (sorted by announcement date) were read, which reach back to late Jul 2026 (Morrison) and Jul 2026 (Shiu), so
  a CTLN paper by these authors in that window would have appeared.
- Preprint servers other than arXiv (bioRxiv, OSF, Zenodo search, SSRN) were not searched, and the Zenodo and GitHub checks
  were limited to the two known records. The GitHub API search is blocked in this session, so a new unknown repository
  would not be found. WebSearch (US-only, search-engine index) found nothing.
- The text of 2606.07336 was not read beyond its abstract page; the text of 2607.21396 was not re-read in this task
  (the abstract page was re-fetched; its version count was confirmed).
- Whether Geneson has told anyone about plans to map the gap is unknown; nothing public shows it.
- Journal acceptance or publication status of any of these was not checked.

## 6. Sources (all retrieved 2026-10-28 unless noted)

- https://arxiv.org/abs/2607.21396 (v1, 23 Jul 2026)
- https://arxiv.org/abs/2511.05517 (v1, 26 Oct 2025)
- https://arxiv.org/abs/2606.07336 (v1 5 Jun 2026, v2 28 Sep 2026)
- arXiv search result pages for the queries in section 2 (https://arxiv.org/search/ with the query terms listed there)
- https://zenodo.org/records/23089682
- https://github.com/ipitchford/ctln-six-is-minimal (HEAD e28a309), https://github.com/ChickenNungets/CTLN (HEAD 52ad000)
- https://pith.science/paper/2607.21396
