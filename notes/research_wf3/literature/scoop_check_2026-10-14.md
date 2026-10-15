# Scoop check 2026-10-14 (arXiv/WebFetch/WebSearch; shell curl not used)

## Relevant items
- 2607.21396 Geneson, 23 Jul 2026, v1 only, math.CO/cs.DM. 6-neuron counterexample (delta=29eps/25) + proof for q >= n-2-(n-3)eps/2. No open-problems section; does not pose the gap question. Refs: CGM 2024, Curto-Morrison Notices AMS, Geneson 2511.05517, MDIC.
- 2511.05517 Geneson, Oct 2025: 3-neuron competitive TLN with stable non-minimal fixed point (general TLN, not CTLN). Background only.
- 2510.05098 Curto, Oct 2025 (domination). No revision seen.
- Geneson listing (arXiv search by author): after 2607.21396 only unrelated combinatorics (2608.12604, 2608.17294, 2609.25107, 2610.02707).
- Curto listing: 2606.00373 (Zang-Curto, E-I TLN sequential chaotic oscillations, v2 3 Oct 2026). Not about stable CTLN supports.
- 2606.07336 Londono Alvarez (v1 5 Jun, v2 28 Sep 2026): low-rank gluing rules, fixed point compositionality. No stable-support / target-free content.
- Morrison / Shiu: nothing new on CTLN since 2025 (Shiu results are chemical reaction network papers).
- 2507.10383 Cohen-Lengyel (v2 Sep 2026): dense attractor stability, not CTLN. Irrelevant.
- No arXiv paper found citing Geneson 2607.21396 (searched; not exhaustive: no citation-index tool).

## OVERLAP: non-arXiv, unrefereed
- "Six Is Minimal: Exact Certification and Infinite Families of Stable Nonclique CTLNs", author "Anonymous" (GitHub ipitchford/ctln-six-is-minimal; Zenodo 10.5281/zenodo.23089682), 1 Oct 2026, v1.1.0-candidate, CC0/MIT, explicitly unrefereed, AI-review-responsive.
  Claims: (1) min stable nonclique support size is exactly 6 (exact exhaustive check of graphs <=5 vertices: 9,846 iso classes / 1,052,741 labelled, full legal domain); (2) six-neuron witness = Geneson's graph (not claimed new); (3) infinite family on 5r+6 vertices (6, 11, 16,...) with shrinking margin, no 7-vertex member; (4) full-support phase region for seed graph via roots of a sextic and a quartic.
  Does NOT: classify all n=6 graphs, count n=6 counterexamples, treat n=7, give q-thresholds, improve Geneson's q bound.
  Relation to us: overlaps Stage 1 / sizes<=5 (we have n<=5 full-support plus known supports<=4), and the seed (Geneson) graph's full-support region (our region2d). Not overlapping: our complete n=6 map (K6 + 3 non-clique), twin partner and idx 1387528, n=7 limit candidates (123), finite-eps-only n=7 candidate, Lemma R/I/J, "no q*>1" statement.
  Caveats: their "up to five" counts (9,846 / 1,052,741) differ from classic 9,608 total digraphs on 5 nodes (A000273) and n=5 labelled 2^20=1,048,576; 9,846 and 1,052,741 look like they include something else (e.g. loops or counted differently) -- unverified, worth a diff if relevant; could not inspect code (only README/Zenodo page fetched).
  NOT verified by me: any of their claims; PDF/code not read.

## Verdict
No published (arXiv/journal) result found overlapping: complete n=6 map, n=7 candidates, finite-eps-only counterexample, or q-threshold improvement below n-2. Partial overlap with an anonymous unrefereed Oct 2026 Zenodo/GitHub preprint (n<=5 exclusion, Geneson seed phase region, 5r+6 family). Repeat check periodically; consider reading their CLAIMS.md / certificates for n=5 cross-validation and seed-region comparison.
Limits: search tools are web-based, coverage of Oct 2026 listings incomplete (search pages truncated at ~40 hits); no arXiv listing for math.CO/q-bio.NC/math.DS day-by-day was scanned.
