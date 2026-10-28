# Figure referee findings: what was applied

Source: `notes/paper/figures/referee_figures.json` (verdict HOLD: 0 FATAL, 0 MAJOR, 3 MINOR, 7 NOTE).
Files touched: `notes/paper/figures/make_figures.py` (re-run; all five PDFs and `figures_data.json`
regenerated, every assertion passing), `notes/paper/paper.tex` (two captions). Nothing else was modified;
git was not run. The paper rebuilds (pdflatex, bibtex, pdflatex x2) with 0 errors, 0 undefined references,
0 overfull boxes; build artefacts removed. Final page count: **33** (unchanged by these fixes). Figures sit
on pp. 11 (Fig. 1), 17 (Fig. 2), 18 (Fig. 3), 22 (Fig. 4), 24 (Fig. 5).

## MINOR (all applied)

1. **Caption of Fig. 3 (`fig:hist8`): "the bins refine the histogram of Section 6.4".**
   The first figure bin (1, 21/20] merges three Section 6.4 bins (37 + 541 + 2,777 = 3,355), so "refine" was
   wrong below 21/20. Changed to "the counts are consistent with the histogram of Section 6.4 (coarser below
   21/20, finer above)". No counts changed.

2. **Fig. 5 (`fig:thresholds`): unexplained white streaks in the n = 7 and n = 8 stacks.**
   - `make_figures.py` now counts the rows whose lower end is of kind `root` (interval starts above q = 1),
     asserts the counts (n = 6: 0, n = 7: 6, n = 8: 971) and records them in `figures_data.json`
     (`fig2_thresholds.pdf.white_gaps`).
   - The n = 8 envelope is now drawn with `fill_betweenx(..., step="mid")`: one horizontal band per row
     between its own drawn lower and upper end, no interpolation between neighbouring rows, so a gap is a
     rectangle at that row's height rather than a wedge across its neighbours. (The referee's alternative,
     a LineCollection of 9,133 sub-point segments, was not used: at a 0.011 pt row pitch it renders as the
     same envelope but with 9,133 anti-aliased hairlines, which is less faithful in print.)
   - Caption now says: rows are sorted by upper end "(ties by lower end)"; the envelope is "one horizontal
     band per row between its lower and upper end"; and "White gaps inside a stack are rows whose interval
     starts above q = 1 (lower end a root rather than beta = 2): 6 of the 123 rows at n = 7 and 971 of the
     9,133 at n = 8, grouped because rows with equal upper ends are sub-sorted by lower end."

3. **In-figure text sizes (base 8 pt, ticks 7.5 pt, bar counts 6.3 pt, funnel annotation 6.8 pt).**
   Figures are placed 1:1 at `\textwidth`, so these printed as-is beside a 10.95 pt body. Sizes are now
   named constants: base/labels 9 pt (`FS_BASE`), panel titles 9.5 pt (`FS_TITLE`), ticks and legend
   8 pt (`FS_SMALL`), and nothing below 7.2 pt (`FS_MIN`: the bar counts of Fig. 3 and the two funnel
   annotations). The printed counts on the left histogram were kept rather than dropped (the referee's
   alternative), since the right-hand bars with counts 1 to 16 are sub-pixel and rest on the labels alone;
   at 7.2 pt all 32 rotated labels are legible at 150 dpi. Figure heights were raised slightly (Fig. 1
   3.55 -> 3.65 in, Fig. 3 funnel 2.35 -> 2.45 in, Fig. 4 3.75 -> 3.85 in, histograms 2.5 -> 2.7 in) and
   the funnel's label margin widened (0.30 -> 0.33) so the larger text fits; every figure was rasterised and
   checked for clipping. One consequence fixed on the way: at 9.5 pt the right-hand title of Fig. 4 ran
   off the panel, so it was shortened from "29-arc graph (edge bound attained), upper end ~1.0131, 29 arcs"
   to "edge bound attained, upper end ~1.0131, 29 arcs" (the arc count was printed twice).

## NOTE (applied where cheap; the rest left, with reasons)

4. **File names vs figure numbers** (fig3_funnel.pdf is Figure 2, fig5_histograms.pdf is Figure 3,
   fig2_thresholds.pdf is Figure 5). Left as is: the names are cited in `figures_data.json`,
   `referee_figures.json` and the captions, and figure numbers depend on float placement; renaming would
   churn every reference for no gain in correctness.

5. **Fig. 5 right margin (q_8*(infinity) label 1.3 mm from the edge, no `bbox_inches='tight'`).** Applied:
   the n = 8 label is now placed to the left of its tick (`ha="right"`), so its length can never reach the
   margin. `bbox_inches='tight'` was deliberately not added, so the PDFs keep the exact `\textwidth` width.

6. **q* ticks at the midpoint of the stored isolating interval (q_6* interval of width 3.5e-4, tick 1.4e-4
   right of the record rows' drawn ends).** Applied: each q_n* isolating interval is now narrowed by the
   same exact bisection used for the rows (`narrow_root`, width < 1e-9, sign change asserted, the narrowed
   interval asserted to still contain the largest drawn upper end), the tick is placed at the midpoint of
   the narrowed interval, both intervals are recorded in `figures_data.json` (`q_star_ticks.*.narrowed_interval`),
   and the caption says so. No mathematical content changes: the narrowed interval is a sub-interval of the
   certified one bracketing the same simple root.

7. **`pitch_pt` used 5.9 data units for an axis spanning 6.3.** Applied: the pitch is now derived from the
   axes box, figure height and `YLIM` (named constants), i.e. 6.3 units. The n = 7 line width stays at the
   0.5 pt clamp (true pitch 0.52 pt), so the drawing is unchanged.

8. **Funnel bar length is log10(count).** Left: the caption states it, every count is printed, and the
   referee notes the ordinal ramp on an ordered funnel is permitted and reads in grayscale.

9. **Underfull hbox at the Fig. 1 caption (badness 1057, unbreakable `\texttt` path).** Left: harmless,
   and breaking the file path would make it harder to copy.

10. **Fonts: STIXGeneral/DejaVu Sans Mono in the figures against a Computer Modern body.** Left:
    `text.usetex` would require a full TeX install on the figure-building machine; both fonts are embedded
    as Type 42 subsets and the figures are self-contained.

11. **Fig. 3 `[t]` float above the "6.3 The 9,133 graphs" heading.** Left: `[t]` is the standard choice
    and the figure is referenced on the same page; forcing `[h]` would risk a worse break.
