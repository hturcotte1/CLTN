#!/usr/bin/env python3
"""Build the five vector figures of notes/paper/paper.tex from the project's certificate files.

Re-runnable:  python3 notes/paper/figures/make_figures.py   (from anywhere; paths are resolved from this file).

Reads ONLY these data files (nothing else is consulted, nothing is recomputed from the dynamics):
  certificates/stage0_geneson_2607.21396.json          Geneson's six-node graph (A_ours)
  certificates/n6_summary.json                         the three n = 6 counterexamples and their limit intervals
  certificates/n6_survivor_equivalence_classes.json    which of the three share Geneson's class (classes_with_hits)
  certificates/qstar/qstar_n6.json                     q_6^* exact (minimal polynomial + isolating interval)
  certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl  123 n = 7 limit counterexamples with their intervals
  certificates/qstar/qstar_n7limit.json                q_7^*(limit) = plastic number, isolating interval
  certificates/n8_summary.json                         n = 8 aggregate: per-record q rows with refined isolating
                                                       intervals, edge counts, upper-end histogram, q_8^*(limit)
  n8/n8_limit_HITS.jsonl                               the 9,133 n = 8 records (d6, edges, stable_intervals)
  n8/n8_prefilter_TOTAL.txt                            prefilter totals (accounting chain)
  n8/n8_limit_tally.json                               limit-stage tally (accounting chain)
  n8/k8decomp/g6tools.py                               the project's reference digraph6 decoder (cross-check)

Writes:  fig1_counterexamples.pdf  fig2_thresholds.pdf  fig3_funnel.pdf  fig4_n8_graphs.pdf  fig5_histograms.pdf
         figures_data.json   (exactly what each figure draws: strings, vertex orders, arc lists, rational
                              endpoints, counts)

Conventions (project conventions file, section 2): A[i][j] = 1 iff arc i -> j.  digraph6: '&' + chr(63+n) + the n*n adjacency
bits row by row, bit (i,j) = arc i -> j, six bits per character, char = value + 63, most significant bit first.

Exact arithmetic (fractions.Fraction) for every decision: bin membership, sorting of interval ends, the
narrowing of isolating intervals.  Floats appear only when a rational is handed to matplotlib for pixel placement.
"""
import json
import os
import sys
from fractions import Fraction as F

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "n8", "k8decomp"))
import g6tools  # noqa: E402  (reference decoder; used only to cross-check)

P = lambda *a: os.path.join(ROOT, *a)  # noqa: E731

# ----------------------------------------------------------------------------------------------- style
INK = "#0b0b0b"        # primary text / single arcs
INK2 = "#52514e"       # secondary text
MUTED = "#898781"      # axis labels
AXIS = "#c3c2b7"       # baselines
GRID = "#e1e0d9"       # hairline grid
BLUE = "#2a78d6"       # the one series colour (categorical slot 1, validated)
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]  # ordinal blue ramp (validated --ordinal)
MISSING = "#c3c2b7"    # light dashed non-adjacent pairs
SURFACE = "#ffffff"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.linewidth": 0.5,
    "axes.edgecolor": AXIS,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "text.color": INK,
    "axes.labelcolor": INK2,
    "savefig.dpi": 300,
})
TEXTWIDTH_IN = 6.7  # 8.5in - 2 * 0.9in (paper geometry); figures are placed at \textwidth, i.e. 1:1
# In-figure text sizes (pt).  The figures are not scaled when placed, so these print as given next to an
# 10.95 pt body: nothing below 7 pt (referee), titles and labels at 9 pt.
FS_BASE = 9.0
FS_TITLE = 9.5
FS_SMALL = 8.0
FS_MIN = 7.2

DATA = {"generator": "notes/paper/figures/make_figures.py", "figures": {}}


# ----------------------------------------------------------------------------------------------- digraph6
def decode_d6_mine(s):
    """Own decoder, written from the convention in the docstring; cross-checked against g6tools.decode_d6."""
    assert s[0] == "&", "digraph6 strings start with '&'"
    n = ord(s[1]) - 63
    assert 1 <= n < 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, "bad digraph6 character"
        bits.extend(((v >> (5 - b)) & 1) for b in range(6))
    assert len(bits) >= n * n
    assert all(b == 0 for b in bits[n * n:]), "padding bits must be zero"
    A = [bits[i * n:(i + 1) * n] for i in range(n)]
    assert all(A[i][i] == 0 for i in range(n)), "no loops"
    return n, A


def encode_d6_mine(A):
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    bits += [0] * ((6 - len(bits) % 6) % 6)
    out = "&" + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out += chr(63 + v)
    return out


def decode_checked(s):
    """Decode with both decoders and assert agreement (and round trip through both encoders)."""
    n1, A1 = decode_d6_mine(s)
    n2, A2 = g6tools.decode_d6(s)
    assert n1 == n2 and A1 == A2, "the two digraph6 decoders disagree on %r" % s
    assert encode_d6_mine(A1) == s == g6tools.encode_d6(n1, A1), "round trip failed for %r" % s
    return n1, A1


def encode_checked(A):
    s = encode_d6_mine(A)
    assert s == g6tools.encode_d6(len(A), A), "the two digraph6 encoders disagree"
    n, B = decode_checked(s)
    assert B == A
    return s


def classify_pairs(A):
    n = len(A)
    bidir, single, missing = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] and A[j][i]:
                bidir.append((i, j))
            elif A[i][j]:
                single.append((i, j))
            elif A[j][i]:
                single.append((j, i))
            else:
                missing.append((i, j))
    arcs = sum(map(sum, A))
    assert arcs == 2 * len(bidir) + len(single)
    return bidir, single, missing, arcs


# ----------------------------------------------------------------------------------------------- exact roots
def poly_eval(coeffs_hi_to_lo, x):
    v = F(0)
    for c in coeffs_hi_to_lo:
        v = v * x + c
    return v


def sign(v):
    return (v > 0) - (v < 0)


def narrow_root(coeffs_hi_to_lo, lo, hi, width=F(1, 10 ** 9)):
    """Exact bisection of an isolating interval [lo, hi] of the stored minimal polynomial.

    The certificate states that [lo, hi] contains exactly one root (Sturm count 1) of this irreducible integer
    polynomial.  A root of an irreducible polynomial is simple, so the polynomial changes sign across it; we
    check that sign change at the stored ends and then bisect exactly.  Returns (lo', hi', status): the narrowed
    interval, still containing the same root, or the stored interval unchanged with status 'stored' if the
    sign check fails (it never does for the data used here, but the drawing must never be longer than certified).
    """
    lo, hi = F(lo), F(hi)
    coeffs = [F(c) for c in coeffs_hi_to_lo]
    slo, shi = sign(poly_eval(coeffs, lo)), sign(poly_eval(coeffs, hi))
    if slo == 0:
        return lo, lo, "exact_rational_root"
    if shi == 0:
        return hi, hi, "exact_rational_root"
    if slo == shi:
        return lo, hi, "stored"
    while hi - lo > width:
        mid = (lo + hi) / 2
        sm = sign(poly_eval(coeffs, mid))
        if sm == 0:
            return mid, mid, "exact_rational_root"
        if sm == slo:
            lo = mid
        else:
            hi = mid
    return lo, hi, "narrowed"


def brackets_root(coeffs_hi_to_lo, a, b):
    """True if [a, b] is an exact rational root (a == b, p(a) == 0) or the polynomial changes sign across [a, b]."""
    c = [F(x) for x in coeffs_hi_to_lo]
    if a == b:
        return poly_eval(c, a) == 0
    return sign(poly_eval(c, a)) * sign(poly_eval(c, b)) == -1


def fstr(x):
    return str(F(x))


# ----------------------------------------------------------------------------------------------- drawing helpers
def draw_digraph(ax, A, label_math, d6, arcs, node_r=0.13):
    """Circular layout; vertex k at angle 90 deg - 360 k/n (clockwise from the top), labelled k."""
    import math
    n = len(A)
    pos = []
    for k in range(n):
        ang = math.pi / 2 - 2 * math.pi * k / n
        pos.append((math.cos(ang), math.sin(ang)))
    bidir, single, missing, arcs2 = classify_pairs(A)
    assert arcs2 == arcs

    def trim(p, q, pad):
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy)
        ux, uy = dx / L, dy / L
        return (p[0] + ux * pad, p[1] + uy * pad), (q[0] - ux * pad, q[1] - uy * pad)

    pad = node_r + 0.015
    for (i, j) in missing:
        p, q = trim(pos[i], pos[j], pad)
        ax.plot([p[0], q[0]], [p[1], q[1]], color=MISSING, lw=0.7, ls=(0, (2.5, 2.2)), solid_capstyle="butt",
                zorder=1)
    for (i, j) in bidir:
        p, q = trim(pos[i], pos[j], pad)
        ax.plot([p[0], q[0]], [p[1], q[1]], color=BLUE, lw=2.2, solid_capstyle="butt", zorder=2)
    for (i, j) in single:  # i -> j : head at the receiver j
        p, q = trim(pos[i], pos[j], pad)
        ax.annotate("", xy=q, xytext=p,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.18", lw=0.8, color=INK,
                                    shrinkA=0, shrinkB=0, mutation_scale=8, joinstyle="miter"), zorder=3)
    for k, (x, y) in enumerate(pos):
        ax.add_patch(Circle((x, y), node_r, facecolor=SURFACE, edgecolor=INK, lw=0.7, zorder=4))
        ax.text(x, y, str(k), ha="center", va="center", fontsize=FS_BASE, color=INK, zorder=5)
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.32, 1.32)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(label_math + ", %d arcs" % arcs, fontsize=FS_TITLE, color=INK, pad=15)
    ax.text(0, 1.30, d6, ha="center", va="bottom", fontsize=FS_BASE, color=INK, family="monospace")
    return {
        "digraph6": d6, "n": n, "arcs": arcs,
        "vertex_order_clockwise_from_top": list(range(n)),
        "adjacency_rows": ["".join(map(str, row)) for row in A],
        "bidirected_pairs": [list(p) for p in bidir],
        "single_arcs_sender_to_receiver": [list(p) for p in single],
        "non_adjacent_pairs": [list(p) for p in missing],
        "counts": {"bidirected_pairs": len(bidir), "single_arcs": len(single), "non_adjacent_pairs": len(missing)},
    }




def line_legend(fig, y, x0=0.5):
    """Manual legend of the three line kinds (arrowheads cannot be shown by the automatic legend)."""
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_zorder(-1)
    items = [("bidirected pair", "bidir"), ("single arc (head at the receiver)", "single"),
             ("non-adjacent pair", "missing")]
    widths = [0.18, 0.32, 0.20]
    total = sum(widths)
    x = x0 - total / 2
    for (txt, kind), w in zip(items, widths):
        xa, xb = x, x + 0.045
        if kind == "bidir":
            ax.plot([xa, xb], [y, y], color=BLUE, lw=2.2, solid_capstyle="butt", transform=ax.transAxes)
        elif kind == "single":
            ax.annotate("", xy=(xb, y), xytext=(xa, y), xycoords="axes fraction", textcoords="axes fraction",
                        arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.18", lw=0.8, color=INK,
                                        shrinkA=0, shrinkB=0, mutation_scale=8))
        else:
            ax.plot([xa, xb], [y, y], color=MISSING, lw=0.7, ls=(0, (2.5, 2.2)), transform=ax.transAxes)
        ax.text(xb + 0.008, y, txt, va="center", ha="left", fontsize=FS_SMALL, color=INK2, transform=ax.transAxes)
        x += w


def fmt_int(k):
    return "{:,}".format(k)


# =============================================================================================== data loading
def load_json(*parts):
    with open(P(*parts)) as f:
        return json.load(f)


def load_jsonl(*parts):
    with open(P(*parts)) as f:
        return [json.loads(line) for line in f if line.strip()]


# =============================================================================================== fig 1
def fig1():
    gen = load_json("certificates", "stage0_geneson_2607.21396.json")
    A_g = gen["A_ours"]
    d6_g = encode_checked(A_g)
    assert d6_g == "&EZjtwWE", d6_g  # the project's digraph6 form of Geneson's graph
    n6 = load_json("certificates", "n6_summary.json")
    classes = load_json("certificates", "n6_survivor_equivalence_classes.json")["classes_with_hits"]
    hits = {h["idx"]: h for h in n6["hits"]}
    assert len(hits) == 3 and n6["limit"]["limit_stable_noncomplete"] == 3
    geneson_idx = [i for i, h in hits.items() if h["features"]["isomorphic_to_geneson"]]
    assert len(geneson_idx) == 1
    gclass = [c for c in classes if geneson_idx[0] in c]
    assert len(gclass) == 1
    outside = sorted(i for i in hits if i not in gclass[0])
    assert outside == [1387528], outside  # H_3 of the paper: the counterexample outside Geneson's twin class
    h3 = hits[1387528]
    A_h = h3["A"]
    d6_h = encode_checked(A_h)
    assert h3["features"]["edges"] == sum(map(sum, A_h)) == 18

    fig = plt.figure(figsize=(TEXTWIDTH_IN, 3.65))
    axs = [fig.add_axes([0.03, 0.14, 0.44, 0.73]), fig.add_axes([0.53, 0.14, 0.44, 0.73])]
    rec_g = draw_digraph(axs[0], A_g, "Geneson's graph (stage-0 certificate)", d6_g, 19)
    rec_h = draw_digraph(axs[1], A_h, r"$H_3$ (index 1,387,528, outside Geneson's twin class)", d6_h, 18)
    line_legend(fig, 0.055)
    out = os.path.join(HERE, "fig1_counterexamples.pdf")
    fig.savefig(out)
    plt.close(fig)
    DATA["figures"]["fig1_counterexamples.pdf"] = {
        "sources": ["certificates/stage0_geneson_2607.21396.json", "certificates/n6_summary.json",
                    "certificates/n6_survivor_equivalence_classes.json"],
        "layout": "circular, vertex k at angle 90deg - 360k/n (clockwise from the top), labelled k in digraph6 order",
        "line_kinds": {"bidirected_pair": "thick solid, no heads", "single_arc": "thin solid, arrowhead at receiver",
                       "non_adjacent_pair": "light dashed"},
        "panels": [dict(rec_g, source="certificates/stage0_geneson_2607.21396.json (A_ours)"),
                   dict(rec_h, source="certificates/n6_summary.json hits[idx=1387528].A",
                        n6_summary_idx=1387528, n6_summary_key=h3["key"],
                        twin_class_note="classes_with_hits = %s; 1387528 is the hit outside the class of the "
                                        "Geneson-isomorphic hit 1000255" % json.dumps(classes))],
        "decoder_cross_check": "own decoder and n8/k8decomp/g6tools.py decode_d6 agree on every string; encoders agree",
    }
    return out


# =============================================================================================== fig 2
def collect_intervals():
    """One (lo_inner, hi_inner) rational pair per limit counterexample, plus provenance, for n = 6, 7, 8."""
    out = {}
    # ---- n = 6: three hits in n6_summary.json
    n6 = load_json("certificates", "n6_summary.json")
    rows = []
    for h in n6["hits"]:
        ivs = h["limit_stable_intervals"]
        assert len(ivs) == 1
        rows.append(interval_from_cert(ivs[0], "idx%d" % h["idx"]))
    out[6] = rows
    # ---- n = 7: 123 records
    rows = []
    for r in load_jsonl("certificates", "COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl"):
        ivs = r["cert"]["stable_intervals"]
        assert len(ivs) == 1 and r["verdict"] == "stable_somewhere" and not r["complete"]
        rows.append(interval_from_cert(ivs[0], r["d6"]))
    assert len(rows) == 123
    out[7] = rows
    # ---- n = 8: 9,133 rows of n8_summary.json (refined isolating intervals), cross-checked with the HITS file
    s8 = load_json("certificates", "n8_summary.json")
    hits8 = {r["d6"]: r for r in load_jsonl("n8", "n8_limit_HITS.jsonl")}
    assert len(hits8) == 9133 == s8["hits"] == len(s8["rows"])
    rows = []
    for r in s8["rows"]:
        h = hits8[r["d6"]]
        iv = h["cert"]["stable_intervals"]
        assert len(iv) == 1
        iv = iv[0]
        hi_stored = [F(x) for x in iv["hi"]["q_interval"]]
        hi_ref = [F(x) for x in r["hi_q_interval_refined"]]
        assert hi_stored == [F(x) for x in r["hi_q_interval_stored"]]
        assert hi_stored[0] <= hi_ref[0] <= hi_ref[1] <= hi_stored[1] and hi_ref[1] - hi_ref[0] < F(1, 10 ** 9)
        # the refined interval must bracket the same simple root of the stored polynomial
        c = iv["hi"]["q_minpoly_hi_to_lo"]
        assert brackets_root(c, hi_ref[0], hi_ref[1]), r["d6"]
        if r["lo"]["kind"] == "beta_min":
            assert iv["lo"]["kind"] == "beta_min" and F(iv["lo"]["q"]) == 1 == F(r["lo"]["q"])
            lo_stored = [F(1), F(1)]
            lo_ref = [F(1), F(1)]
            lo_kind = "beta_min"
        else:
            assert iv["lo"]["kind"] == "root"
            lo_stored = [F(x) for x in iv["lo"]["q_interval"]]
            lo_ref = [F(x) for x in r["lo"]["refined"]]
            assert lo_stored[0] <= lo_ref[0] <= lo_ref[1] <= lo_stored[1] and lo_ref[1] - lo_ref[0] < F(1, 10 ** 9)
            c = iv["lo"]["q_minpoly_hi_to_lo"]
            assert c == r["lo"]["q_minpoly_hi_to_lo"]
            assert brackets_root(c, lo_ref[0], lo_ref[1]), r["d6"]
            lo_kind = "root"
        assert lo_ref[1] < hi_ref[0]
        rows.append({"label": r["d6"], "arcs": r["arcs"], "lo_kind": lo_kind, "hi_vanishing": iv["hi"]["vanishing"],
                     "lo_stored": [fstr(x) for x in lo_stored], "hi_stored": [fstr(x) for x in hi_stored],
                     "lo_refined": [fstr(x) for x in lo_ref], "hi_refined": [fstr(x) for x in hi_ref],
                     "refinement": "certificates/n8_summary.json rows[].lo.refined / hi_q_interval_refined",
                     "drawn_lo": lo_ref[1], "drawn_hi": hi_ref[0]})
    out[8] = rows
    return out


def interval_from_cert(iv, label):
    """n = 6 / n = 7 certificate interval -> narrowed inner bounds (the stored isolating intervals are coarse)."""
    hi = iv["hi"]
    assert hi["kind"] == "root"
    hi_stored = [F(x) for x in hi["q_interval"]]
    a, b, st = narrow_root(hi["q_minpoly_hi_to_lo"], hi_stored[0], hi_stored[1])
    assert st == "narrowed", (label, st)
    if iv["lo"]["kind"] == "beta_min":
        assert F(iv["lo"]["q"]) == 1
        lo_stored = [F(1), F(1)]
        la, lb, lst = F(1), F(1), "beta_min"
    else:
        lo = iv["lo"]
        assert lo["kind"] == "root"
        lo_stored = [F(x) for x in lo["q_interval"]]
        la, lb, lst = narrow_root(lo["q_minpoly_hi_to_lo"], lo_stored[0], lo_stored[1])
        assert lst == "narrowed", (label, lst)
    assert lb < a, (label, lb, a)
    return {"label": label, "lo_kind": iv["lo"]["kind"], "hi_vanishing": hi["vanishing"],
            "lo_stored": [fstr(x) for x in lo_stored], "hi_stored": [fstr(x) for x in hi_stored],
            "lo_refined": [fstr(la), fstr(lb)], "hi_refined": [fstr(a), fstr(b)],
            "refinement": "exact bisection of the stored minimal polynomial inside the stored isolating interval "
                          "(width < 1e-9), done by make_figures.py",
            "drawn_lo": lb, "drawn_hi": a}


def fig2():
    ivs = collect_intervals()
    q6 = load_json("certificates", "qstar", "qstar_n6.json")["n6_full_support"]["aggregate"]["q_star"]
    q7 = load_json("certificates", "qstar", "qstar_n7limit.json")["n7_limit_CANDIDATES_lower_bound"]["max_beta_c_minus_1"]
    q8 = load_json("certificates", "n8_summary.json")["q_star_limit"]
    qstar = {
        6: {"minpoly_hi_to_lo": q6["minpoly_hi_to_lo"], "isolating_interval": q6["isolating_interval"],
            "label": r"$q_6^*\approx1.1664$"},
        7: {"minpoly_hi_to_lo": q7["minpoly_hi_to_lo"], "isolating_interval": q7["isolating_interval"],
            "label": r"$q_7^*=\rho\approx1.3247$"},
        8: {"minpoly_hi_to_lo": q8["minpoly_hi_to_lo"], "isolating_interval": q8["isolating_interval"],
            "label": r"$q_8^*(\infty)\approx1.7107$"},
    }
    assert q7["minpoly_hi_to_lo"] == [1, 0, -1, -1]  # the plastic number
    assert q8["degree"] == 15 and q6["degree"] == 6
    # q* must equal the largest drawn upper end's root: check the isolating intervals overlap the max row
    for n in (6, 7, 8):
        a, b = (F(x) for x in qstar[n]["isolating_interval"])
        mx = max(ivs[n], key=lambda r: r["drawn_hi"])
        assert F(mx["hi_refined"][0]) <= b and a <= F(mx["hi_refined"][1]), n
        # narrow the stored isolating interval of q_n^* exactly (same bisection as the rows) so that the tick
        # coincides with the drawn end of the record rows instead of sitting at the midpoint of a coarse interval
        na, nb, st = narrow_root(qstar[n]["minpoly_hi_to_lo"], a, b)
        assert st == "narrowed" and a <= na < nb <= b and nb - na < F(1, 10 ** 9), (n, st)
        assert F(mx["hi_refined"][0]) <= nb and na <= F(mx["hi_refined"][1]), n
        qstar[n]["narrowed_interval"] = [fstr(na), fstr(nb)]
        qstar[n]["tick"] = (na + nb) / 2
    # rows whose interval starts above q = 1 (lower end of kind root) show as white gaps inside a stack
    gaps = {n: sum(1 for r in ivs[n] if r["lo_kind"] == "root") for n in (6, 7, 8)}
    assert gaps == {6: 0, 7: 6, 8: 971}, gaps

    FIG_H = 4.1
    AX_BOX = [0.075, 0.11, 0.9, 0.85]
    YLIM = (-0.25, 6.05)
    fig = plt.figure(figsize=(TEXTWIDTH_IN, FIG_H))
    ax = fig.add_axes(AX_BOX)
    xlo, xhi = 1.0, 1.8
    # strip geometry (data units on a dimensionless y axis)
    strips = {8: (0.0, 2.6), 7: (3.0, 4.6), 6: (5.0, 5.6)}
    for n in (6, 7, 8):
        rows = sorted(ivs[n], key=lambda r: (r["drawn_hi"], r["drawn_lo"]))
        y0, y1 = strips[n]
        k = len(rows)
        ys = [y0 + (y1 - y0) * (i + 0.5) / k for i in range(k)]
        if n == 8:
            # dense stack rendered as its envelope: row i (from the bottom) spans the i-th smallest upper end.
            # step="mid" gives every row its own horizontal band (no interpolation between neighbouring rows),
            # so a row whose interval starts above q = 1 leaves a rectangular white gap at its own height.
            los = [float(r["drawn_lo"]) for r in rows]
            his = [float(r["drawn_hi"]) for r in rows]
            ax.fill_betweenx(ys, los, his, color=BLUE, lw=0, zorder=2, step="mid")
        else:
            segs = [[(float(r["drawn_lo"]), y), (float(r["drawn_hi"]), y)] for r, y in zip(rows, ys)]
            pitch_pt = (y1 - y0) / k * (AX_BOX[3] * FIG_H * 72 / (YLIM[1] - YLIM[0]))  # true row pitch in pt
            lw = min(1.6, max(0.5, pitch_pt * 0.8))
            ax.add_collection(LineCollection(segs, colors=[BLUE], linewidths=lw, capstyle="butt", zorder=2))
        # strip baseline & label
        ax.plot([xlo, xhi], [y0 - 0.08, y0 - 0.08], color=AXIS, lw=0.5, zorder=1)
        ax.text(xlo + 0.004, y1 + 0.02, r"$n=%d$: %s interval%s" % (n, fmt_int(k), "" if k == 1 else "s"),
                ha="left", va="bottom", fontsize=FS_BASE, color=INK)
        # q_n^* tick at the midpoint of the narrowed isolating interval (width < 1e-9; float for placement only)
        xq = float(qstar[n]["tick"])
        ax.plot([xq, xq], [y0 - 0.08, y1 + 0.02], color=INK, lw=0.7, zorder=3)
        if n == 8:  # label to the left of the tick: keeps it clear of the right margin whatever its length
            ax.text(xq - 0.008, y1 - 0.02, qstar[n]["label"], ha="right", va="top", fontsize=FS_BASE, color=INK)
        else:
            ax.text(xq + 0.008, y1 - 0.02 if n != 6 else y1 + 0.02, qstar[n]["label"], ha="left",
                    va="top" if n != 6 else "bottom", fontsize=FS_BASE, color=INK)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(*YLIM)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_position(("data", -0.2))
    ax.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
    ax.set_xlabel(r"$q=\delta(1-\varepsilon)/\varepsilon$  (limit $\varepsilon\to0$: $q=\beta-1$)", fontsize=FS_BASE)
    ax.tick_params(axis="x", length=2.5)
    out = os.path.join(HERE, "fig2_thresholds.pdf")
    fig.savefig(out)
    plt.close(fig)

    def rows_json(rows):
        return [{k: (fstr(v) if isinstance(v, F) else v) for k, v in r.items()} for r in rows]

    DATA["figures"]["fig2_thresholds.pdf"] = {
        "sources": ["certificates/n6_summary.json", "certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl",
                    "certificates/n8_summary.json", "n8/n8_limit_HITS.jsonl", "certificates/qstar/qstar_n6.json",
                    "certificates/qstar/qstar_n7limit.json"],
        "q_axis": [1, 1.8],
        "rule": "each interval is drawn from the UPPER bound of the isolating interval of its lower end (or q = 1 "
                "when the interval starts at beta = 2) to the LOWER bound of the isolating interval of its upper "
                "end, so nothing is drawn longer than certified; rows sorted by upper end (then lower end), bottom "
                "to top",
        "rendering": {"n6": "LineCollection of 3 thin segments", "n7": "LineCollection of 123 thin segments",
                      "n8": "the stack of 9,133 rows rendered as its envelope (fill_betweenx with step='mid': "
                            "one horizontal band per row between its drawn lower and upper end, no interpolation "
                            "between neighbouring rows), which is what 9,133 segments of sub-point thickness "
                            "would show"},
        "white_gaps": {"meaning": "a white gap inside a stack is a row whose interval starts above q = 1 (lower "
                                  "end of kind root); such rows are grouped because rows with equal upper ends "
                                  "are sub-sorted by lower end",
                       "rows_with_root_lower_end": {str(n): gaps[n] for n in (6, 7, 8)}},
        "counts": {str(n): len(ivs[n]) for n in (6, 7, 8)},
        "q_star_ticks": {str(n): {"minpoly_hi_to_lo": qstar[n]["minpoly_hi_to_lo"],
                                  "isolating_interval": qstar[n]["isolating_interval"],
                                  "narrowed_interval": qstar[n]["narrowed_interval"],
                                  "drawn_at": "midpoint of the narrowed isolating interval (exact bisection of "
                                              "the stored minimal polynomial inside the stored isolating "
                                              "interval to width < 1e-9; float for placement only)",
                                  "label": qstar[n]["label"]} for n in (6, 7, 8)},
        "intervals": {str(n): rows_json(ivs[n]) for n in (6, 7, 8)},
    }
    return out


# =============================================================================================== fig 3
def fig3():
    tot = open(P("n8", "n8_prefilter_TOTAL.txt")).read().split()
    tot = {tot[i]: int(tot[i + 1]) for i in range(0, len(tot), 2)}
    tally = load_json("n8", "n8_limit_tally.json")
    s8 = load_json("certificates", "n8_summary.json")["counts"]
    assert tot["total"] == 801558826511 == s8["digraphs_read"]
    assert tot["domfree"] == 113697991934 == s8["domination_free"]
    assert tot["survivors"] == 417951273 == s8["spectral_survivors"]
    assert tot["boundary"] == 417603634 == s8["boundary"]
    assert tot["survivors"] - tot["boundary"] == 347639 == tally["S"] == s8["strict"]
    assert tally["stable_noncomplete"] == 9133 == s8["candidates"] and tally["complete"] == 1
    assert tally["S"] == tally["never_positivity"] + tally["never_full"] + tally["stable_somewhere"] + tally["escalate"]
    stages = [
        ("digraphs read (at least 29 arcs)", tot["total"], ""),
        ("domination-free", tot["domfree"], ""),
        ("spectral survivors", tot["survivors"], "of which %s boundary lines, dropped" % fmt_int(tot["boundary"])),
        ("strict survivors (exact limit engine)", 347639, ""),
        ("limit counterexamples", 9133, "plus $K_8$"),
    ]
    fig = plt.figure(figsize=(TEXTWIDTH_IN, 2.45))
    ax = fig.add_axes([0.33, 0.19, 0.53, 0.78])
    import math
    ys = list(range(len(stages)))[::-1]
    for (name, cnt, note), y, col in zip(stages, ys, RAMP):
        ax.barh(y, math.log10(cnt), left=0, height=0.42, color=col, lw=0, zorder=2)
        ax.text(math.log10(cnt) + 0.12, y, fmt_int(cnt), va="center", ha="left", fontsize=FS_BASE, color=INK,
                zorder=3)
        if note:
            ax.text(math.log10(cnt) + 0.12, y - 0.36, note, va="center", ha="left", fontsize=FS_MIN, color=INK2,
                    zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([s[0] for s in stages], fontsize=FS_BASE, color=INK)
    ax.set_xlim(0, 14.3)
    ax.set_xticks(range(0, 13, 3))
    ax.set_xticklabels([r"$10^{%d}$" % k for k in range(0, 13, 3)], fontsize=FS_SMALL)
    ax.set_xlabel("count (logarithmic scale)", fontsize=FS_BASE)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(stages) - 0.4)
    out = os.path.join(HERE, "fig3_funnel.pdf")
    fig.savefig(out)
    plt.close(fig)
    DATA["figures"]["fig3_funnel.pdf"] = {
        "sources": ["n8/n8_prefilter_TOTAL.txt", "n8/n8_limit_tally.json", "certificates/n8_summary.json (counts)"],
        "bars_top_to_bottom": [{"stage": s, "count": c, "annotation": n} for s, c, n in stages],
        "bar_length": "log10(count), bars start at 1",
        "colour": "single-hue ordinal blue ramp, light to dark",
    }
    return out


# =============================================================================================== fig 4
def fig4():
    hits8 = {r["d6"]: r for r in load_jsonl("n8", "n8_limit_HITS.jsonl")}
    s8 = load_json("certificates", "n8_summary.json")
    rec, sparse = "&GMx|\\Aw^qz]w", "&GM?WUZeEie^O"
    assert rec in s8["q_star_limit"]["attained_by"] and s8["q_star_limit"]["arcs_of_attaining_graphs"] == [40]
    assert sparse in s8["min_arcs_hits"] and s8["edge_counts"][0] == [29, 2]
    rows = {r["d6"]: r for r in s8["rows"]}
    n1, A1 = decode_checked(rec)
    n2, A2 = decode_checked(sparse)
    assert n1 == n2 == 8
    assert sum(map(sum, A1)) == 40 == hits8[rec]["edges"] == rows[rec]["arcs"]
    assert sum(map(sum, A2)) == 29 == hits8[sparse]["edges"] == rows[sparse]["arcs"]
    fig = plt.figure(figsize=(TEXTWIDTH_IN, 3.85))
    axs = [fig.add_axes([0.03, 0.13, 0.44, 0.74]), fig.add_axes([0.53, 0.13, 0.44, 0.74])]
    r1 = draw_digraph(axs[0], A1, r"record upper end $q_8^*(\infty)\approx1.7107$", rec, 40, node_r=0.12)
    r2 = draw_digraph(axs[1], A2, r"edge bound attained, upper end $\approx1.0131$", sparse, 29, node_r=0.12)
    assert r1["counts"]["non_adjacent_pairs"] == 0 and r1["counts"]["bidirected_pairs"] == 12
    assert r2["counts"]["non_adjacent_pairs"] == 2 and r2["counts"]["bidirected_pairs"] == 3
    line_legend(fig, 0.05)
    out = os.path.join(HERE, "fig4_n8_graphs.pdf")
    fig.savefig(out)
    plt.close(fig)
    DATA["figures"]["fig4_n8_graphs.pdf"] = {
        "sources": ["n8/n8_limit_HITS.jsonl", "certificates/n8_summary.json (q_star_limit, min_arcs_hits, rows)"],
        "layout": "circular, vertex k at angle 90deg - 360k/8 (clockwise from the top), labelled k in digraph6 order",
        "panels": [dict(r1, upper_end_refined=rows[rec]["hi_q_interval_refined"], lower_end=rows[rec]["lo"],
                        tie="&GVrm]A{^pZ]w (same class, same root)"),
                   dict(r2, upper_end_refined=rows[sparse]["hi_q_interval_refined"], lower_end=rows[sparse]["lo"])],
        "decoder_cross_check": "own decoder and n8/k8decomp/g6tools.py decode_d6 agree; arc counts equal the "
                               "'edges' field of the HITS record and 'arcs' of the summary row",
    }
    return out


# =============================================================================================== fig 5
def fig5():
    s8 = load_json("certificates", "n8_summary.json")
    hits = load_jsonl("n8", "n8_limit_HITS.jsonl")
    assert len(hits) == 9133
    # --- arcs: recount from the HITS file and assert against the summary
    arc_counts = {}
    for h in hits:
        n, A = decode_checked(h["d6"]) if h["edges"] in (29, 45) else g6tools.decode_d6(h["d6"])
        assert sum(map(sum, A)) == h["edges"]
        arc_counts[h["edges"]] = arc_counts.get(h["edges"], 0) + 1
    assert sorted(arc_counts.items()) == [tuple(x) for x in s8["edge_counts"]]
    assert sum(arc_counts.values()) == 9133 and min(arc_counts) == 29 and max(arc_counts) == 45
    # --- upper ends: certified upper bound of the refined isolating interval, exact bins of width 1/20
    edges = [F(1) + F(k, 20) for k in range(0, 16)]  # 1.00, 1.05, ..., 1.75
    bins = [0] * 15
    straddle = 0
    for r in s8["rows"]:
        a, b = (F(x) for x in r["hi_q_interval_refined"])
        assert F(1) < a <= b <= F(7, 4)  # a == b only for the two exact rational ends q = 5/3
        ib = next(i for i in range(15) if edges[i] < b <= edges[i + 1])
        ia = next(i for i in range(15) if edges[i] < a <= edges[i + 1])
        if ia != ib:
            straddle += 1
        bins[ib] += 1
    assert sum(bins) == 9133
    # cross-check against the summary's own histogram where the bin edges coincide
    H = {d["q_in"]: d["count"] for d in s8["q_data"]["upper_end_histogram"]}
    assert bins[0] == H["(10001/10000, 1001/1000]"] + H["(1001/1000, 101/100]"] + H["(101/100, 21/20]"] + H["(1, 10001/10000]"]
    assert bins[1] == H["(21/20, 11/10]"]
    assert bins[2] + bins[3] == H["(11/10, 6/5]"]
    assert bins[4] + bins[5] == H["(6/5, 13/10]"]
    assert bins[6] + bins[7] == H["(13/10, 7/5]"]
    assert bins[8] + bins[9] == H["(7/5, 3/2]"]
    assert bins[10] + bins[11] == H["(3/2, 8/5]"]
    assert bins[12] + bins[13] == H["(8/5, 17/10]"]
    assert bins[14] == H["(17/10, 7/4]"]

    fig = plt.figure(figsize=(TEXTWIDTH_IN, 2.7))
    axL = fig.add_axes([0.08, 0.19, 0.39, 0.70])
    axR = fig.add_axes([0.585, 0.19, 0.39, 0.70])
    gap_pt = 1.5  # ~2 px at 96 dpi

    def bars(ax, xs, hs, width_units, labels):
        ax_w_pt = ax.get_position().width * TEXTWIDTH_IN * 72
        span = ax.get_xlim()[1] - ax.get_xlim()[0]
        gap_units = gap_pt / ax_w_pt * span
        w = width_units - gap_units
        ax.bar(xs, hs, width=w, color=BLUE, lw=0, zorder=2, align="center")
        ymax = max(hs)
        for x, h, lab in zip(xs, hs, labels):
            ax.text(x, h + ymax * 0.015, lab, ha="center", va="bottom", rotation=90, fontsize=FS_MIN, color=INK2,
                    zorder=3)
        ax.set_ylim(0, ymax * 1.30)
        for s in ("left", "right", "top"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)

    ks = list(range(29, 46))
    axL.set_xlim(28.3, 45.7)
    bars(axL, ks, [arc_counts[k] for k in ks], 1.0, [fmt_int(arc_counts[k]) for k in ks])
    axL.set_xticks([29, 31, 33, 35, 37, 39, 41, 43, 45])
    axL.set_yticks([0, 500, 1000, 1500])
    axL.set_xlabel("number of arcs (of 56)", fontsize=FS_BASE)
    axL.set_title("counterexamples by arc count", fontsize=FS_TITLE, color=INK, loc="left")

    centres = [float((edges[i] + edges[i + 1]) / 2) for i in range(15)]
    axR.set_xlim(1.0, 1.75)
    bars(axR, centres, bins, 0.05, [fmt_int(b) for b in bins])
    axR.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
    axR.set_yticks([0, 1000, 2000, 3000])
    axR.set_xlabel("upper end of the stable $q$-interval (bins of width 0.05)", fontsize=FS_BASE)
    axR.set_title("counterexamples by upper end of the $q$-interval", fontsize=FS_TITLE, color=INK, loc="left")
    out = os.path.join(HERE, "fig5_histograms.pdf")
    fig.savefig(out)
    plt.close(fig)
    DATA["figures"]["fig5_histograms.pdf"] = {
        "sources": ["n8/n8_limit_HITS.jsonl (arc counts recounted from the decoded adjacency and the edges field)",
                    "certificates/n8_summary.json (edge_counts asserted equal; rows[].hi_q_interval_refined for "
                    "the upper ends; q_data.upper_end_histogram asserted consistent)"],
        "arc_count_histogram": {str(k): arc_counts[k] for k in ks},
        "upper_end_histogram": [{"bin": "(%s, %s]" % (fstr(edges[i]), fstr(edges[i + 1])), "count": bins[i]}
                                for i in range(15)],
        "bin_rule": "a record is binned by the certified UPPER bound of the refined isolating interval of its "
                    "upper end (exact rational comparison); records whose lower refined bound falls in a "
                    "different bin: %d" % straddle,
        "bar_gap": "1.5 pt (about 2 px) of surface between adjacent bars",
        "count_labels": "printed above every bar at %.1f pt (the figure is placed 1:1 at the text width)" % FS_MIN,
    }
    return out


def main():
    outs = [fig1(), fig2(), fig3(), fig4(), fig5()]
    with open(os.path.join(HERE, "figures_data.json"), "w") as f:
        json.dump(DATA, f, indent=1)
    for o in outs:
        print("wrote", os.path.relpath(o, ROOT), os.path.getsize(o), "bytes")
    print("wrote notes/paper/figures/figures_data.json")


if __name__ == "__main__":
    main()
