import json
res = json.load(open('results.json'))
lifts = json.load(open('lift_123_to_8.json'))
sections = open('sections.md').read()
recipe = open('recipe_check.py').read()
hdr = r'''# HANDCHECK kit: six n = 8 limit counterexample candidates, exact, hand-verifiable

Agent `handcheck_kit`, second verification campaign (`notes/verify_n8_2/handcheck_kit/`), 2026-10-28.
Everything below was computed by my own exact implementation (`engine.py`: Python `fractions.Fraction` only,
Faddeev-LeVerrier characteristic polynomial cross-checked against determinants, Hurwitz minors, three-valued verdict),
and independently reproduced by the 10-line sympy recipe in section 2 (identical exact charpolys and minors on all six).
No floating point is used anywhere in a decision; decimals in parentheses are approximations for orientation only.

## 1. Definitions used (CLAUDE.md section 2)

* Digraph on n = 8 nodes, `A[i][j] = 1` iff there is an arc `i -> j` (row i = who i SENDS to).
* digraph6: `&` + chr(n + 63) + the n*n adjacency bits row by row, 6 bits per character, MSB first, char = value + 63.
* CTLN with theta = 1:  `W[i][j] = 0` if i = j;  `-1 + eps` if `A[j][i] = 1` (j sends to i, weak inhibition);  `-1 - delta` otherwise.
  Equivalently  **I - W = -delta I + (1 + delta) J - (eps + delta) A^T**  (J = all-ones): diagonal 1, off-diagonal `1 - eps` where `A[j][i] = 1`, else `1 + delta`.
* The full support {0,...,7} is a stable fixed point iff (1) `(I - W) x = 1` has a solution with every `x_i > 0`, and (3) every eigenvalue of
  `-I + W` has negative real part; condition (2) is vacuous for the full support.  (3) is decided by the Routh-Hurwitz criterion on the
  characteristic polynomial `det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8`: all eight leading principal minors `D1..D8` of the Hurwitz matrix
  `H[i][j] = c_{2j - i}` (1-indexed, c_0 = 1, c_k = 0 outside 0..8) positive <=> all roots in the open left half plane.
  A vanishing minor would be UNDECIDED (never happens below: all 48 minors are strictly positive).
* Legal CTLN: eps > 0, delta > 0, **q := delta (1 - eps) / eps > 1**.
* A non-complete digraph (some arc missing) whose full support is a stable positive fixed point at a legal point is a counterexample to the
  target-free clique conjecture (the full support of a non-complete graph is not a clique, hence not a target-free clique).
* Limit ray of the project: eps = 1/(t + 1), delta = (beta - 1) eps, so q = (beta - 1)(1 - eps) -> beta - 1 as t -> infinity.
  The stored `stable_intervals` are intervals in beta (and q = beta - 1) of the eps -> 0 limit object; the witness points below are at finite t.

## 2. How to check one of these yourself (10 lines)

**Plain Python + sympy (exact rationals)** -- save as `recipe_check.py`, run `python3 recipe_check.py '<digraph6>' <eps> <delta>` (quote the string: it contains `\` and `|`):

```python
''' + recipe + r'''```

**Sage**: `A = matrix(ZZ, 8, 8, bits)`; `IW = -delta*identity_matrix(8) + (1+delta)*ones_matrix(8) - (eps+delta)*A.transpose()` with `eps, delta` in `QQ`;
`x = IW.solve_right(vector(QQ, [1]*8))`; `c = (-IW).charpoly().coefficients(sparse=False)[::-1]`; build `H[i,j] = c[2j-i+1]` (1-indexed, out of range 0) and check `all(H[:k,:k].det() > 0 for k in 1..8)`.

**Mathematica**: `IW = -delta IdentityMatrix[8] + (1+delta) ConstantArray[1,{8,8}] - (eps+delta) Transpose[A]` (eps, delta exact `Rational`);
`x = LinearSolve[IW, ConstantArray[1,8]]`; `c = CoefficientList[CharacteristicPolynomial[-IW, s], s] // Reverse`; Hurwitz minors as above with `Det`.
`Positive[x]` and all minors positive is the claim. (A float `Eigenvalues[N[-IW, 50]]` is a fine sanity check but is not the decision.)

Every quantity in section 3 (A, eps, delta, I - W, x, c1..c8, D1..D8) can be regenerated this way in well under a minute per graph.

## 3. The six candidates

Summary (all six: non-complete, legal q > 1, x > 0 exactly, D1..D8 > 0 exactly => **stable**, no undecided minor):

| tag | digraph6 | arcs | hits index | eps | delta | q (exact) | q (approx) | verdict at point |
|---|---|---|---|---|---|---|---|---|
'''
rows = ''
for tag, r in res.items():
    rows += f"| {tag} | `{r['d6']}` | {r['arcs']} | {r['idx']} | {r['eps']} | {r['delta']} | {r['q']} | {float(__import__('fractions').Fraction(r['q'])):.6f} | positive, {r['stability']} |\n"
tail = r'''

### Appendix A: the 123 n = 7 limit counterexamples lifted to 8 nodes (closed-twin blow-up)

Closed-twin blow-up at vertex v: add a new vertex v' with exactly the in- and out-neighbourhood of v, plus the two arcs v -> v' and v' -> v
(v and v' become bidirected closed twins).  For each of the 123 records of `certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl` I built the
7 blow-ups, canonicalised them and the 9,133 hits with nauty `labelg` (`n7/labelg`, digraph6 in/out), and intersected.
Result (`lift_123_to_8.json`): **52 of the 123 parents have at least one blow-up in the hit list; 61 distinct hits are such blow-ups**
(63 (parent, vertex) pairs land in the hit list; two hits each arise from two different vertices of the same parent, no hit arises from two different parents).  This agrees with the independent count of the `completeness_twins` agent (61).
The remaining 9,072 hits are not single closed-twin lifts of an n = 7 limit counterexample.  None of the six candidates above is such a lift.

Note on the Geneson graph: Geneson's 6-node counterexample is `&EZjtwWE` (re-encoded by me from the Stage-0 certificate matrix; my engine reproduces the
Stage-0 verdicts: positive but unstable at eps = 1/10, 1/50; stable at 1/100, 1/1000 with delta = 29 eps/25).  Its 6 single closed-twin blow-ups
(7 nodes) are not among the 123, and its 21 distinct double blow-ups (8 nodes) are not among the 9,133 hits, so the fallback rule ("any candidate
with 34 arcs") was applied for C5.  (The hit list covers only digraphs whose underlying graph has 15..28 edges, passing the F1 domination-freeness
filter and the spectral prefilter, so absence from the list is not a statement about stability of these blow-ups.)

Lift table (parent index in the n = 7 file, parent digraph6, parent arcs -> blown-up vertex v, hit index, hit digraph6, hit arcs):

'''
for p in lifts:
    if p['lifts_in_hits']:
        tail += f"- n7 #{p['n7_index']} `{p['n7_d6']}` ({p['n7_arcs']} arcs) -> " + '; '.join(
            f"v={l['v']}: hit #{l['hit_index']} `{l['hit_d6']}` ({l['hit_arcs']} arcs)" for l in p['lifts_in_hits']) + '\n'
tail += r'''
### Appendix B: files in this directory

* `engine.py` -- the exact decision engine (own code).  `t0_selftest.py` -- self-test against Stage 0 (Geneson) and K4 / 3-cycle / toy Hurwitz cases.
* `build.py` -- selects the six candidates by the stated rules, runs the engine at the certified legal witness, writes `results.json` and `sections.md`.
* `refine.py` -- exact bisection of the stored isolating intervals (`interval_refine.json`), proving each witness beta is strictly inside the stored stable interval.
* `recipe_check.py` -- the 10-line sympy recipe; `lifts.py` -- Appendix A (`lift_123_to_8.json`, `lift_hit_indices.json`); `float_sanity.json` -- non-decision numpy eigenvalues.
* Re-run: `python3 -I t0_selftest.py; python3 -I lifts.py; python3 -I build.py; python3 -I refine.py; python3 -I build.py; python3 -I assemble.py` (about 1 minute total, single process).
'''
open('HANDCHECK.md', 'w').write(hdr + rows + '\n' + sections + tail)
print('written', len(hdr + rows + sections + tail), 'bytes')
