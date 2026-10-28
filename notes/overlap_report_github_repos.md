# Overlap report: two unrefereed GitHub repositories vs. the CTLN target-free-clique project

Date of retrieval: 2026-10-13. Retriever: `literature` agent.

## 0. Handling and trust

- Both repos were cloned (`git clone --depth 50`) into separate fresh directories under
  `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/repos/`.
  Repo 1 HEAD `e28a309` (2026-10-01, Ian Pitchford, single visible commit "Publish reviewed CTLN sharp-minimum candidate v1.1").
  Repo 2 HEAD `52ad000` (2026-08-08, Maximiliano Vazquez; 5 commits 2026-08-06..08).
- Nothing from either repo was executed (no scripts, notebooks, binaries, no Z3 on their queries). Files were read with
  `cat`/`sed`; PDFs with `pdftotext`; the one tarball (Repo 2) was listed first (1757 regular entries, no absolute paths,
  no `..`, no symlinks) and then extracted into its own fresh directory.
- Every number marked "MY CHECK" below was produced by scripts I wrote, kept OUTSIDE the clones, in
  `.../scratchpad/scripts/` and run as `python3 -I` from `/tmp`. They use `fractions.Fraction` and exact sympy only; no float
  enters a decision. Logs: `check_seed.log`, `grid.py` output (quoted below), `near_one.log`, `family_nd.log`,
  `repo2_counts.log`, `small_order_sanity.log`.
- The GitHub API (`gh api`) is NOT enabled for these repos in this session (HTTP 403), so CI run status, issues, forks, stars,
  releases and the Zenodo DOI page were NOT checked. I did not attach the repos with credentials because the task did not need it.

## 1. Headline findings (read this first)

**H1. A third-party claim that lands inside our gap, which I reproduced exactly. Per CLAUDE.md s3 this is a "claimed
counterexample": it should go to a human for hand-verification before anything is built on it.**
Repo 1 claims Geneson's 6-node graph is a stable full-support non-clique at `eps = 1/12, delta = 11/120`
(`q = delta(1-eps)/eps = 121/120 ~ 1.00833`), not only on Geneson's ray `delta = 29 eps/25` (`q -> 29/25`).
MY CHECK, in OUR convention (W built from A per CLAUDE.md s2; exact Fractions; full-support so condition 2 is vacuous):
legal (`eps = 1/12 < delta/(1+delta) = 11/131`); `x* = (2034120, 3972000, 3972000, 4806120, 922800, 7326120)/22726801 > 0`
(equals their vector); all six Hurwitz determinants of `det(zI + (I-W))` are strictly positive (none zero, so no
three-valued ambiguity); all 63 principal minors and all 192 uniform-input column-replacement determinants of `I-W` are nonzero
(globally nondegenerate). Their seed graph is, in our convention, EXACTLY the matrix you gave me (identity permutation).
This is the known Geneson graph at a NEW parameter point, so it is not a new topology; but it is a decided counterexample
at q within 0.84% of 1. Please re-run `ctln.analyze_support` on your own A at (1/12, 11/120) before relying on my script.

**H2. If their phase theorem (Thm 5.1) is right, the lower end of the n=6 gap is already closed for this graph.**
With `u = delta/eps - 1` and `t = (1-eps)/eps` one has `q = (1+u)(1-eps)` and legality `q>1 <=> t*u>1`.
Their claim: stable full-support iff `0 < u < u*` and `t > max(1/u, T(u))`, `u* = 0.166355027849...`. Hence every counterexample
of this graph has `1 < q < 1 + u* = 1.16635502784...`, and for `u <= 1/10` EVERY legal `(eps, delta)` works (their Cor 5.2), so
counterexamples exist for q arbitrarily close to 1. MY CHECK: exact stability at q-1 as small as 1e-22 (u = 1e-5, t just above 1/u)
and a 1,353-point exact grid over (u, t) with 0 mismatches against the claimed region (details in s2(c)). I did NOT re-prove the
theorem; I verified the polynomials it rests on and sampled it. Consequence for our plan (if confirmed): for this graph there is
no "q* below which counterexamples appear" (inf q = 1); the live n=6 question becomes (i) are there other 6-node graphs, and
(ii) how high in q can any counterexample reach (this graph stops at 1.1664; Geneson's sufficient bound is q >= 4 - 1.5 eps at n=6).

**H3. What is and is not scooped.**
- Repo 1 claims (all-legal-parameter, computer-assisted, exact-integer): (i) no nonclique full support at ANY order <= 5 (so
  our Stage 1 and the "supports of size <= 5" part of Stage 4 would be covered, for full support of the induced subgraph);
  (ii) the complete legal full-support phase region of Geneson's graph (our Stage 3 for that one graph);
  (iii) an infinite family at n = 5r+6 (r >= 0) with size-dependent parameters (n = 6, 11, 16, ...; NOT n = 7).
- Repo 2 claims only n = 5 full support (Z3, 875 classes); strictly weaker than Repo 1 (i) and worse-certified.
- NOT claimed by either: the complete n=6 map over all 1,540,944 digraphs (other 6-node counterexamples, or proof there are
  none); any n=7 statement; any classification of which graphs fail. Repo 1 says so explicitly (CLAIMS.md "Explicit exclusions";
  paper "Scope"). Our Stage 2/3/5 programme at n=6 and any n=7 work are not preempted by what is in these repos.
- Authorship: Repo 1 is "Scholarly creator: Anonymous", AI-assisted, "Evidence Press", maintainer Ian Pitchford; Repo 2 is
  Maximiliano Vazquez (GitHub ChickenNungets), apparently agent-driven (has an `AGENTS.md`). Neither is refereed.

**H4. Things in Repo 2 that are directly useful to us.** Its `docs/CONVENTIONS.md` uses OUR conventions exactly (source-major
`edge[u][v]=1 iff u->v`; `N[i][j]=1 iff j-/->i`; `q = delta(1-eps)/eps`; legal iff `q>1`), so its masks and nonedge matrices are
directly diff-able against ours. It quotes Curto's graphical-domination Definition 1.1 (and I confirmed the wording against the PDF
copy it bundles; see s4). Its `manifest.json` lists the 875 canonical masks + orbit sizes, and MY CHECK reproduced every filter
count, the 92,645 survivors, the 875 classes and every orbit size exactly from the stated definitions.

## 2. Repo 1: ipitchford/ctln-six-is-minimal (v1.1, 1 Oct 2026)

### 2(a) CLAIMS.md and STATUS.md, verbatim

Note the internal version inconsistency: CLAIMS.md says "Release: v1.0", STATUS.md/README/paper say v1.1 (HEAD is the single v1.1 commit; FINAL_CHECKS.md also says "Release v1.0 ... PDF: 14 pages" while the shipped PDF has 15 pages).

#### CLAIMS.md
~~~~
# Claim and evidence ledger

Release: v1.0, 1 October 2026. “Proved” below denotes the candidate proof and executed exact certificate, not completed peer review.

| ID | Claim | Exact evidence | Scope / status |
|---|---|---|---|
| C1 | No permitted positive-stable nonclique active matrix of order at most five exists at legal parameters. | Theorem 3.1; `small_order_obstructions.json`; `verify_small_order.py`; exact complete labelled orbit cover. | Computer-assisted proof. Covers the full legal continuum, not a parameter grid. No nondegeneracy needed for this algebraic statement. |
| C2 | Minimum stable nonclique support size is six. | C1 plus rational seed witness and all 255 global nondegeneracy determinant checks. | For local exponential stability in **nondegenerate** legal CTLNs, as in the source conjecture. |
| C3 | The seed graph is stable at epsilon=1/12, delta=11/120. | Integer matrix `10I+21N+110J`, exact equilibrium, quintic, four positive Hurwitz determinants, nonzero determinant enumerations. | Existing graph due to Geneson; new rational certificate in this candidate, not a new graph. |
| C4 | Stable nonclique full supports exist for every size 5r+6, integer r>=0, at explicit legal rational parameters. | Theorem 4.3; five-class quotient; full symbolic permission and Hurwitz polynomials; parity identity. | All-r theorem, not an extrapolation from the five regression sizes. Does not claim every integer size. |
| C5 | Global nondegeneracy of the family holds for all subsets and all sizes. | `D K = D I +(2D+2)N+2D^2 J = I (mod 2)`, with D odd; column-replacement argument. | Analytic all-size proof. The modular device is not claimed as a wholly new general method. |
| C6 | The seed's complete legal full-support stability region is 0<u<u_star and t>max(1/u,T(u)). | Theorem 5.1; coefficient bounds; unique roots of Hurwitz determinants; ordering t2<t3<T; exact rational isolators. | Necessary and sufficient for this **one graph**, all legal parameters. Does not assert that all other subsets are nondegenerate at every such point. |
| C7 | Every legal seed parameter with 1<delta/epsilon<=11/10 is stable. | Corollary 5.2; positive lower bound 2256/25 for the legal-boundary polynomial. | Whole cone, not just one parameter ray. |
| C8 | Outside-neuron activity is governed by the weighted missing-source inequality. | Proposition 6.1; exact substitution of the equilibrium into the outside input. | Direct reformulation of the standard off-neuron test. Not claimed as a previously unknown general criterion. |
| C9 | The rational seed has exactly ten inclusion-minimal missing-source sets that silence an outside neuron. | Exact 64-subset check of weights with common denominator 1613071. | Specific to beta=21/10; no equality cases at that beta. |
| C10 | All included verification scripts passed; nine intentionally invalid certificates were rejected. | `verification/verification_summary.json`, `regression_result.json`, script logs. | An executed software fact, not external peer review or a formal proof-assistant certificate. |

## New versus inherited

**Inherited and attributed:** the CTLN formulation and support criterion; target-free-clique sufficiency; classification through size four; Geneson's seed topology and its counterexample status; rank-one permission reduction; earlier modular nondegeneracy device and large-t spectral principle.

**Candidate new contributions:** all-parameter five-neuron exclusion and sharp six-neuron minimum; the explicit twin-cluster family and all-r rational parameter theorem; the complete seed stability region and whole-cone corollary; a compact rational witness and a computed external antichain.

**Standard tools, not novelty claims:** determinant multilinearity, Cramer's rule, rank-one identities, Routh–Hurwitz, Descartes' rule of signs, parity, and invariant-subspace factorization. The quintic square identity is used as a reusable elementary certificate; no unsupported claim that the algebraic identity itself is unprecedented is made.

## Explicit exclusions from the claim set

There is no claimed general graph-level classification of all CTLNs; no solution for every size-seven-and-up graph; no biological experimental validation; no assertion of global dynamical convergence or nonlinear Hopf bifurcation; no proof-assistant formalization; no external independent referee report; no guaranteed priority; no self-assigned REF rating; and no publication, authorship, or endorsement claim on behalf of the commissioning user or Evidence Press.
~~~~

#### STATUS.md
~~~~
# Status and scope

Version 1.1, 1 October 2026. Unrefereed mathematical candidate.

The sharp minimum is six for stable nonclique supports in nondegenerate legal
CTLNs with uniform positive drive. The six-neuron graph is Geneson's prior
construction; this package supplies the exact lower bound through five and the
stated rational, all-size and phase-region certificates.

The infinite family has sizes 5r+6 with r>=0 and size-dependent parameters.
It is not a construction at every size or at one common parameter pair. Its
transverse decay rate tends to zero. No biological validation, global convergence,
formal proof-assistant verification or external human peer review is claimed.

See ASSURANCE.md, REPLAY_RECEIPT.md, SOURCES.md and REVISION_RESPONSE.md.
~~~~

### 2(b) The seed graph

**As given by the repo** (three places agree; I asserted equality of the paper TeX matrix and the certificate JSON matrix):
paper eq. (15) `N_0` in `paper/main.tex`, `certificates/seed_phase_and_witness.json["nonedge_matrix"]`, and
`code/verify_phase_and_witness.py` variable `N`. It is a NONEDGE matrix (not an adjacency matrix):

```
N0 = [[0,0,0,0,1,1],
      [0,0,0,0,0,1],
      [0,0,0,0,0,1],
      [1,1,0,0,1,0],
      [0,0,1,1,0,0],
      [0,0,0,1,1,0]]
```

**Stated convention (verbatim):** README.md "Matrices encode **missing** directed edges, not present ones: `N[i][j] = 1` means
source `j` does not point to target `i`." Paper: "`N_{ii}=0` and `N_{ij}=1` precisely when the source `j` does not point to the
target `i`. Thus `sigma` is a bidirected clique exactly when `N=0`." Weights: `W_ij = 0 (i=j); -1+eps (j->i); -1-delta (j-/->i)`;
legal range `delta>0, 0<eps<delta/(1+delta)` (identical to ours). Parameters: `beta = 1+delta/eps`, `t = (1-eps)/eps`,
`K = eps^{-1}(I-W)_sigma = I + beta N + t J`; Jacobian `-eps K`. Node order in the paper's family description is `(a, A, B, c, d, e)`.
Provenance they claim: Geneson arXiv:2607.21396 eq. (19) / Theorem 3.2. MY CHECK against the Geneson PDF copy bundled in Repo 2
(`literature/2607.21396v1.pdf`, sha256 1bb6d9a3...466fe3, provenance of that copy unverified): its eq. (19) is the same 6x6 matrix, entry for entry.

**Conversion to OUR convention** (`A[i][j]=1 iff i->j`): their `N[i][j]=1 iff j-/->i`, so for `i != j`, `A[j][i] = 1 - N[i][j]`,
i.e. `A = J - I - N^T`. Result (nodes 1..6 in their order):

```
A = [[0,1,1,0,1,1],
     [1,0,1,0,1,1],
     [1,1,0,1,0,1],
     [1,1,1,0,0,0],
     [0,1,1,0,0,0],
     [0,0,0,1,1,0]]
```
Edge list (1-indexed, i->j): (1,2) (1,3) (1,5) (1,6) (2,1) (2,3) (2,5) (2,6) (3,1) (3,2) (3,4) (3,6) (4,1) (4,2) (4,3) (5,2) (5,3) (6,4) (6,5).
19 directed edges: 5 reciprocal pairs, 9 one-way pairs, 1 non-adjacent pair; out-degrees [4,4,4,3,2,2], in-degrees [3,4,4,2,3,3].

**Isomorphism to the matrix you gave me:** MY CHECK (brute force over all 720 permutations `p`, test `A[p[i]][p[j]] == GEN[i][j]`):
**ISOMORPHIC, and in fact IDENTICAL: the only isomorphism found is the identity permutation** (automorphism group of the seed is
trivial, size 1). So the converted seed equals your `A = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]`
with the same node labelling. Convention-error controls (all FALSE, as they should be): converted seed vs `GEN^T`; raw `N0` read as
adjacency vs `GEN`; raw `N0^T` vs `GEN`; `J-I-N0` (complement without transpose) vs `GEN`. So a transpose slip would have been detected.

### 2(c) Seed phase region (their C6/C7): exact polynomials, verbatim from `paper/main.tex` (LaTeX source)

**Definitions used (paper s2 and s5):**
`beta = 1 + delta/eps`, `t = (1-eps)/eps`, `u = beta - 2 = delta/eps - 1` (so `delta = (1+u) eps`), `v = u t - 1`; legal iff `u>0, v>0` iff `t u > 1`.
`det(zI_6 + I_6 + beta N_0 + t J_6) = (z+1) f_{beta,t}(z)`, `f = z^5 + a_1 z^4 + ... + a_5`, and `Delta_4(u,t)` is the 4th Hurwitz
determinant of `f` (Hurwitz convention `Delta_j = det(a_{2c-r+1})_{r,c=0}^{j-1}`; the 5th is `a_5 Delta_4`), written as a polynomial in `t`
after `beta = 2+u`:

```
\Delta_4(u,t)=\sum_{j=0}^4 B_j(u)t^j,      B_4(u)=6F(u),
F(u)=797u^6+3096u^5+5684u^4+4712u^3+1536u^2-1184u+128.          (the "u_star sextic")
```
```
B_3 = 1536-21056u+33344u^2+123312u^3+162280u^4+107488u^5+37866u^6+5158u^7
B_2 = -256-25728u+49728u^2+207584u^3+299504u^4+231240u^5+102944u^6+24250u^7+2301u^8
B_1 = -2048-16896u+32000u^2+156160u^3+249728u^4+220528u^5+118400u^6+38136u^7+6722u^8+496u^9
B_0 = -1024-5120u+6400u^2+42752u^3+77440u^4+77856u^5+49136u^6+19872u^7+4988u^8+706u^9+43u^{10}
```
(`T(u)` = "the unique positive root of eq. (phaseH4)", i.e. the unique positive zero in `t` of the quartic `sum_j B_j(u) t^j`.)
Characteristic coefficients (eq. seedcoefficients), `beta`-form:
```
a_1 = 5+6t
a_2 = 10-beta^2+(24-11 beta)t
a_3 = 10-3beta^2+4beta^3+(36-33beta+16beta^2)t
a_4 = 5-3beta^2+8beta^3-3beta^4+(24-33beta+32beta^2-10beta^3)t
a_5 = 1-beta^2+4beta^3-3beta^4+beta^5+(6-11beta+16beta^2-10beta^3+3beta^4)t
```
Other quoted polynomials: `d(beta)=det(I+beta N_0)=beta^5-3beta^4+4beta^3-beta^2+1`; permission numerators
`(beta-1)(beta^3-2beta^2+beta-1), beta^2-beta+1, beta^2-beta+1, (beta-1)^2(beta^2-beta+1), -beta^3+3beta^2-2beta+1, beta^4-3beta^3+4beta^2-2beta+1`;
`Delta_2 = (12-66u)t^2+(12-110u-22u^2)t-(56u+26u^2+4u^3)`; `G(u)=32-162u-141u^2-116u^3` (leading coefficient of `Delta_3` in `t` is `6G`);
whole-cone polynomial `u^4 Delta_4(u,1/u) = 768-5568u-12096u^2+33840u^3+189224u^4+415320u^5+574334u^6+561586u^7+406070u^8+220506u^9+89573u^10+26594u^11+5484u^12+706u^13+43u^14`.

**Theorem 5.1 (Sharp seed phase region), verbatim:**
> The polynomial F has a unique zero u_* in (0,2/11), with 4/25 < u_* < 1/6, 0.166355027849 < u_* < 0.166355027850.
> For each 0<u<u_*, the quartic Delta_4(u,t) has exactly one positive zero T(u). At legal parameters, Geneson's graph (Nseed) has a
> stable positive full-support fixed point if and only if  0 < u < u_*,  t > max{1/u, T(u)}.
> Equivalently, with delta=(1+u)eps,  0<u<u_*,  0<eps<min{ u/(1+u), 1/(1+T(u)) }.

**Corollary 5.2 (whole legal cone), verbatim:** "For Geneson's graph, every legal pair satisfying 1 < delta/eps <= 11/10 has a stable positive
full-support fixed point." Also (paper, after Cor 5.2): "T(1/10) ~ 4.28377051836, T(4/25) ~ 52.93123957158 ... the previously used ray
delta/eps = 29/25 has its full-support threshold explicitly resolved by u = 4/25."

**Your sextic vs their `F`.** MY CHECK (sympy, exact): substituting `u = b - 2` into their `F` gives exactly
`797 b^6 - 6468 b^5 + 22544 b^4 - 44440 b^3 + 53280 b^2 - 38016 b + 13824`, i.e. **identical to the sextic you quote (difference = 0)**,
with `b = 1 + delta/eps = beta`. Real roots of that sextic: b = 2.166355027849143532380273 (= 2 + u_*, the smallest root > 2) and
b = 2.197503670655442785160828 (u = 0.19750...). The second real root is NOT a stability point: the paper excludes `u >= 2/11 = 0.1818` by
`Delta_2 < 0` there (Delta_2 has nonpositive coefficients and negative constant term for u >= 2/11) and `F>0` again beyond 0.19750 is irrelevant.
Positivity: their fifth permission numerator is `-beta^3+3beta^2-2beta+1 = 1-2u-3u^2-u^3` (they only need it >0 on u<2/11; value 707/1331 at u=2/11).
MY CHECK: `b^3-3b^2+2b-1` has the single real root 2.32471795724474602596 and equals `r^3-r-1` under `b = 1+r`, so the root is `1 + rho`
(rho = plastic number 1.3247...): **your positivity statement matches**. Since 2 + u_* = 2.1664 < 2.3247 the positivity bound is not binding.
So in the eps -> 0 limit (t -> infinity, any fixed u) their region reduces to `0 < u < u_*`, i.e. `2 < b < 2.16636`, **exactly your stated result**.

Independent derivation (MY CHECK): I rebuilt `W` from OUR rule applied to the converted `A` (not from their formulas), formed `K = eps^{-1}(I-W)` with
`eps=1/(1+t), delta=(beta-1)eps`, confirmed `K = I + beta N_0 + t J` symbolically, divided `det(zI+K)` by `(z+1)`, computed the 4th Hurwitz determinant,
and substituted `beta = 2+u`. **`B_0, B_1, B_2, B_3` and `B_4 = 6F` all equal their listed coefficient vectors exactly**, and `a_1..a_5` equal their printed formulas.
`T(1/10)=4.28377051835558`, `T(4/25)=52.9312395715791` (my sympy real roots; theirs 4.28377051836 and 52.93123957158). The Cramer numerators equal theirs.
Pointwise exact tests (OUR pipeline, Fractions): u=1/10 stable at t = 10+1e-12 (q-1 = 9e-15), t=11, t=1e6; u=4/25: t=49 UNSTABLE (<T=52.93), t=53 and t=54 STABLE;
u = 0.17, 1/6, 0.2 UNSTABLE even at t = 1e6..1e9. Grid: u in {k/200 : k=1..80}, ~19 t-values each (just-above-1/u, multiples of 1/u, 3..1e9), legal points only, 1,353 points
compared against "claimed region" (u>=u_* upper isolator => unstable; else stable iff Delta_4(u,t)>0): **0 mismatches, 454 stable points, 0 points with a vanishing Hurwitz minor.**
Near q=1 (exact): stable with x>0 for u in {1e-1, 1e-2, 1e-3, 1e-5}, t = 1/u + 1e-12 .. 1e12 (q-1 from 1e-22 to 0.1).
NOT verified by me: the proof of Thm 5.1 itself (Descartes arguments, ordering t2<t3<T, coefficient bounds), only its polynomial inputs and a sample of its conclusion.

### 2(d) Infinite family on 5r+6 vertices

Verbatim from `paper/main.tex` (section "Graph construction and explicit parameters"):
~~~~
\subsection{Graph construction and explicit parameters}
For each integer $r\geq0$, introduce vertices $a,c,d,e$ and a twin pool
$T=A\sqcup B$ with
\[
 |A|=3r+1,\qquad |B|=2r+1,\qquad |T|=5r+2.
\]
There are $n=5r+6$ vertices. Specify the \emph{missing incoming sources} as follows; every other off-diagonal directed edge is present:
\begin{center}
\begin{tabular}{ll}\toprule
Target & Missing incoming sources\\\midrule
$a$ & $\{d,e\}$\\
each vertex of $T$ & $\{e\}$\\
$c$ & $\{a,d\}\cup A$\\
$d$ & $\{c\}\cup B$\\
$e$ & $\{c,d\}$\\\bottomrule
\end{tabular}
\end{center}
Let this graph be $G_r$ and its nonedge matrix be $N_r$. At $r=0$ it is \eqref{eq:Nseed}, under the ordering $(a,A,B,c,d,e)$. The graph contains the clique $T\cup\{a\}$ of size $n-3$, but is itself not a clique.

For the classes $(\{a\},T,\{c\},\{d\},\{e\})$, the quotient data are
\begin{equation}\label{eq:familyquotient}
 M_r=\begin{pmatrix}
0&0&0&1&1\\
0&0&0&0&1\\
1&3r+1&0&1&0\\
0&2r+1&1&0&0\\
0&0&1&1&0
\end{pmatrix},\qquad
 w_r=(1,5r+2,1,1,1)^T.
\end{equation}

\begin{theorem}[Unbounded stable nonclique supports]\label{thm:family}
For every integer $r\geq0$, set
\begin{equation}\label{eq:familyparams}
 D=30r+31,\quad
 \eps=\frac1{2D+1},\quad
 \delta=\frac{D+2}{D(2D+1)}.
\end{equation}
For every drive $\theta>0$, the CTLN on $G_r$ at these parameters is legal and globally nondegenerate, and has a stable fixed point with nonclique full support of size $5r+6$.
\end{theorem}
~~~~
**Parameters (Thm 4.3, verbatim):** for every integer r >= 0, `D = 30r+31`, `eps = 1/(2D+1)`, `delta = (D+2)/(D(2D+1))`; "For every drive theta>0, the CTLN on G_r at
these parameters is legal and globally nondegenerate, and has a stable fixed point with nonclique full support of size 5r+6." Then `beta = 2+2/D`, `t = 2D`, legal margin `t(beta-2)=4`.
In our q: `q = delta(1-eps)/eps = 2(D+2)/(2D+1) = 1 + 3/(2D+1)` (r=0: 22/21; r=1: 42/41; r=2: 62/61; r=3: 82/81), i.e. q -> 1+ as r grows.
Caveats stated by the repo itself (STATUS.md; paper Remark "Size-dependent stability"): parameters depend on r; the 5r+1 transverse eigenvalues of the physical Jacobian are
`-eps_r = -1/(60r+63)` -> 0 (no uniform decay margin); NOT a construction at every size and NOT at one common (eps,delta). Sizes are 6, 11, 16, 21, ...; **nothing at n = 7.**
The graph is the clique expansion of the seed with multiplicities (1, 3r+1, 2r+1, 1, 1, 1) for (a, A, B, c, d, e); the paper credits composite-graph theory (CGM 2024 s2.2, Lemma 2.6, s5) for the reduction and claims only the multiplicities/schedule/certificates.

**Missing sources, in plain terms** (target <- sources that do NOT send to it; every other ordered pair is an edge):
a <- {d,e}; every vertex of A and of B <- {e}; c <- {a, d} plus all of A; d <- {c} plus all of B; e <- {c, d}.

**The r = 1 graph (n = 11) in OUR convention** (`A[i][j]=1 iff i->j`). Node order (0-indexed): 0=a, 1..4=A0..A3, 5..7=B0..B2, 8=c, 9=d, 10=e.
(|A| = 3r+1 = 4, |B| = 2r+1 = 3.) It was built by me from the missing-source table above; I checked that r=0 reproduces the 6-node seed matrix exactly and that the class
counts match their quotient matrix `M_r`.

Adjacency matrix A (rows = sources):
   [0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1]
   [1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1]
   [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1]
   [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
   [1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1]
   [1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]
   [1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1]
   [1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1]
   [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
   [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
   [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0]

Edge list i->j (89 edges of 110 ordered pairs), 0-indexed:
   [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 9), (0, 10), (1, 0), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 9), (1, 10), (2, 0), (2, 1), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 9), (2, 10), (3, 0), (3, 1), (3, 2), (3, 4), (3, 5), (3, 6), (3, 7), (3, 9), (3, 10), (4, 0), (4, 1), (4, 2), (4, 3), (4, 5), (4, 6), (4, 7), (4, 9), (4, 10), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 6), (5, 7), (5, 8), (5, 10), (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 7), (6, 8), (6, 10), (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 8), (7, 10), (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6), (9, 7), (10, 8), (10, 9)]

The 21 ABSENT ordered pairs (src -/-> tgt):
   [('a', 'c'), ('A0', 'c'), ('A1', 'c'), ('A2', 'c'), ('A3', 'c'), ('B0', 'd'), ('B1', 'd'), ('B2', 'd'), ('c', 'd'), ('c', 'e'), ('d', 'a'), ('d', 'c'), ('d', 'e'), ('e', 'a'), ('e', 'A0'), ('e', 'A1'), ('e', 'A2'), ('e', 'A3'), ('e', 'B0'), ('e', 'B1'), ('e', 'B2')]

**MY CHECK of the family claim in our convention** (exact Fractions, W from A per CLAUDE.md s2, full support so no outside condition):

| r | n | D | eps | delta | q | x*>0 | all Hurwitz minors >0 |
|---|---|---|---|---|---|---|---|
| 0 | 6 | 31 | 1/63 | 11/651 | 22/21 | yes | yes (none zero) |
| 1 | 11 | 61 | 1/123 | 21/2501 | 42/41 | yes | yes (none zero) |
| 2 | 16 | 91 | 1/183 | 31/5551 | 62/61 | yes | yes (none zero) |
| 3 | 21 | 121 | 1/243 | 41/9801 | 82/81 | yes | yes (none zero) |

Global nondegeneracy (all principal minors and all column-replacement determinants of I-W nonzero), exact: r=0 (63 + 192 determinants) and r=1 (2047 + 11264 determinants): **zero vanishing**.
Not checked: global nondegeneracy for r>=2 (their parity-mod-2 argument is stated for all r; I did not re-derive it); the all-r symbolic Hurwitz/permission certificates (`certificates/infinite_family.json`, `verify_family.py`).

### 2(e) The n <= 5 exclusion: what is certified, how, and whether the certificates are machine-checkable without their code

Verbatim from `paper/main.tex`, section 3 (LaTeX source; Theorem 3.1, certificate method, certificate-type table, final exceptional cases, Lemma 3.2, proof):
~~~~
\section{The five-neuron exclusion and the sharp minimum}\label{sec:minimum}
\begin{theorem}[Exact small-support classification]\label{thm:small}
For every nonzero zero-one matrix $N$ with zero diagonal and order $m\leq5$, and every legal pair $(\beta,t)$, the matrix
\[
 I_m+\beta N+tJ_m
\]
cannot be both permitted and positive stable. Consequently, the supports of size at most five satisfying the algebraic stable-support conditions are exactly the target-free cliques. The same conclusion holds for stable fixed points of any nondegenerate legal CTLN, regardless of the ambient network size.
\end{theorem}

\subsection{Finite certificate method}
The finite part of the proof checks one canonical representative of each isomorphism class, while verifying coverage by explicitly marking its whole vertex-permutation orbit. Off-diagonal entries are encoded in row-major order, with the first entry in bit zero. Thus an integer mask determines $N$, with no external graph catalogue.

For each representative, the checker reconstructs $\det C$, all Cramer numerators, and \eqref{eq:characteristic} using direct signed-permutation determinant expansions. The separate certificate generator uses Newton identities. The certificate records only the graph mask and an obstruction instruction; it does not provide trusted polynomial coefficients.

The allowed obstructions are the following.
\begin{enumerate}[label=(\roman*)]
\item A determinant or Cramer numerator has all coefficients nonpositive after $\beta=2+u$. This contradicts \eqref{eq:cramernecessary} for every $u>0$.
\item The characteristic polynomial vanishes identically at $z=\beta-1$, giving the negative eigenvalue $1-\beta$ of $K$.
\item A characteristic coefficient, or a Hurwitz determinant, becomes coefficientwise nonpositive after substituting \eqref{eq:uv} and multiplying by a sufficient positive power of $u$.
\item An exact factorization supplies an unstable quadratic or a quotient with a nonpositive coefficient.
\item The quintic split argument proved below supplies a sign obstruction on the two regions $0<u<2$ and $u\geq2$.
\end{enumerate}
These are pointwise algebraic implications on the entire parameter region, not numerical grid tests. Empty polynomials represent the identically zero obstruction, which also contradicts required strict positivity.

\begin{table}[ht]
\centering
\caption{Exhaustive graph coverage verified by explicit permutation orbits.}
\begin{tabular}{rrr}\toprule
Order & Isomorphism classes & Labelled matrices\\\midrule
1&1&1\\2&3&4\\3&16&64\\4&218&4,096\\5&9,608&1,048,576\\\midrule
Total&9,846&1,052,741\\\bottomrule
\end{tabular}
\end{table}

\begin{table}[ht]
\centering
\caption{Disjoint certificate types at order five. The clique is the sole permitted exception to exclusion.}
\begin{tabular}{lr}\toprule
Certificate type & Number of classes\\\midrule
Clique ($N=0$)&1\\
Necessary determinant/Cramer sign&9,203\\
Eigenvalue $1-\beta$&242\\
Characteristic-coefficient sign&134\\
Hurwitz-determinant sign&22\\
Unstable directed-cycle quadratic&2\\
Quadratic factor followed by quotient sign&1\\
Quintic discriminant split&3\\\midrule
Total&9,608\\\bottomrule
\end{tabular}
\end{table}

\subsection{The final exceptional cases}
After the determinant, coefficient, eigenvalue, and direct Hurwitz screens, six five-vertex masks remain:
\[
74883,\quad76900,\quad76966,\quad79298,\quad79491,\quad109731.
\]
For masks 79491 and 109731 the characteristic polynomial has the exact factor
\begin{equation}\label{eq:cyclefactor}
 (z+1)^2-\beta(z+1)+\beta^2.
\end{equation}
Its linear coefficient $2-\beta$ is negative, so it is not Hurwitz stable.
For mask 76966 a factor is $(z+1)^2+\beta^2$ and the quotient is
\begin{align}\label{eq:quotient}
 &(z+1)^3-\beta^2(z+1)+2\beta^3\notag\\
 &\hspace{13mm}+t\bigl(5(z+1)^2-8\beta(z+1)+7\beta^2\bigr).
\end{align}
Its coefficient of $z$ is $3-\beta^2+(10-8\beta)t<0$ for $\beta>2$, $t>0$.

The remaining three cases use a reusable identity.
\begin{lemma}[Quintic discriminant obstruction]\label{lem:quintic}
For a monic real quintic
$p=z^5+a_1z^4+a_2z^3+a_3z^2+a_4z+a_5$, put
$X=a_1a_2-a_3$ and $Y=a_1a_4-a_5$. Then
\begin{align}
 \Delta_3&=a_3X-a_1Y,\label{eq:d3}\\
 \Delta_4&=(a_3a_4-a_2a_5)X-Y^2,\label{eq:d4}\\
 4a_1^2\Delta_4&=(a_3^2-4a_1a_5)X^2-(2a_1Y-a_3X)^2.\label{eq:square}
\end{align}
In particular, if $a_1>0$ and $a_3^2-4a_1a_5\leq0$, then $p$ is not Hurwitz stable.
\end{lemma}
\begin{proof}
Expand the third and fourth Hurwitz determinants to obtain \eqref{eq:d3}--\eqref{eq:d4}. Substituting $Y=a_1a_4-a_5$ and completing the square gives \eqref{eq:square}. The last assertion follows from $\Delta_4\leq0$.
\end{proof}

For each mask $h\in\{74883,76900,79298\}$ let $p_h$ be its characteristic quintic, with $a_1=5(1+t)>0$. Set
\[
 P_h(u,v)=-u^2(a_3^2-4a_1a_5)\bigm|_{\beta=2+u,\,t=(1+v)/u}.
\]
In each case every nonzero coefficient of $P_h$ is strictly positive except the coefficient $-1$ of $u^8$. The coefficient of $u^7$ is respectively $8,14,20$. Hence $P_h>0$ for $0<u\leq2$, $v>0$, because $c u^7-u^8=u^7(c-u)>0$ and all other nonzero terms are positive. On $u\geq2$, writing $u=2+x$, the polynomial
\[
 -u^2\Delta_2(p_h)\bigm|_{\beta=2+u,\,t=(1+v)/u,\,u=2+x}
\]
has every nonzero coefficient strictly positive, including its constant term. Thus $\Delta_2<0$ there. The explicit coefficient vectors appear in Appendix~\ref{app:hard}. The checker also validates the bounded interval through $u=2x/(1+x)$, an algebraically independent sign presentation of the same inequality.

\begin{proof}[Proof of Theorem~\ref{thm:small}]
The certificate file and checker establish the complete orbit cover in Table~1, with one valid obstruction for every nonzero representative. Each obstruction is sufficient by Lemma~\ref{lem:permission}, the Hurwitz criterion, the factor tests, or Lemma~\ref{lem:quintic}. This proves the matrix assertion. If $N=0$, then $K=I+tJ$ is permitted and positive stable. The external-input test in Section~\ref{sec:external} reduces exactly to target-freeness. The dynamical assertion follows from the nondegenerate support criterion of Section~\ref{sec:model}.
\end{proof}
~~~~

Their summary of the checker (README.md "Reproduce every mathematical check", verbatim excerpt): "The only finite graph certificate consumed as instructions is `certificates/small_order_obstructions.json`. Polynomial coefficients are reconstructed, not trusted from that file. The family and phase scripts recompute their polynomial certificates from the displayed matrices and verify sign conditions; their JSON files are inspectable outputs, not an external oracle."

**What is certified (their claim).** For every nonzero zero-diagonal 0/1 nonedge matrix N of order m <= 5 (9,846 isomorphism classes; 1,052,741 labelled graphs; per order 1,3,16,218,9608
classes, which match A000273; Repo 2's ledger attributes the 9,608 figure to Curto-Langdon-Morrison 2019 Table 1 p. 6, which I have not checked against Curto's n=5 package), the matrix `K = I + beta N + t J` is never both permitted (`K^{-1}1 > 0`) and positive stable for ANY legal `(beta,t)`, i.e. `beta = 2+u, t = (1+v)/u, u>0, v>0`.
Because stability conditions 1 and 3 depend only on the induced subgraph on the support, this is a statement about the full support of every induced subgraph of order <= 5; combined with the clique-target computation it gives "stable supports of size <= 5 are exactly target-free cliques" (nondegenerate setting). Condition 2 (outside nodes) is not needed for nonclique supports.
**How:** one canonical mask per isomorphism class (row-major off-diagonal nonedge bits, bit 0 first, `N_ij = 1` means j does not point to i); a pre-computed "obstruction instruction" per class of one of the kinds
`clique, necessary_sign, eigenvalue_1_minus_beta, coefficient, hurwitz, cycle_quadratic, quadratic_quotient, quintic_split`; the checker rebuilds `det C`, the Cramer numerators and the characteristic polynomial by direct signed-permutation (Leibniz) expansion in integer polynomials in (beta,t), applies the substitution `beta=2+u` and checks the sign pattern of the coefficient vector, and marks the whole `S_n`-orbit to prove complete labelled coverage.
Counts by type at n=5 (from the JSON): clique 1, necessary_sign 9203, eigenvalue_1_minus_beta 242, coefficient 134, hurwitz 22, cycle_quadratic 2, quadratic_quotient 1, quintic_split 3 (total 9608).

**Machine-checkable without their code? No, not as shipped.** `certificates/small_order_obstructions.json` (1,004,822 bytes) holds only records like
`{"n":5,"mask":74883,"kind":"quintic_split","unit_factors":0,"split_u":2}` or `{"n":3,"mask":25,"kind":"coefficient","unit_factors":0,"coefficient":2}` or `{"n":2,"mask":1,"kind":"necessary_sign","polynomial":1}`.
There are NO polynomials, coefficient vectors, or identities in it; a checker must itself reconstruct `det C`, the Cramer numerators and `det(zI+K)` for each mask and then test the stated sign pattern. The paper's printed appendix gives coefficient vectors only for the 3 hardest masks (74883, 76900, 79298). So an independent party can re-check it only by writing a new implementation from the paper's description (the format is simple and documented: mask encoding + kind + index), which I did NOT do for the all-parameter obstructions. The repo itself concedes the shared-primitive risk (generator and verifier share `exact_polys.py`; "not a formal proof"). Their own logs (`verification/*.log`) say all five scripts passed on Python 3.13.5; I did not run them (instruction) and did not see the GitHub Actions run.

**MY independent sanity checks of the n <= 5 claim (partial; these are NOT a re-proof):**
1. Using only the list of masks from their JSON as data and MY OWN implementation of the vertex-permutation action on their documented bit encoding: for n = 1..5 the listed representatives (1, 3, 16, 218, 9608) have pairwise-disjoint orbits that cover all 1, 4, 64, 4096, 1,048,576 labelled graphs. PASS. (So their class counts and "complete orbit cover" are confirmed independently.)
2. Pointwise exact grid in OUR convention (Fractions; `A[j][i] = 1 - N[i][j]`; `eps = 1/(1+t)`, `delta = (1+u) eps`): u in {1/20, 1/10, 1/5, 1/2, 1, 3, 10} times t*u in {101/100, 2, 10, 100, 10^4}, i.e. 35 legal points, over every nonclique canonical graph of order 2..5 (2 + 15 + 217 + 9,607 = 9,841 graphs). **No graph at any grid point has both `x* = (I-W)^{-1} 1 > 0` and all Hurwitz minors > 0 (0 violations).** For n=5: ~336k evaluations, 23,945 permitted-only points, 1,457 stable-only points, 0 both; 163 points with a vanishing Hurwitz minor and 6 with singular `I-W` (undecided, counted separately and not classified either way, per the three-valued rule; n=4: 3 vanishing + 9 singular). This is consistent with their theorem but is a 35-point sample per graph and says nothing about the thin regions between sample points that the algebraic proof is meant to cover.
### 2(f) The paper PDF

`paper/six_is_minimal_ctln_v1_1.pdf` (sha256 174c14ce61ab9f0604864f2d315a3051cddd88bbdd817eb0785a4f3585dc89a3; 15 pages; pdfTeX 1.40.29; created 2026-10-01 19:30 UTC; PDF metadata Author: "Anonymous"; Subject: "Evidence Press research candidate; computer-assisted mathematical proof"). The PDF text for Thm 3.1, Thm 4.3, Thm 5.1, Cor 5.2 agrees with `main.tex` (I compared the extracted text).

**Title:** Six Is Minimal: Exact Certification and Infinite Families of Stable Nonclique CTLNs
**Author line:** Anonymous. **Date line:** 1 October 2026 | Version 1.1 (candidate). Running footer: "Prepared for Evidence Press". Keywords: Threshold-linear network; stable fixed point; target-free clique; exact enumeration; Routh-Hurwitz criterion; equitable reduction. MSC 2020: 05C20, 37N25, 92B20.

**Abstract (pdftotext, verbatim):**
~~~~

                                   Abstract
We resolve the five-neuron gap left by the six-neuron counterexample to the target-free
clique conjecture for combinatorial threshold-linear networks (CTLNs). For every legal
parameter pair, a support of size at most five satisfies the algebraic stable-support conditions
if and only if it is a target-free clique. Consequently, among nondegenerate legal CTLNs, the
minimum size of a stable nonclique support is exactly six. The lower bound is a computer-
assisted proof covering all 1,052,741 labelled loopless digraphs of orders one through five,
with independently reconstructed integer-polynomial obstructions on the entire unbounded
parameter region. The final exceptional cases are excluded by a quintic discriminant identity,
rather than numerical sampling.
Beyond minimality, we give an exact parameter classification for Geneson’s six-neuron
graph and a twin-cluster construction with stable nonclique full supports of every size
5r + 6, r ≥ 0. The latter has explicit rational legal parameters, an all-size polynomial
stability certificate, and a parity proof of global nondegeneracy. We also state the weighted
missing-input rule for implanting these supports in larger networks. The accompanying
standard-library Python checker uses only integers and exact rational fractions. This is a
candidate mathematical result, not a claim of completed external peer review or proof-assistant
~~~~

**Theorem statements, verbatim from `main.tex` (LaTeX source).**

Lemma 2.1 (Rank-one permission reduction):
~~~~
\begin{lemma}[Rank-one permission reduction]\label{lem:permission}
Let $C=I+\beta N$, where $N\geq0$ and $\beta,t\geq0$. Then $C+tJ$ is permitted if and only if $C$ is permitted. When this holds, with $y=C^{-1}\ones$ and $s=\ones^Ty$,
\begin{equation}\label{eq:permission}
 (C+tJ)^{-1}\ones=\frac{y}{1+ts},\qquad
 \det(C+tJ)=\det C\,(1+ts).
~~~~
Theorem 3.1 (Exact small-support classification):
~~~~
\begin{theorem}[Exact small-support classification]\label{thm:small}
For every nonzero zero-one matrix $N$ with zero diagonal and order $m\leq5$, and every legal pair $(\beta,t)$, the matrix
\[
 I_m+\beta N+tJ_m
\]
cannot be both permitted and positive stable. Consequently, the supports of size at most five satisfying the algebraic stable-support conditions are exactly the target-free cliques. The same conclusion holds for stable fixed points of any nondegenerate legal CTLN, regardless of the ambient network size.
\end{theorem}
~~~~
Lemma 3.2 (Quintic discriminant obstruction):
~~~~
\begin{lemma}[Quintic discriminant obstruction]\label{lem:quintic}
For a monic real quintic
$p=z^5+a_1z^4+a_2z^3+a_3z^2+a_4z+a_5$, put
$X=a_1a_2-a_3$ and $Y=a_1a_4-a_5$. Then
\begin{align}
 \Delta_3&=a_3X-a_1Y,\label{eq:d3}\\
 \Delta_4&=(a_3a_4-a_2a_5)X-Y^2,\label{eq:d4}\\
 4a_1^2\Delta_4&=(a_3^2-4a_1a_5)X^2-(2a_1Y-a_3X)^2.\label{eq:square}
\end{align}
In particular, if $a_1>0$ and $a_3^2-4a_1a_5\leq0$, then $p$ is not Hurwitz stable.
\end{lemma}
~~~~
Corollary 3.3 (Sharp minimum):
~~~~
\begin{corollary}[Sharp minimum]\label{cor:minimum}
Among nondegenerate legal CTLNs, the minimum cardinality of a stable nonclique fixed-point support is exactly six.
\end{corollary}
~~~~
Lemma 4.1 (Identical-row cluster reduction):
~~~~
\begin{lemma}[Identical-row cluster reduction]\label{lem:twin}
Partition the vertices into $k$ classes such that rows of $N$ agree within each class. Let $m_b$ be the size of class $b$, and let $M_{ab}$ be the number of nonedges in a row of class $a$ whose source is in class $b$. Write $w=(m_1,\ldots,m_k)^T$. Then
\begin{equation}\label{eq:quotientfactor}
 \det(zI_n+I_n+\beta N+tJ_n)
 =(z+1)^{n-k}\det(zI_k+Q),
 \qquad Q=I_k+\beta M+t\ones_k w^T.
\end{equation}
Permission is equivalent to permission of $I_k+\beta M$, with its solution replicated inside each class. Positive stability is equivalent to positive stability of $Q$.
\end{lemma}
~~~~
Theorem 4.3 (Unbounded stable nonclique supports):
~~~~
\begin{theorem}[Unbounded stable nonclique supports]\label{thm:family}
For every integer $r\geq0$, set
\begin{equation}\label{eq:familyparams}
 D=30r+31,\quad
 \eps=\frac1{2D+1},\quad
 \delta=\frac{D+2}{D(2D+1)}.
\end{equation}
For every drive $\theta>0$, the CTLN on $G_r$ at these parameters is legal and globally nondegenerate, and has a stable fixed point with nonclique full support of size $5r+6$.
\end{theorem}
~~~~
Theorem 5.1 (Sharp seed phase region): see 2(c). Corollary 5.2 (whole legal cone):
~~~~
\begin{corollary}[An explicit whole legal cone]\label{cor:cone}
For Geneson's graph, every legal pair satisfying
\[
 1<\delta/\eps\leq11/10
~~~~
Proposition 6.1 (Weighted missing-source rule):
~~~~
\begin{proposition}[Weighted missing-source rule]\label{prop:external}
Let an active motif satisfy $y=(I+\beta N)^{-1}\ones>0$ and let $s=\ones^Ty$. For an outside neuron $k$, let $M_k\subseteq\sigma$ be the active sources with no edge to $k$. Then
\begin{equation}\label{eq:external}
 x^*_\sigma=\frac{\theta y}{\eps(1+ts)},\qquad
 \theta+\sum_{j\in\sigma}W_{kj}x_j^*
 =\frac{\theta\bigl(1-\beta\sum_{j\in M_k}y_j\bigr)}{1+ts}.
\end{equation}
Consequently $k$ is strictly inactive exactly when
\begin{equation}\label{eq:weighted}
 \beta\sum_{j\in M_k}y_j>1.
\end{equation}
If the motif is positive stable and every added neuron satisfies this strict inequality, the embedded equilibrium is locally exponentially stable, whatever the edges back to the support or among the inactive neurons.
\end{proposition}
~~~~
The proposition is a rephrasing of the standard off-neuron test (they say so: CLAIMS.md C8). Its formula for the outside input,
`theta + sum_j W_kj x_j = theta(1 - beta*sum_{j in M_k} y_j)/(1 + t s)` with `y = (I+beta N)^{-1} 1`, `s = 1^T y`, is relevant to our Stage 4 (it reduces condition 2 for any permitted nonclique motif to a weighted missing-source inequality `beta*sum_{j in M_k} y_j > 1`); I did not re-derive it.

**Reproducibility surface (all unexecuted by me):** `code/verify_all.py --regenerate --read-only`, stdlib-only Python >= 3.9, assertions required; `MANIFEST.sha256` (not checked by me); `.github/workflows/verify.yml`; Zenodo DOI 10.5281/zenodo.23089682 (not resolved by me); `REPLAY_RECEIPT.md` defers to an external `REPLAY_SUMMARY.json` that is NOT in the repo. Assurance level in their own words (ASSURANCE.md): "unrefereed, AI-assisted computer-assisted mathematical candidate"; "not ... proof-assistant formalisation, independent human peer review or institutional independence"; generator and verifier "share low-level routines".

## 3. Repo 2: ChickenNungets/CTLN (Aug 2026)

Repo layout: `README.md`, `AGENTS.md`, `docs/{PROJECT_SPEC,CONVENTIONS,THEOREM_LEDGER,CERTIFICATE_FORMAT}.md`, `docs/plans/`, `src/ctln_search/{graph,filters,isomorphism,model,feasibility,search,audit}.py`, `tests/` (6 files), `legacy/ctln_small_searchF.py`, `literature/` (10 PDFs, see s4), and `certificates/experiments/five_neuron_oracle_v1/`. No LICENSE file (their plan doc says "A repository license requires an owner choice and is not inferred here"). Last commit "markdown math fix" 2026-08-08.

### 3(a) README.md, in full, verbatim
~~~~
# Five-neuron CTLN exhaustive search

This repository checks the following statement:

> A globally nondegenerate five-neuron combinatorial threshold-linear network
> cannot have a stable full-support fixed point unless its graph is a
> bidirected clique.

The completed search supports the statement: all 875 graph classes that
required algebraic testing were infeasible. There were no counterexamples and
no solver timeouts or unknown results.

This is a transparent exact-oracle computation, not a proof-assistant
verification. The graph enumeration and formulas use exact integer and
rational arithmetic, but the final infeasibility decisions trust Z3's solver
for nonlinear real arithmetic.

## What the search does

### 1. Cover every graph

A loopless directed graph on five labeled vertices has 20 possible directed
edges, so the program enumerates all

$$
2^{20}=1{,}048{,}576
$$

labeled graphs. It applies the following reviewed graph-theoretic exclusions in
order. The theorem IDs and their assumptions are stated precisely in
[`docs/THEOREM_LEDGER.md`](docs/THEOREM_LEDGER.md).

| Exclusion | Labeled graphs removed | Scope |
|---|---:|---|
| the full five-clique | 1 | outside the nonclique claim |
| graphical domination | 898,560 | all legal parameters |
| uniform in-degree three | 1,024 | all legal parameters |
| maximum in-degree at most two | 55,991 | globally nondegenerate CTLNs |
| a bidirected four-clique | 295 | globally nondegenerate CTLNs |
| oriented graph with no sink | 60 | globally nondegenerate CTLNs |

These are first-match counts: each graph is recorded under the first applicable
exclusion. In the graphical-domination row, a vertex `k` dominates `j` when
`j -> k`, not `k -> j`, and every other input to `j` is also an input to `k`;
such a `j` cannot be active. The other filters are the degree and subgraph
conditions shown in the table. Together they leave 92,645 labeled graphs.

The remaining graphs are grouped up to relabeling. For each graph, the program
tries all 120 vertex permutations and chooses the smallest bit-mask
representative. This gives 875 isomorphism classes. The coverage file retains
an explicit relabeling permutation for every one of the 92,645 labeled graphs.

### 2. Replace dynamics by exact inequalities

For a graph with receiving-row nonedge matrix $N$, the program constructs

$$
B(\varepsilon,q)=\varepsilon(1-\varepsilon)I
 +\varepsilon(1-\varepsilon+q)N+(1-\varepsilon)^2J.
$$

A stable full-support fixed point exists exactly when

$$
B^{-1}\mathbf 1>0
\quad\text{and}\quad
B\text{ has eigenvalues with positive real parts}.
$$

The first condition becomes five strict polynomial inequalities by Cramer's
rule. The second becomes five strict polynomial inequalities using the
degree-five Routh--Hurwitz criterion. No eigenvalues are approximated
numerically.

A reviewed parameter-cutoff theorem handles the region
$q\geq3-\varepsilon$. Thus the algebraic search only needs the open triangle

$$
0<\varepsilon<1,
\qquad 1<q<3-\varepsilon.
$$

For each of the 875 graph classes, the program asks whether all four domain
inequalities, five Cramer inequalities, and five stability inequalities can
hold simultaneously. The exact SMT-LIB queries contain integer polynomial
constraints; this is an all-parameter satisfiability question, not a grid or
random sample.

### 3. Check and replay the result

Z3 4.13.3 returned `unsat` for all 875 systems. The audit program then:

- rebuilds the graph enumeration and theorem-filter counts;
- checks every relabeling witness and orbit size;
- verifies every saved query hash and raw solver response;
- regenerates every query from its canonical graph and requires byte-for-byte
  equality; and
- optionally reruns all 875 saved queries.

The direct algebraic systems do not impose global nondegeneracy, so each
algebraic `unsat` result is stronger than required for its graph. Some of the
graph-theoretic filters do use the stated global-nondegeneracy assumption.

## Running the search

Python 3.11 or newer and a [Z3](https://github.com/Z3Prover/z3) executable on
`PATH` are required. The `solver` extra installs the recorded Z3 4.13.3
version alongside [SymPy](https://www.sympy.org/). Exact replay should use this
version: solver upgrades can legitimately change `unsat` into a timeout-backed
`unknown` without producing a contradictory `sat` result.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[solver]'
z3 --version
```

The activation step ensures the search and audit find the pinned executable in
`.venv/bin` before any system-wide Z3 installation.

Run the complete, resumable search:

```bash
.venv/bin/ctln-search-all \
  --output-dir /tmp/ctln-exhaustive-v1 \
  --timeout-ms 30000
```

The output directory contains a manifest, graph-coverage records, one result
per graph class, every exact SMT-LIB query, and every raw solver response. An
interrupted run resumes from `results.jsonl`. Existing `sat`, `unsat`, and
`unknown` records are preserved by default. To retry only `unknown` records
with a larger timeout, use:

```bash
.venv/bin/ctln-search-all \
  --output-dir /tmp/ctln-exhaustive-v1 \
  --timeout-ms 600000 \
  --retry-unknown
```

Retries use attempt-specific artifact names and atomically replace the old
JSONL record after the new evidence has been written.

On the recorded run, constructing the 875 symbolic systems took about 26
minutes and solving them took about 26 seconds; graph enumeration adds further
time. A complete run should be expected to take tens of minutes, depending on
the machine.

Audit the completed output and rerun every query:

```bash
.venv/bin/ctln-audit-search \
  --output-dir /tmp/ctln-exhaustive-v1 \
  --rebuild-space \
  --replay-solver
```

The audit always regenerates queries, which repeats the symbolic construction
work. Replay checks continue through all 875 files. The JSON summary separates
logical `sat`/`unsat` contradictions from inconclusive mismatches involving
`unknown`; the command exits with status 4 or 5, respectively.

Run the test suite:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

## Code guide

Reading the files in this order follows the mathematical argument.

- `src/ctln_search/graph.py` defines directed graphs, the 20-bit encoding, and
  the conversion from source-target edges to the receiving-row nonedge matrix.
- `src/ctln_search/filters.py` implements the reviewed graph exclusions. Every
  match records its theorem ID, assumptions, and a concrete graph witness.
- `src/ctln_search/isomorphism.py` relabels graphs and chooses canonical
  representatives by trying every vertex permutation.
- `src/ctln_search/model.py` constructs $B$, the five Cramer determinants,
  the characteristic polynomial, and the Routh--Hurwitz determinants using
  exact SymPy expressions.
- `src/ctln_search/feasibility.py` converts rational polynomials to SMT-LIB,
  calls Z3, and preserves `sat`, `unsat`, and `unknown` without numerical
  tolerances. For `sat`, it keeps the decision run separate from the model
  query and transcript.
- `src/ctln_search/search.py` joins the preceding pieces into the exhaustive,
  resumable search and writes the evidence files.
- `src/ctln_search/audit.py` regenerates graph-specific queries, checks a
  completed search directory, and can replay its saved queries.
- `src/ctln_search/__init__.py` only identifies the Python package.

The six files in `tests/` cover the seven substantive modules above. They test
edge orientation, exact matrix formulas, stability boundaries, theorem-filter
hypotheses, graph relabeling, solver serialization, positive and negative
controls, resumption, auditing, and artifact corruption.

## Mathematical documentation and evidence

- [`docs/PROJECT_SPEC.md`](docs/PROJECT_SPEC.md) states the theorem, its
  reduction to polynomial feasibility, and the exact status of the result.
- [`docs/CONVENTIONS.md`](docs/CONVENTIONS.md) fixes graph orientation, matrices,
  parameters, stability signs, and the nondegeneracy convention.
- [`docs/THEOREM_LEDGER.md`](docs/THEOREM_LEDGER.md) gives the precise approved
  results and assumptions used by the search.
- [`docs/CERTIFICATE_FORMAT.md`](docs/CERTIFICATE_FORMAT.md) describes the saved
  computation and its trust boundary.
- [`certificates/experiments/five_neuron_oracle_v1/`](certificates/experiments/five_neuron_oracle_v1/)
  contains the completed manifest, results, audit summary, checksum, and a
  compressed archive with all exact queries and raw responses.

The most important remaining caveat is that both search and audit are ordinary
Python programs. Query regeneration deliberately shares the reviewed model and
serialization implementation with the search, and Z3's `unsat` answers are not
converted into independently checked real-algebraic proof certificates.
~~~~

### 3(b) The graph exclusions, exactly, with their stated justification (docs/THEOREM_LEDGER.md, verbatim and complete)
~~~~
# Theorem ledger

Status: **reviewed for the version-1 exhaustive search**.

Only results used by the completed search appear here. An approved result may
be used only with every listed hypothesis and degeneracy restriction checked.

- **APPROVED-LITERATURE** means checked against the cited source.
- **APPROVED-INTERNAL** means proved in this file and tested in the project.
- **No global ND** means global CTLN nondegeneracy is not a hypothesis.

The target theorem assumes global nondegeneracy. A filter depending on that
assumption eliminates a graph only within the target scope. The direct
algebraic systems omit the assumption and hence test a stronger statement.

## Literature results

### CTLN-DEFINITION-01 — matrices and legal parameters

- **Status:** APPROVED-LITERATURE
- **Source:** [Geneson 2026](../literature/2607.21396v1.pdf), Section 2,
  equations (2)--(6), printed pp. 2--3.
- **Statement used:** For a loopless graph,
  `W[i,j] = -1 + epsilon` when `j -> i`, `W[i,j] = -1 - delta`
  otherwise, and `W[i,i]=0`. The legal range is
  `delta>0` and `0<epsilon<delta/(1+delta)`.
- **Coordinates used:** With receiving-row nonedge matrix `N`,
  `A=I-W=(1-epsilon)J+epsilon I+(delta+epsilon)N`. For
  `q=delta(1-epsilon)/epsilon`, legality is equivalent to
  `0<epsilon<1` and `q>1`.
- **Drive:** A positive uniform drive only rescales a fixed point and does not
  change support or stability.
- **Degeneracy:** No global ND.

### PARAM-CUTOFF-05 — all-graph parameter cutoff

- **Status:** APPROVED-LITERATURE
- **Source:** Geneson 2026, Theorem 4.1, printed pp. 9--11.
- **Statement used:** For a nonzero `m`-by-`m` zero-one nonedge matrix, if the
  associated legal CTLN matrix is permitted and positive stable, then

  $$
  q<m-2-\frac{m-3}{2}\varepsilon.
  $$

- **Five-neuron consequence:** Every nonclique candidate satisfies
  `q<3-epsilon`; the complementary region needs no graph-specific algebra.
- **Degeneracy:** No global ND for the algebraic statement.

### DOMINATION-ACTIVE-01 — graphical domination excludes full support

- **Status:** APPROVED-LITERATURE
- **Source:** [Curto 2025](../literature/2510.05098v1.pdf), Definition 1.1,
  Corollary 3.6, and Lemma 3.8, printed pp. 2 and 20--22.
- **Graph hypotheses:** Distinct vertices `j,k` satisfy `j->k`, not `k->j`,
  and, for every other `i`, `i->j` implies `i->k`.
- **Conclusion used:** No fixed point can have `x_j>0`; therefore the graph
  cannot support a full-support fixed point.
- **Network hypotheses:** Positive gCTLN edge/nonedge parameters and uniform
  drive as in the source. Legal CTLNs satisfy them.
- **Degeneracy:** No global ND.

### SPARSE-STABILITY-01 — small maximum in-degree

- **Status:** APPROVED-LITERATURE under target global ND.
- **Source:** [Curto--Geneson--Morrison 2023](../literature/1909.02947v3.pdf),
  Theorems 1.7 and 3.3, printed pp. 8 and 16--17.
- **Statement used:** A graph on `m` active vertices cannot support a stable
  fixed point if its maximum in-degree is less than `m/2`.
- **Five-neuron predicate:** Maximum in-degree at most two.
- **Degeneracy:** The paper adopts global nondegeneracy as a standing
  assumption in Section 1.5.

### NEAR-CLIQUE-STABILITY-01 — clique of size four

- **Status:** APPROVED-LITERATURE under target global ND.
- **Source:** Curto--Geneson--Morrison 2023, Theorem 1.8 and Section 4, printed
  pp. 8 and 18--19.
- **Statement used:** A nonclique graph containing a bidirected clique of size
  one less than the support cannot support a stable fixed point.
- **Five-neuron predicate:** A bidirected clique on four vertices in a
  nonclique full graph.
- **Degeneracy:** The paper's standing global-ND assumption applies.

### ORIENTED-NO-SINK-01 — oriented graph without sinks

- **Status:** APPROVED-LITERATURE under target global ND.
- **Source:** [Morrison--Degeratu--Itskov--Curto 2023](../literature/1605.04463v4.pdf),
  Theorem 3.4, printed p. 8.
- **Statement used:** A nondegenerate legal CTLN whose graph is oriented and has
  no sinks has no stable fixed point.
- **Predicate:** No vertex pair is bidirected, and every vertex has positive
  out-degree.

### ROUTH-HURWITZ-05 — exact positive-stability test

- **Status:** APPROVED-LITERATURE
- **Source:** [Holtz 2005](../literature/math_0512591v1.pdf), Theorems 2--3,
  printed pp. 2--3.
- **Statement used:** For
  `p(s)=det(sI+B)=s^5+a1*s^4+...+a5`, all eigenvalues of `B` have positive real
  part exactly when all leading principal minors of

  ```text
  [a1  a3  a5   0   0]
  [ 1  a2  a4   0   0]
  [ 0  a1  a3  a5   0]
  [ 0   1  a2  a4   0]
  [ 0   0  a1  a3  a5]
  ```

  are positive.
- **Convention check:** The roots of `det(sI+B)` are the negatives of the
  eigenvalues of `B`; left-half-plane stability of the polynomial is positive
  stability of `B`.
- **Degeneracy:** Strict positivity excludes imaginary-axis roots; no global
  ND.

### DIGRAPH-COUNT-05 — unlabeled graph-count regression

- **Status:** APPROVED-LITERATURE for regression only.
- **Source:** [Curto--Langdon--Morrison 2019](../literature/1902.10270v3.pdf),
  Table 1, printed p. 6.
- **Statement used:** There are 9,608 isomorphism classes of loopless directed
  graphs on five vertices.
- **Restriction:** This known count is a check on canonicalization, not the
  proof of labeled graph coverage.

## Project derivations

### TARGET-SCOPE-01 — the full clique is outside the claim

- **Status:** APPROVED-INTERNAL
- **Contract:** The desired implication assumes a nonclique graph. The full
  bidirected clique is recorded as outside the target, not as an eliminated
  counterexample.

### FULL-SUPPORT-01 — equilibrium and stability criterion

- **Status:** APPROVED-INTERNAL
- **Contract:** A stable full-support fixed point exists exactly when `A=I-W`
  is positive stable and `A^-1 1>0`.
- **Proof:** At full support all thresholds are locally inactive, so the
  equilibrium equation is `Ax=theta*1` and the Jacobian is `-A`. Stability
  makes `A` positive stable and therefore invertible. Conversely, if `A` is
  positive stable and `A^-1 1>0`, then `theta*A^-1 1` is an interior fixed
  point with Hurwitz Jacobian `-A`.
- **Degeneracy:** No global ND. Active nonsingularity follows from stability.

### POLYNOMIAL-SCALING-01 — denominator-free matrix

- **Status:** APPROVED-INTERNAL
- **Contract:** Substituting
  `delta=q*epsilon/(1-epsilon)` and multiplying `A` by `1-epsilon` gives

  $$
  B=\varepsilon(1-\varepsilon)I
   +\varepsilon(1-\varepsilon+q)N+(1-\varepsilon)^2J.
  $$

  Since `1-epsilon>0`, `A` and `B` have the same positive-stability status and
  their inverse-times-one vectors have the same entrywise signs.

### CRAMER-POSITIVITY-01 — polynomial fixed-point inequalities

- **Status:** APPROVED-INTERNAL
- **Contract:** If real `B` is positive stable, then `det(B)>0`, and
  `B^-1 1>0` exactly when all five determinants obtained by replacing one
  column of `B` by `1` are positive.
- **Proof:** Real eigenvalues of a positive-stable real matrix are positive;
  nonreal eigenvalues occur in conjugate pairs with positive product. Hence
  `det(B)>0`, after which Cramer's rule gives the claim.
- **Warning:** Without positive stability, the replacement determinants must
  be compared with the sign of `det(B)`.

### HURWITZ-LAST-MINOR-01 — reduced stability list

- **Status:** APPROVED-INTERNAL
- **Contract:** For a monic degree-`n` polynomial with constant coefficient
  `a_n`, `Delta_n=a_n*Delta_(n-1)`. Consequently strict positivity of all
  Hurwitz minors is equivalent to positivity of
  `Delta_1,...,Delta_(n-1),a_n`.
- **Proof:** Cofactor expansion along the final column of the Hurwitz matrix.
  Strict positivity of `Delta_(n-1)` gives the equivalence without admitting
  boundary zeros.

### UNIFORM-CODEGREE-ONE-01 — uniform in-degree `m-2`

- **Status:** APPROVED-INTERNAL
- **Origin:** Corrected and independently checked from
  [Notes on Uniform Degree Conjecture](../literature/uniform_degree_conjecture%20%281%29.pdf),
  Section 2. The local proof, not the unpublished note by itself, is approved.
- **Contract:** A nonclique graph on `m>=2` vertices with uniform in-degree
  `m-2` cannot have a positive-stable full active matrix at legal parameters.
- **Proof:** Its nonedge matrix has one `1` in each row and zero diagonal.
  Following the unique entry from each row produces a cycle of length
  `ell>=2`. On the corresponding row-vector subspace, right multiplication by
  the nonedge matrix contains the cyclic permutation matrix. Choose a
  nontrivial `ell`th root of unity `zeta` with `Re(zeta)<=-1/2`. Its row
  eigenvector has coordinate sum zero, so the `J` term vanishes. The active
  matrix has eigenvalue
  `epsilon+(delta+epsilon)zeta`, whose real part is at most
  `(epsilon-delta)/2<0`; legality implies `delta>epsilon`.
- **Degeneracy:** No global ND.
- **Five-neuron predicate:** Uniform in-degree three and nonclique.

## Version-1 usage record

The first-match filter pass used the following entries, in application order:

| Theorem ID | Labeled masks |
|---|---:|
| `TARGET-SCOPE-01` | 1 |
| `DOMINATION-ACTIVE-01` | 898,560 |
| `UNIFORM-CODEGREE-ONE-01` | 1,024 |
| `SPARSE-STABILITY-01` | 55,991 |
| `NEAR-CLIQUE-STABILITY-01` | 295 |
| `ORIENTED-NO-SINK-01` | 60 |

The remaining 875 canonical classes were tested using
`FULL-SUPPORT-01`, `POLYNOMIAL-SCALING-01`, `CRAMER-POSITIVITY-01`,
`ROUTH-HURWITZ-05`, and `HURWITZ-LAST-MINOR-01` on the parameter triangle
supplied by `PARAM-CUTOFF-05`. Every exact SMT-LIB formula returned `unsat`.

For a formal publication, the most important human checks are the recent
preprint statements `PARAM-CUTOFF-05` and `DOMINATION-ACTIVE-01`, the standing
global-ND scope of the 2023 filters, and the internal
`UNIFORM-CODEGREE-ONE-01` proof.
~~~~

### 3(c) Adjacency / matrix / parameter conventions (docs/CONVENTIONS.md, verbatim and complete)
~~~~
# Conventions

Status: **normative for the version-1 exhaustive search**.

This file fixes meanings that must not be inferred independently in code.
Statements marked **Literature fact** are taken from the cited source;
statements marked **Project derivation** are proved here from those definitions.

## 1. Vertices and directed graphs

- Code uses vertices `0, ..., n - 1`. Papers generally use `[n] = {1, ..., n}`.
- Graphs are loopless directed graphs. For each unordered pair, neither, either,
  or both directed edges may be present.
- `u -> v` means an edge whose source is `u` and whose target is `v`.
- A code-facing adjacency array is source-major:
  `edge[u, v] == 1` if and only if `u -> v`. Thus rows are sources and columns
  are targets.
- A clique is a vertex set with both directed edges between every distinct pair.
  A full graph is nonclique exactly when at least one ordered edge is absent.

### Labeled graph masks

For `n` vertices, bit positions enumerate ordered pairs `(source, target)` with
`source != target` in lexicographic order. Bit `k` is one exactly when the
corresponding directed edge is present. Thus a five-neuron graph uses 20 bits,
and the full bidirected clique has mask `2^20 - 1`. This is a project convention
and is tested by exact graph-to-mask round trips; it is distinct from the
receiving-row matrix orientation below.

**Literature fact.** The CTLN papers write matrix rows as receiving vertices and
columns as sending vertices: `W[i,j]` is the influence `j -> i`. See
[Geneson 2026](../literature/2607.21396v1.pdf), equation (2), printed p. 2.

## 2. Source, proper source, sink, and isolated vertex

- A **source** has in-degree zero.
- A **proper source** has in-degree zero and at least one outgoing edge.
- A **sink** has out-degree zero.
- An **isolated vertex** has both in-degree and out-degree zero. It is a source
  and a sink, but it is not a proper source.

**Literature fact.** These definitions are stated in
[Curto--Geneson--Morrison 2018](../literature/1804.00794v2.pdf), Section 4.2,
printed p. 18. Its Rule 6 excludes supports containing a *proper* source; Rule
12 treats an isolated vertex differently.

**Correctness warning.** A filter for “has a source” is unsound if it also
removes isolated vertices. For full support, the justified filter is “has a
proper source.”

## 3. CTLN dynamics and matrices

For a loopless directed graph `G`, the CTLN is

$$
  \dot x_i=-x_i+\left[\sum_j W_{ij}x_j+\theta\right]_+,
  \qquad \theta>0,
$$

with

$$
W_{ij}=\begin{cases}
0&i=j,\\
-1+\varepsilon&j\to i,\\
-1-\delta&j\not\to i.
\end{cases}
$$

The legal parameter range is

$$
  \delta>0,\qquad 0<\varepsilon<\frac{\delta}{1+\delta}.
$$

These are **literature facts** from
[Geneson 2026](../literature/2607.21396v1.pdf), equations (2)--(3), printed p. 2.

Define the receiving-row nonedge matrix

$$
  N_{ii}=0,\qquad N_{ij}=1\iff j\not\to i\quad(i\ne j).
$$

Thus `N[i,j] = 1 - edge[j,i]` off the diagonal: `N` is aligned with `W`, not
with the source-major adjacency array. Let `J = 11^T` and

$$
  A=I-W.
$$

Then

$$
  A=(1-\varepsilon)J+\varepsilon I+(\delta+\varepsilon)N.
$$

This identity is both a **literature fact** (Geneson 2026, equations (5)--(6))
and a required convention test.

## 4. Parameter coordinates

Use

$$
 t=\frac{1-\varepsilon}{\varepsilon},\qquad
 \beta=1+\frac{\delta}{\varepsilon},\qquad
 q=\frac{\delta(1-\varepsilon)}{\varepsilon},
$$

and

$$
 K=\varepsilon^{-1}A=I+\beta N+tJ.
$$

**Literature fact.** The legal condition is equivalent to `q > 1`; see
[Geneson 2026](../literature/2607.21396v1.pdf), equations (5)--(6), printed
pp. 2--3.

For exact polynomial calculations, the preferred coordinates are
`(epsilon, q)` with

$$
  0<\varepsilon<1,\qquad q>1,
  \qquad \delta=\frac{q\varepsilon}{1-\varepsilon}.
$$

Clear the positive denominator with

$$
  B=(1-\varepsilon)A
   =\varepsilon(1-\varepsilon)I
    +\varepsilon(1-\varepsilon+q)N
    +(1-\varepsilon)^2J.
$$

The expression for `B` is a **project derivation**. Multiplication by the
positive scalars `epsilon` or `1 - epsilon` preserves entrywise positivity of
the inverse solution up to a positive factor and preserves the signs of real
parts of eigenvalues. Consequently `A`, `K`, and `B` are interchangeable for
the two full-support tests when their corresponding right-hand-side scaling is
handled consistently.

## 5. Full support and stability

“Full support” means `x_i > 0` for every vertex. There are no off-neuron or
outside-support inequalities.

We reserve:

- **positive stable** for a real matrix whose eigenvalues all have positive
  real part;
- **Hurwitz** for a real matrix whose eigenvalues all have negative real part;
- **stable full-support fixed point** for the dynamical equilibrium.

At a full-support point the threshold map is locally the identity, so the
Jacobian is `-A`. Hence the following is a **project derivation**:

$$
\begin{split}
&\text{a stable full-support fixed point exists}\\
&\quad\Longleftrightarrow\quad
A\text{ is positive stable and }A^{-1}\mathbf1>0\\
&\quad\Longleftrightarrow\quad
B\text{ is positive stable and }B^{-1}\mathbf1>0.
\end{split}
$$

The value of `theta > 0` only multiplies the equilibrium and does not affect
support or stability. This is stated explicitly in Geneson 2026, printed p. 2.

## 6. Degeneracy

The project distinguishes three notions.

1. **Global CTLN nondegeneracy:** for every nonempty vertex subset `rho`, the
   principal matrix `A_rho` and every uniform-input column-replacement
   determinant are nonzero. This is the standing assumption in several older
   CTLN papers; see
   [Curto--Geneson--Morrison 2023](../literature/1909.02947v3.pdf), Definition
   2.1, printed p. 11.
2. **Active nonsingularity:** the matrix for the candidate support is
   nonsingular.
3. **Boundary degeneracy:** a determinant or strict inequality used in a
   support test vanishes at a parameter value.

The target theorem **does assume global nondegeneracy**, matching the standing
convention in the CTLN literature. A certificate or theorem filter that relies
on this premise must say so explicitly.

**Project derivation.** Positive stability of the full active matrix implies
that zero is not an eigenvalue, so active nonsingularity is automatic. Strict
positivity of `A^{-1}1` excludes zero Cramer coordinates. Consequently, a
direct full-support infeasibility calculation may ignore global
nondegeneracy and prove a stronger statement. Global nondegeneracy remains
relevant because literature filters may reason about proper supports and
because a feasible point from the stronger search is a counterexample to the
project theorem only if all proper-support nondegeneracy conditions also hold.

Literature theorems under a global-nondegeneracy convention are approved for
use when their remaining graph and parameter hypotheses are verified. The
filter metadata and manifest scope counts record that they cover the globally
nondegenerate target, not all CTLNs.

## 7. Exact inequalities and boundaries

All domain, positivity, and stability inequalities are strict. Boundary points
must be assigned by proved weak/strict inequalities, never by floating-point
tolerances.

[Geneson 2026](../literature/2607.21396v1.pdf), Theorem 4.1, proves that for a
nonzero five-by-five nonedge matrix, permission and positive stability imply

$$
  q<3-\varepsilon.
$$

Thus the only parameter region requiring per-graph feasibility work is the open
triangle

$$
  0<\varepsilon<1,\qquad 1<q<3-\varepsilon.
$$

The complementary boundary `q >= 3 - epsilon` is covered by that theorem,
without a global-nondegeneracy assumption.

## 8. Graphical domination

For distinct vertices `j,k`, vertex `k` graphically dominates `j` in the full
graph when

1. for every `i` other than `j,k`, `i -> j` implies `i -> k`;
2. `j -> k`; and
3. `k -/> j`.

This is **literature fact** from
[Curto 2025](../literature/2510.05098v1.pdf), Definition 1.1, printed p. 2.
For a CTLN it implies input domination and therefore no fixed point can have
`x_j > 0` (Lemma 3.8 and Corollary 3.6, printed pp. 20--22). Those two results
do not require global nondegeneracy.

**Project derivation.** Every proper source is dominated: choose one of its
out-neighbors as `k`; the incoming-edge implication is vacuous. Therefore a
proper source rules out full support for every legal parameter value. An
isolated source has no such out-neighbor and is not covered by this argument.
~~~~

### 3(d) Summary of Repo 2's method, region, convention, and evidence (my reading, with checks)

- **Claim (docs/PROJECT_SPEC.md):** for a globally nondegenerate CTLN on 5 vertices, if G is not a bidirected clique there is no stable fixed point with full support (all 5). "Nothing is claimed about proper supports, degenerate CTLNs, or networks with six or more neurons." Strictly weaker in scope than Repo 1's order <= 5 claim, and weaker in rigor (Z3 `unsat` is trusted, not a proof certificate).
- **Method:** (1) enumerate all 2^20 labelled loopless digraphs on 5 vertices; (2) first-match graph exclusions (table in the README, justifications in the ledger above); (3) canonicalize survivors by min mask over all 120 relabelings; (4) for each of the 875 canonical survivors build `B(eps,q) = eps(1-eps) I + eps(1-eps+q) N + (1-eps)^2 J` (= (1-eps)(I-W) with delta = q eps/(1-eps)); assert, in Z3 `QF_NRA`, the 4 domain inequalities (`eps>0, 1-eps>0, q-1>0, 3-eps-q>0`), the 5 Cramer numerators `>0`, and the 5 reduced Routh-Hurwitz conditions `Delta_1..Delta_4, a_5 > 0` (14 strict polynomial constraints); Z3 4.13.3 returned `unsat` on all 875 (30 s timeout; total solver time 26.0 s, max 0.23 s per query; symbolic construction 1575 s).
- **The exclusions, exactly (first-match counts, labelled graphs):** full 5-clique 1 (outside the claim); graphical domination 898,560 (all legal parameters, no global-ND assumption; Curto 2510.05098 Def 1.1 + Cor 3.6 + Lemma 3.8); uniform in-degree 3 on 5 nodes 1,024 (all legal parameters; "UNIFORM-CODEGREE-ONE-01", proved in their ledger from an unpublished note, flagged by them as needing human review); maximum in-degree <= 2 55,991 (CGM 2023 Thms 1.7/3.3, requires global ND); contains a bidirected 4-clique 295 (CGM 2023 Thm 1.8, requires global ND); oriented graph with no sink 60 (Morrison-Degeratu-Itskov-Curto, Thm 3.4, requires global ND). Total removed 955,931; survivors 92,645 labelled = 875 classes.
  Their definition of domination (README): "a vertex k dominates j when j -> k, not k -> j, and every other input to j is also an input to k"; ledger: "Distinct vertices j,k satisfy j->k, not k->j, and, for every other i, i->j implies i->k." This is word-for-word the content of Definition 1.1 in the PDF copy of 2510.05098 (s4).
- **Parameter region searched:** the open triangle `0 < eps < 1, 1 < q < 3 - eps` (q = delta(1-eps)/eps), "all-parameter satisfiability, not a grid or random sample". The complementary region `q >= 3 - eps` is NOT searched: it is discharged by Geneson's Theorem 4.1 (n=5 case of `q < m - 2 - (m-3)eps/2`), whose statement I confirmed in the PDF copy of 2607.21396 (s4) but whose proof I did not check. No global-nondegeneracy constraint in the direct query (so `unsat` is stronger for the 875 graphs); the 4 global-ND-dependent filters (56,346 labelled graphs) only eliminate within the nondegenerate scope.
- **Stated adjacency convention:** OURS. "A code-facing adjacency array is source-major: `edge[u, v] == 1` if and only if `u -> v`. Thus rows are sources and columns are targets." Masks: "bit positions enumerate ordered pairs (source, target) with source != target in lexicographic order. Bit k is one exactly when the corresponding directed edge is present ... the full bidirected clique has mask 2^20 - 1." Receiving-row nonedge matrix `N[i,j] = 1 - edge[j,i]` off the diagonal (their `N` = Geneson/Repo-1 `N`, i.e. the `J-I-A^T` of ours). `W[i,j]=-1+eps if j->i; -1-delta else`, legal iff `q>1`. So masks, edge lists and nonedge matrices in Repo 2 can be diffed against ours with no transposition; Repo 1's masks use a DIFFERENT encoding (nonedge bits, row-major, `N_ij`), do not mix them.
- **Are the unsat certificates and the 875 graphs in the repo, readable?** Yes for the graphs and the queries; no for proof objects.
  - `certificates/experiments/five_neuron_oracle_v1/` has `manifest.json` (all 875 canonical masks + orbit sizes + filter counts, plain JSON), `results.jsonl` (875 records: canonical mask, receiving-row nonedge matrix, orbit size, status, SHA-256 of the query, timings), `audit.json`, `audit-query-provenance.json`, `SHA256SUMS`, and `ctln-exhaustive-v1.tar.gz` (1.87 MB; expands to 29 MB: `coverage.jsonl` with the 92,645 labelled survivors and their relabeling permutations, 875 `queries/<mask>.smt2`, 875 `solver_outputs/<mask>.txt`). MY CHECK: `sha256sum -c SHA256SUMS` OK for all five files.
  - Each `.smt2` is a plain SMT-LIB 2 `QF_NRA` problem in the two real variables epsilon and q with integer polynomial constraints (only options present in any file: `(set-logic QF_NRA)` and `(set-option :timeout 30000)`). Each `solver_outputs/*.txt` is the single word `unsat` (875/875). So anyone with Z3 can replay; but there are NO proof certificates (Z3 gives no checkable real-algebraic proof here), as their `CERTIFICATE_FORMAT.md` states: "They do not provide an independently checked real-algebraic proof of each `unsat` result."
  - The 875 graphs are directly readable from `manifest.json` / `results.jsonl` (e.g. canonical mask 4375, orbit size 20, `nonedge_matrix_receiving_row = [[0,0,0,0,1],[0,0,1,1,1],[0,1,0,1,1],[0,1,1,0,1],[1,1,1,1,0]]`, status `unsat`).
- **MY CHECKS (independent of their code; stated definitions only):** re-enumerating all 2^20 masks with their stated bit convention and applying the filters in their stated first-match order I get exactly 1 / 898,560 / 1,024 / 55,991 / 295 / 60, 92,645 survivors, 875 canonical classes (min mask over 120 relabelings), and every canonical mask and every orbit size equals their manifest (symmetric difference 0). For all 875 records, `nonedge_matrix_receiving_row` equals `1 - edge[j][i]` decoded from the canonical mask under their stated convention. This confirms the bookkeeping and that the stated definitions give their numbers; it does NOT validate the soundness of the five cited filter theorems, nor re-solve any query.
- **Soundness flags for anyone reusing their filters:** (i) four filters hold only under global nondegeneracy; (ii) `UNIFORM-CODEGREE-ONE-01` rests on a locally written proof derived from an unpublished note; (iii) `PARAM-CUTOFF-05` is a one-month-old preprint theorem; (iv) their own CONVENTIONS warn that "has a source" is unsound as a domination filter (isolated vertices) and the justified one is "has a proper source".


## 4. Side findings: primary-source text obtained from the PDF copies bundled in Repo 2 (`literature/`)

CAUTION: these are copies inside an untrusted repo, text-extracted with `pdftotext -layout`; PDF Producer is `pikepdf 8.15.1` (re-processed), so I cannot certify they are byte-identical to
arXiv. The first-page header text reads `arXiv:2607.21396v1 [math.CO] 23 Jul 2026` and `arXiv:2510.05098v1 [q-bio.NC] 6 Oct 2025`. SHA-256: 2607.21396v1.pdf `1bb6d9a37513a36790d6d01125e59a2741c88410c0e7dc896a917f009e466fe3`;
2510.05098v1.pdf `4eee2b0f9b4c8bf6e9e136258764b6f3faf32eb6b02e13282e7463363c3a4b8d`; 2301.12638v1.pdf `39f7410c87506821da072e6cca977a6257f057e692ef1dd160d388c3748ad0bf` (present, NOT read by me).
The directory also holds 1409.1534v1, 1605.04463v4, 1804.00794v2, 1902.10270v3, 1909.02947v3, math_0512591v1 and "uniform_degree_conjecture (1).pdf" (an unpublished note). Treat anything below as "confirmed against a copy", to be re-confirmed against arXiv by a trusted fetch.

**Geneson 2607.21396v1 (abstract and s2, verbatim, pdftotext):**
> "For every sufficiently small eps > 0, the CTLN defined by one fixed graph at delta = 29 eps/25 is nondegenerate and has a stable fixed point with nonclique full support. Its values of q = delta(1 - eps)/eps tend to 29/25. In the complementary direction, for any CTLN on n >= 3 vertices, we prove that in the parameter range  q >= n - 2 - (n-3)/2 * eps,  no nonclique support can satisfy both the fixed-point positivity and linear stability conditions. ... In particular, this holds when eps <= delta/(delta + n - 2)."
> Legal range (3): `delta > 0, 0 < eps < delta/(1+delta)`. (5): `t = (1-eps)/eps, beta = 1 + delta/eps, q = delta(1-eps)/eps`; (6): `A_sigma := (I_n - W)_sigma, K := eps^{-1} A_sigma = I_m + beta N + t J_m`; "The legal condition is equivalent to q > 1, and also to t(beta - 2) > 1."
> Nonedge matrix: "N_aa = 0, N_ab = 1 <=> v_b -/-> v_a (a != b)".
> Theorem 3.2. Let N = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]] (19). There is an eps_0 > 0 such that, for every 0 < eps < eps_0 and every theta > 0, the CTLN with full nonedge matrix N and delta = (29/25) eps is legal and nondegenerate, its full support is permitted, and the corresponding fixed point is stable. The full support is not a clique, and q = (29/25)(1 - eps) -> 29/25 (eps -> 0).
> Theorem 4.1. Let delta > 0 and 0 < eps < delta/(1+delta), and define t, beta, q by (5). Let m >= 3 and let N in {0,1}^{m x m} be nonzero with zero diagonal. If K = I_m + beta N + t J_m is permitted and stable, then q < m - 2 - (m-3)/2 * eps. (24) Consequently ... every CTLN on a graph with n >= 3 vertices has the algebraic target-free clique classification whenever q >= n - 2 - (n-3)/2 * eps. (25)
So the "q ~ 29/25" and "q >= n - 2 - (n-3) eps/2" you were told secondhand are confirmed IN THIS COPY, and `q` there is `delta(1-eps)/eps`, same as CLAUDE.md. Note 29/25 is only the limit along ONE ray chosen by Geneson (his Theorem 3.2 says nothing about the lower end of the stable region); the ray `delta/eps = 29/25` is `u = 4/25` in Repo 1's coordinates.

**Curto 2510.05098v1 (Definition 1.1 and Theorems 1.4, 1.6, verbatim, pdftotext):**
> Definition 1.1. Let j, k in [n] be vertices of G. We say that k graphically dominates j in G, and write k > j, if the following two conditions hold: (i) For each vertex i in [n] \ {j, k}, if i -> j then i -> k. (ii) j -> k and k -/-> j. If there exists a k that graphically dominates j, we say that j is a dominated node (or dominated vertex) of G. If G has no dominated nodes, we say that it is domination free.
> Theorem 1.4 (Theorem 1). Suppose j is a dominated node in a directed graph G. Then the fixed points of a gCTLN constructed from G satisfy FP(G) = FP(G|_{[n]\j}), for any choice of gCTLN parameters {eps_i, delta_i, theta}.
> Theorem 1.6 (Theorem 2). Let G be a directed graph on n nodes. Then the reduced graph G~ is unique.   [Definition 1.5: G~ is the reduced graph if (1) G~ = G|_tau is an induced subgraph, (2) G~ is domination free, (3) G~ can be obtained from G by iteratively removing dominated nodes.]
This is the real definition for replacing the `graphs.dominates()` placeholder (in this copy). It matches Repo 2's implementation. Orientation: here `i -> j` means edge from i to j, i.e. `A[i][j] = 1` in our convention: `dominates(k, j)` iff `A[j][k] and not A[k][j] and all(A[i][k] for i not in (j,k) if A[i][j])`.
MY CHECK (regression targets for a replacement `dominates`, computed from Definition 1.1 exactly, `domfree_counts.log`): number of DOMINATION-FREE digraphs, labelled / up to isomorphism:
n=2: 2 / 2 (of 4 / 3); n=3: 16 / 6 (of 64 / 16); n=4: 716 / 47 (of 4096 / 218); n=5: 150,016 / 1,434 (of 1,048,576 / 9,608). (n=5: 898,560 labelled graphs have a dominated node, matching Repo 2's count.)
The preservation lemma is Theorem 1.4 (FP(G) = FP(G restricted to [n]\j) for every gCTLN parameter choice, so also CTLNs); I did not read its proof.

Further statements in the same copy that matter for pruning (verbatim): "Corollary 3.6. Suppose k input dominates j in the TLN (W, b). Then there can be no fixed point x* of (W, b) with x*_j > 0." (no nondegeneracy hypothesis);
"Theorem 3.7. Let (W, b) be a nondegenerate TLN on n nodes, and suppose k input dominates j. Then, FP(W, b) = FP(W|_{[n]\j}, b|_{[n]\j})."; "Lemma 3.8. Consider a directed graph G on n nodes. Suppose k graphically dominates j for some j, k in [n]. Then, for any gCTLN or E-I TLN with graph G, k input dominates j."
CAUTION for the soundness gate in CLAUDE.md s6: what these results say is about the SET of fixed points (supports/values). For a FULL-SUPPORT sweep that is already enough to discard any graph with a dominated node (Cor 3.6: no fixed point has x_j > 0, so the full vertex set cannot be a support), with no stability argument needed. For pruning a graph to its reduced graph while tracking STABILITY of fixed points, the quoted theorems do not state preservation of Jacobian eigenvalues; that would need its own lemma (plausible via the block-triangular structure when the inactive input is strictly negative, but I did not verify it and the paper text I read does not assert it).


## 5. What I could NOT verify (explicit gaps)

1. **The proofs.** I verified polynomial inputs and sampled conclusions, but did not re-derive or formally check: Repo 1's Theorem 3.1 (all-parameter n <= 5 exclusion; I did not write an independent obstruction checker and did not run theirs), Theorem 4.3 (all-r family; only r = 0..3 pointwise, nondegeneracy r = 0,1), Theorem 5.1 (uniqueness/ordering arguments; only polynomial identity and 1,353 exact grid points), Proposition 6.1, or the parity-mod-2 nondegeneracy argument for r >= 2.
2. **Execution claims.** Nothing from either repo was run. Their logs (`verification/*.log`, `summary.json`) say PASS on Python 3.13.5; the GitHub Actions run, `REPLAY_SUMMARY.json` (referenced but not in the repo), the Zenodo DOI 10.5281/zenodo.23089682 and `MANIFEST.sha256` were not checked. GitHub API metadata (issues, forks, stars, CI status, who else has contacted these authors) is unavailable: `gh api` returns HTTP 403 "GitHub access to this repository is not enabled for this session".
3. **Repo 2's Z3 results.** I did not re-run Z3 on any query and did not independently rebuild the 875 polynomial systems; `unsat` is Z3's word and there are no proof objects. I only checked the graph bookkeeping (filter counts, classes, orbit sizes, nonedge matrices) and file integrity (`sha256sum -c`).
4. **The cited filter theorems** (CGM 2023 Thms 1.7/1.8/3.3, MDIC Thm 3.4, Geneson Thm 4.1 proof, Repo 2's own UNIFORM-CODEGREE-ONE-01 proof): not read or checked, apart from Geneson Thm 4.1's and Curto's Def 1.1/Cor 3.6/Lemma 3.8/Thm 3.7 statements as printed in the bundled PDF copies.
5. **PDF copies in Repo 2's `literature/`:** not confirmed identical to arXiv (re-processed with pikepdf). Re-fetch from arXiv (WebFetch) before citing. `2301.12638v1.pdf` (Curto-Morrison, the 14 open questions) was NOT read in this task.
6. **Global nondegeneracy of the seed on its whole phase region** is not claimed by Repo 1 either (their C6 "Does not assert that all other subsets are nondegenerate at every such point"); I only verified nondegeneracy at the specific points (1/12, 11/120) and the family r = 0,1 points. For the dynamical "counterexample" label at a general point of the region, nondegeneracy (and hence strict outside inequalities, vacuous here) still has to be checked per point.
7. **Authorship / contact.** Repo 1's creator is "Anonymous" (AI-assisted; maintainer Ian Pitchford, "Evidence Press"); Repo 2 commits are by "Maximiliano Vazquez". I did not look for e-mail addresses or whether either has contacted Geneson/Curto/Shiu (CLAUDE.md s12).

## 6. Files

- Report: `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/repos/overlap_report.md`
- Clones (read-only data): `.../scratchpad/repos/ipitchford_six/repo`, `.../scratchpad/repos/chickennungets_ctln/repo`; extracted Repo-2 archive: `.../scratchpad/repos/chickennungets_ctln/archive_extract/ctln-exhaustive-v1/`
- Extracted PDF text: `.../scratchpad/pdftext/{six_is_minimal,2607.21396v1,2510.05098v1}.txt`
- My scripts and logs (all outside the clones, run as `python3 -I` from `/tmp`): `.../scratchpad/scripts/`: `check_seed.py` (+`check_seed.log`), `grid.py`, `near_one.py` (+log), `family_nd.py` (+log), `quotient_check.py`, `repo2_counts.py` (+log), `domfree_counts.py` (+log), `small_order_sanity.py` (+log), `helpers.py` (exact Fraction helpers: LeVerrier characteristic polynomial, Hurwitz determinants, Gauss solve), `graphs_out.json` (converted seed and r=1 matrices).
