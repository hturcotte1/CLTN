# Enumerating the digraphs on K_8 by their bidirected subgraph

Purpose: the n = 8 limit prefilter enumerates, for every underlying simple graph G on 8 vertices with 15 to 28
edges, the digraphs with underlying graph G and at least 29 arcs (nauty `directg -e29:56`), one representative
per isomorphism class, and runs `filter8` on each.  For G = K_8 the automorphism group has order 40,320 and
`directg`'s isomorph rejection runs at a few hundred digraphs per second, so the 575,009,339 classes would take
days.  The lemma below replaces that single enumeration by 12,345 enumerations with small automorphism groups.
It is a re-enumeration of exactly the same set of isomorphism classes, not a prune.

## Setting

Vertex set V = {0, ..., n-1}, K = K_n, P = the set of the C = n(n-1)/2 pairs.  A digraph D with underlying graph
K assigns to every pair {u,v} one of three states: the single arc u -> v, the single arc v -> u, or both arcs.
Write B(D) for the graph on V whose edges are the pairs in state "both" (the bidirected subgraph), H(D) = K - B(D)
for the complement graph, and O(D) for the orientation of H(D) given by the single arcs.  D is determined by the
pair (B(D), O(D)), and conversely every pair (B, O) with B a graph on V and O an orientation of K - B is such a
digraph.  Permutations pi of V act on digraphs, graphs and orientations in the obvious way.

## Lemma (decomposition by bidirected subgraph)

Let R be a set of representatives of the isomorphism classes of graphs on V, one per class.  Then the map

    (B, [O])  ->  class of the digraph (B, O),     B in R,  [O] an orbit of Aut(B) on the orientations of K - B,

is a bijection from the disjoint union over B in R of the orbit sets onto the set of isomorphism classes of
digraphs with underlying graph K.

Proof.
(i) Equivariance.  For every permutation pi of V, B(pi D) = pi(B(D)) and O(pi D) = pi(O(D)), because pi maps
the pairs in state "both" to pairs in state "both" and single arcs to single arcs.  In particular isomorphic
digraphs have isomorphic bidirected subgraphs.
(ii) Surjectivity.  Given D, pick B in R and a permutation sigma with sigma(B(D)) = B (B(D) is isomorphic to
exactly one member of R).  Then D' = sigma D is isomorphic to D and B(D') = B by (i), so the class of D is the
image of (B, [O(D')]).
(iii) Injectivity.  Let D1 = (B, O1) and D2 = (B', O2) with B, B' in R be isomorphic, D2 = pi D1.  By (i),
B' = pi(B), so B and B' are isomorphic members of R, hence B = B' and pi is in Aut(B); and O2 = pi(O1), so
[O1] = [O2].  Conversely, if O2 = pi(O1) with pi in Aut(B) then pi D1 = (pi(B), pi(O1)) = (B, O2) = D2.
(iv) Aut(K - B) = Aut(B): a permutation of V preserves the edge set of B iff it preserves its complement in P.
So the orbits in the statement are the orbits of Aut(H) on the orientations of H = K - B.  QED.

What nauty's `directg -o H` lists is, in the words of its documentation, the orientations of H with
"isomorphic directed graphs derived from the same input suppressed".  That is the same orbit set: an
isomorphism between two orientations of the same labelled graph H maps H onto itself, so it is an element
of Aut(H), and conversely every element of Aut(H) maps orientations of H to isomorphic orientations.  The
output keeps the input labelling (vertex i of the output is vertex i of H), which is what addarcs relies on
when it adds the arcs of B on the same labels.

## Arc count and the 29-arc bound

|arcs(D)| = C + |E(B(D))| (one arc per pair plus one more per bidirected pair).  For n = 8, C = 28, so the
prefilter's restriction 29 <= arcs <= 56 is exactly "B nonempty", and the enumeration runs over the 12,345
graphs on 8 vertices with 1 to 28 edges (`geng 8 1:28`; 12,346 graphs on 8 vertices minus the empty one).
B = K_8 contributes the single complete bidirected digraph (56 arcs), as in the reference enumeration.

## Counting

By Burnside the number of classes equals sum over B in R, B nonempty, of (1/|Aut(B)|) sum over pi in Aut(B) of
the number of orientations of K - B fixed by pi.  Independently, Burnside over S_8 acting on the 28 pairs with
3 states per pair, minus the tournaments, gives 575,009,339 (`n8/k8decomp/polya_count.py 8`; the run's
`graph_counts.json` entry for `G~~~~{` has the same value from `count_restricted.py`).  The enumeration is
accepted only if the batch totals sum to exactly this number (`n8/check_accounting.py`).

## What the implementation relies on, and how each piece is checked

1. `directg -o H` emits exactly one representative per Aut(H)-orbit of orientations, with the input
   labelling (nauty semantics, see above).  Checked by `n8/k8decomp/validate_decomp.sh n`: for n = 4, 5, 6
   and 7 the canonical forms (labelg) of the decomposed stream and of the reference `directg -e(C+1):(2C)
   K_n` stream were compared as sorted multisets, with a duplicate check, and found equal (38, 570, 21,424
   and 2,141,832 classes; the script prints the counts, and they were compared by hand with the Burnside
   numbers from `polya_count.py`, which agree).  At n = 7 the decomposition took 7 s against 956 s for the
   reference enumeration.  The count-only dry run `count_all_B.py 8` over all 12,345 B sums to exactly
   575,009,339 (`k8_perB_counts.tsv`), and an independent per-B Burnside count over Aut(B) (referee check)
   agrees with the dry run for every B.
   Reproduction note: pass graph6 strings to printf as data, `printf '%s\n' "$H"`, never as the format
   string, since graph6 can contain backslashes.
2. `n8/k8decomp/addarcs B` ORs the bidirected arcs of B into each orientation (digraph6 in, digraph6 out).
   Unit-tested against a pure-Python reference (`test_addarcs.py`, 3,245 cases including nauty round trips).
3. `filter8` is unchanged; its verdict is an isomorphism invariant (characteristic polynomials and the
   domination condition are permutation-invariant), and `validate_decomp.sh` also compares the filter's output
   lines on both streams after canonical labelling (V2).
4. Accounting: the K_8 batches are registered in `n8/jobs/jobs.tsv` and `graph_counts.json` as parts of the
   graph `G~~~~{` with `exact: false`, so `check_accounting.py` sums them against the exact count exactly as
   it does for every `directg -s` split graph.  In addition each batch row carries the batch's exact expected
   count (the sum of the dry-run counts of its B), which the driver checks job by job (a mismatch leaves no
   `.done`), and `check_k8decomp.py` compares the per-B line counts written by addarcs (`job<j>.genlog`) with
   the dry-run ledger B by B after the run.

## Relative form: any underlying graph G (n8/orbitenum)

The same argument works for an arbitrary underlying graph G on V, with Aut(G) in place of S_n.  Fix the
labelled graph G.  A digraph with underlying graph exactly G is a pair (B, O) with B a subset of E(G) (the
edges carrying both arcs) and O an orientation of G - B.  If pi D1 = D2 for two such digraphs then pi maps G
onto G, so pi is in Aut(G), and B(D2) = pi(B(D1)), O(D2) = pi(O(D1)).  Hence, with R_G a set of
representatives of the orbits of Aut(G) on the subsets of E(G), the isomorphism classes of digraphs with
underlying graph isomorphic to G are in bijection with the pairs

    (B in R_G,  orbit of Stab(B) = {g in Aut(G) : g(B) = B} on the orientations of G - B),

by the same surjectivity and injectivity steps as in the lemma: a digraph D whose underlying graph is merely
isomorphic to G is first moved by some vertex permutation tau with tau(U(D)) = G, after which sigma ranges over
Aut(G) and takes B(tau D) to its representative in R_G; B = B' because R_G has one representative per
Aut(G)-orbit; pi lies in Stab(B).  Arcs = |E(G)| + |B|, so the
29..56-arc restriction is |B| >= 29 - |E(G)|.

`n8/orbitenum/orbitenum G min_arcs [p/P]` implements this directly, without directg: Aut(G) by brute force
over the n! vertex permutations; B is a representative iff it is the numerically least member of its
Aut(G)-orbit (edge bitmask, edges of G in lexicographic order); an orientation O of G - B is a representative
iff it is the least member of its Stab(B)-orbit (bit t = 1 iff the t-th free edge (u, v), u < v, is oriented
v -> u).  Least-element representatives give exactly one member per orbit.  Representatives B are numbered
in enumeration order and part p/P takes those with index = p (mod P), so the parts partition the output.
Validation: `validate_orbit.sh` (canonical-form multiset equality with directg's enumeration of the same
graph and the same arc bound, duplicate check, filter8 output equality) passed on K_5 - K_3, K_6 - K_3,
K_6 - 3K_2, K_6 - K_{1,3}, K_7 - K_3, K_7 - 3K_2, K_7 - K_{1,4} and K_7 - (K_2 + K_{1,3}) (up to
10,840,382 classes); the class counts equal directg's on 16 (graph, bound) cases; the count-only runs on the
four n = 8 graphs it replaces (K_8 - 4K_2, K_8 - K_3, K_8 - K_{1,5}, K_8 - (K_2 + K_{1,5})) reproduce the
exact Burnside counts of graph_counts.json (695,586,658; 1,179,317,832; 330,741,720; 84,866,528) in 12, 47,
13 and 3 seconds; the union of the p/P parts equals the unsplit output (K_7 - (K_2 + K_{1,3}), P = 3).
It refuses |Aut(G)| > 8192, so K_8 itself stays with the geng-based batches above.
