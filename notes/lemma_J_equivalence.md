# Lemma J: full-support data equivalence ("in-twin switching" and its generalisation)

Prover notes, 2026-10-13. WORK IN PROGRESS -- saved incrementally.
Conventions exactly CLAUDE.md §2 and notes/lemma_I_spectral.md §0:
`A[i][j] = 1` iff `i -> j`; `N := J - I - A^T` (`N[i][j] = 1` iff `i != j` and `j -/-> i`);
`C := I + beta N`, `K := C + tJ`, `a := adj(C) 1`, `d := det C`, `phi(mu) := det(mu I + C)`,
`omega(mu) := 1^T adj(mu I + C) 1 = det(mu I + C + J) - det(mu I + C)`.
`D(A) := (d, phi, omega, multiset{a_i})`, all as polynomials in `beta` (and `mu`).
Scratch scripts: /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/ (s1..s6 so far).

## 0. First finding: G2 is G1 with ONE edge moved between two in-twins

Relabel G2 by p = (4,3,5,2,1,0) (new node i = old node p[i], graphs.permute). Then
```
G1  = 011011 / 101011 / 110101 / 111000 / 011000 / 000110
G2p = 011011 / 101001 / 110111 / 111000 / 011000 / 000110
```
differ in exactly two entries: the edge 1->4 of G1 is replaced by 2->4. In G1 nodes 1 and 2 are
bidirected (1<->2) and have the same in-neighbours among the others ({0,3,4}). Hence
`N(G2p) = N(G1) + e_4 (e_1 - e_2)^T`, and `x = e_1 - e_2` satisfies `x^T N(G1) = 0`, `x^T 1 = 0`.
Exact check (s2_g1g2.py): `beta_data(G1)` and `beta_data(G2p)` agree entry for entry -- `a` agrees
as a VECTOR in the same labelling, not only as a multiset.
