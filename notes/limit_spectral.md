# The eps -> 0 limit as a spectral condition on the compressed adjacency

Status: orchestrator's derivation, verified numerically on Geneson's graph
(beta_c = 2.16635502785 agrees to 12 digits with the exact sextic root from
fastregion/limit and from symbolic sympy). A rigorous proof is requested from
the prover (Lemma I). Do not cite as a theorem until notes/lemmas.md has it.

Setting (fastregion.py): full support, m nodes, N = receiver-first nonedge
matrix = J - I - A^T, C = I + beta N, K = C + t J, t = (1-eps)/eps,
beta = 1 + delta/eps, q = (beta-1)(1-eps) -> beta - 1 as eps -> 0.
Let P = I - J/m be the orthogonal projector onto 1^perp, and for a matrix X
write X|perp for the compression of X to 1^perp (an (m-1)x(m-1) matrix in any
orthonormal basis of 1^perp). Classical identity: det(lam I - X|perp) =
(1/m) 1^T adj(lam I - X) 1.

Claim 1 (limit eigenvalues). As t -> infinity the eigenvalues of K are one
eigenvalue ~ m t and m-1 eigenvalues converging to the eigenvalues of C|perp.
(Secular equation: det(lam I - K) = det(lam I - C) - t 1^T adj(lam I - C) 1.)

Claim 2 (limit stability). Since J|perp = 0, C|perp = (1 - beta) I - beta A^T|perp,
so the limit eigenvalues are 1 - beta - beta * mu with mu ranging over the
eigenvalues of A|perp (same spectrum as A^T|perp). Hence, with
    gamma := 1 - 1/beta = q/(q+1)  in (1/2, 1)   for legal beta > 2,
K is limit-stable at beta  <=>  every eigenvalue mu of A|perp has Re mu < -gamma
                           <=>  spectral abscissa alpha(A|perp) < -gamma.
The limit-stable beta-range is therefore the interval 2 < beta < beta_c with
    beta_c = 1/(1 + alpha(A|perp)),   q_c = -alpha/(1 + alpha),
nonempty iff alpha(A|perp) < -1/2. For the complete graph A|perp = -I, alpha = -1,
stable for all beta. Geneson: alpha = -0.538395..., beta_c = 2.166355...

Claim 3 (limit positivity). x > 0 in the limit <=> adj(C) 1 is one-signed
<=> C^{-1} 1 is one-signed. Solving C u = 1 with C = (1-beta) I + beta J - beta A^T
gives u proportional to w := (A^T + gamma I)^{-1} 1, i.e. w is the solution of
    gamma w_i + sum_{j -> i} w_j = 1   for every node i.
If w is one-signed it is positive (the left side of each equation would be
negative otherwise), so: limit-positive at beta <=> w(gamma) > 0 entrywise,
where gamma = 1 - 1/beta. (Nondegeneracy: C invertible and 1 - beta 1^T u != 0.)

Consequences.
* trace(A|perp) = -|E|/m, so the mean of the m-1 eigenvalues is -|E|/(m(m-1)) =
  -(edge density); alpha < -1/2 forces edge density > 1/2: limit
  counterexamples are dense (more than half of the m(m-1) ordered pairs are
  edges). Geneson: 19/30.
* Equivalently in terms of the complement adjacency Abar = J - I - A:
  limit-stable at beta <=> all eigenvalues of Abar|perp have Re > -1/beta.
  A 2-node bidirected nonedge {i,j} alone would give eigenvalue -1 on e_i - e_j;
  in a larger graph the compression mixes, which is why Geneson's {4,5} pair
  (no edge either way) is survivable.
* For n = 7 this suggests enumerating complements with few edges rather than all
  880,471,142 digraphs, once the n = 6 limit list shows the pattern.
