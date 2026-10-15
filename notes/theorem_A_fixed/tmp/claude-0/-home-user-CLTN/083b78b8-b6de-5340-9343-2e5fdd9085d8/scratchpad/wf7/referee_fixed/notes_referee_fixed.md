# Referee report: prover_fixed (one fixed legal point, infinitely many non-clique SFP full supports)

All scripts in this directory; exact arithmetic only (Fraction / sympy Rational); ctln imported read-only. No repo file modified.
Every graph with n >= 8 below is a CANDIDATE counterexample pending human hand-verification (CLAUDE.md section 3).

## Verdict: the main claim HOLDS

R = &F^U|m^QTq? (decoded by own decoder `reflib.d6_decode`; rows match the prover's), S = (1,1,0,1,1,1,1),
R_k = R[S, all; k] (`reflib.build_Rk`), (eps, delta) = (1/21, 11/210), q = 22/21 (legal).

PROVED (independent route, `proof_Rk.py` -> `proof_Rk_eps1_21.json`). Not via the prover's Lemma S / Lemma P / D-picture.
Use the equitable partition {0},...,{6},{class} of M = -I + W. Class vectors with zero sum form an M-invariant subspace
with eigenvalue -eps, so charpoly(M) = chi_Q(lam) (lam + eps)^(k-1), where Q is the 8x8 quotient with k symbolic.
- The identity was checked exactly against ctln.charpoly on the explicit graphs, k = 1..6 and k = 30.
- All 8 Hurwitz minors of chi_Q are polynomials in k with all coefficients positive (Sturm: no roots on [1, oo)).
- Positivity: x = P y with y_i = c_i / (3433400 k + 8122137), c_i > 0.
- k = 0 is R itself (ctln).
Hence R_k is a stable positive fixed point (full support, non-clique) for every k >= 0. This agrees with Theorem A,
including the positivity denominator, which matches the prover's Delta_5 constant 13(3433400k+8122137).

CERTIFIED (ctln.analyze_support at the same (1/21, 11/210): legal, x > 0, stable True with all Hurwitz minors positive,
not a clique): k = 0,1,2,3,4,5,6,10,20 (`cert_Rk_small.json`) and k = 41, 60, 100 (n = 48, 67, 107; `cert_Rk_large.json`).
k = 41 is past F_k's death (F_k unstable for k >= 41 at this point).

## Secondary claims checked

- Theorem A+ (beta <= 43/20, every legal point, every k): CHECKED, not re-proved symbolically in (a, k).
  - `region_Aplus.py` -> `region_Aplus.json` (per-point all-k quotient proofs in `pts/`) and `region_k0.py` -> `region_k0.json`.
  - At 28 exact legal points (corners beta = 43/20 with t = 20/3 + 1/1000 and t = 10^4; beta = 21/10, 41/20, 201/100;
    8 random points), the all-k proof succeeds and k = 0 is SFP.
  - Control: at beta = 11/5, t = 5.01, k = 0..3 fail and k = 4, 5 pass, matching Theorem B's k0 = 4 there.
- Theorem SJ (self-join R^{vk}, beta <= 21/10): PROVED independently at the main point (`proof_SJ.py` -> `proof_SJ_eps1_21.json`).
  - Kronecker block-diagonalisation gives charpoly(M) = chi(B0 + (k-1)(eps-1)J) chi(B0 + (1-eps)J)^(k-1).
    This was checked against ctln for k = 1..4.
  - All minors of the k-factor are positive on [1, oo), and the fixed factor is Hurwitz.
  - Positivity follows from Sherman-Morrison, since (I - W_R)^{-1} 1 > 0.
  - Same result at 3 further points with beta <= 21/10. At beta = 43/20 the fixed factor fails (as the prover states).
    See `sj_points.json`.
  - ctln certifies k = 2, 3, 5 (n = 14, 21, 35).
- (1/2, 2) instance of Theorem Lambda: own quotient method (`check_lambda_instance.py` -> `check_lambda_instance.json`).
  - charpoly(Q) is affine in k (checked at k = 0, 1, 2) and has the factor (lam + 1/2)^59. The remaining part has degree 7.
  - Exact threshold k* = 1581: k = 1581 passes with all minors positive, and k = 1580 fails (minor 5 negative).
  - Positivity holds for every k >= 1.
  - This matches the prover's k0 = 1581.
- Y family, Theorem FP and B': own multi-class quotient checker (`quotient_family.py`).
  - At (1/2, 9/8) the threshold is k* = 31, matching the prover.
  - At beta = 24/5 the thresholds are k* = 1413 (t = 5/14 + 1/1000) and 363 (t = 1000). Both are <= 1764, consistent with Theorem B.
  - At beta = 49/10 (t = 3), k* = 1325.
  - At beta = 5 (> beta_Y) the family is not eventually stable, consistent with the ceiling 4.9607.
  - Files: `qf_Y_*.json`.

## Not refereed here (OPEN for this referee)
- Theorem Lambda for EVERY legal point is conditional on the cited Lemma (S) (Rouche rounding) of the parallel Theorem U.
  Its Step-3 extension to arbitrary fixed r was read but not independently re-proved. Only the beta = 5 instance was
  re-verified.
- The general Rouche proof of FP was read and found correct (root sum a(m+P), and one real escaping root). Only instances
  were machine-checked.
