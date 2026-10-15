# Lean formalisation: Remark F, its spectral consequence, Lemma R (forward direction)

Source mathematics: notes/lemmas.md section 7 (Remark F), notes/lemma_M_refuted.md section 3.1 (Lemma R),
ANOMALIES.md entry "2026-10-14 Lemma R". Conventions: CLAUDE.md section 2 and notes/lemmas.md section 0.
Nothing in the repository was modified; everything lives under this directory.

## Toolchain and Mathlib (what worked)

- elan 4.2.4 (pre-installed). Installed `leanprover/lean4:v4.35.0` (the toolchain of Mathlib's `v4.35.0` tag
  and of Mathlib master on 2026-10-15), set as elan default.
- Lake project `CtlnCert/` requires Mathlib at tag `v4.35.0` (commit 56a005855ea17e0b0958f73ec801ae140b0a1679).
- `lake update` fetched Mathlib and its dependencies and its post-update hook downloaded the prebuilt cache from
  https://cache.mathlib.org: 9031 .olean archives, about 95 s. A following `lake exe cache get` reported
  "No files to download". Mathlib was NOT built from source. `.lake` is about 7.3 GB.

## Files

| file | content |
|---|---|
| `CtlnCert/lakefile.toml`, `lean-toolchain`, `lake-manifest.json` | project, pinned toolchain and Mathlib |
| `CtlnCert/CtlnCert/Basic.lean` | `J`, `K N beta t = 1 + beta N + t J`, `Eig` (eigenvector over C), `cplx`, `PosStable`, `HurwitzStable`; `eig_iff_isRoot_charpoly`, `posStable_iff_charpoly` |
| `CtlnCert/CtlnCert/RayIdentity.lean` | Remark F: `affine_identity`; charpoly scaling `charpoly_affine_eval`, `charpoly_affine`, `charpoly_ray`; `eig_affine`, `cplx_affine`, `eig_ray_decomp`; `posStable_ray`, `posStable_convex` |
| `CtlnCert/CtlnCert/IVTStep.lean` | `fR`, `ivt_step` (Mathlib IVT), `curve_point_of_root`; explicit root `lamStar`, `explicit_root_abstract` (abstract sqrt), `explicit_root` (Real.sqrt), `root_unique` |
| `CtlnCert/CtlnCert/Bridge.lean` | `W` exactly as CLAUDE.md section 2 (reads `A j i`), `Nmat A = J - 1 - Aᵀ`; `closed_form`, `IminusW_eq`, `legal_iff`, `hurwitz_iff_posStable` |
| `CtlnCert/CtlnCert/LemmaR.lean` | `lemmaR_forward` (K level), `lemmaR_forward_ctln` (eps, delta level), `curve_shift_legal` |
| `CtlnCert/CtlnCert.lean` | root module importing all of the above |
| `CtlnCert/Audit.lean` | `#print axioms` for all 26 theorems |
| `CtlnCert/build.log`, `CtlnCert/audit.log` | output of the last clean build and of the axiom audit |

## Commands (re-run)

```
export PATH=$HOME/.elan/bin:$PATH
cd <this directory>/CtlnCert
lake exe cache get          # no-op once the cache is present
lake build                  # builds CtlnCert (default target); "Build completed successfully"
lake env lean Audit.lean    # every theorem: [propext, Classical.choice, Quot.sound] (affine_identity: [propext, Quot.sound]); no sorryAx
grep -rn "sorry\|admit\|^axiom" CtlnCert CtlnCert.lean   # empty
```

The only build warnings are Lake's advisory that Mathlib "is designed for use with the module system".

## Statements proved (no sorry)

- `affine_identity` (any commutative ring, any l): `K N (l*beta) (l*t) = (1 - l) • 1 + l • K N beta t`.
- `charpoly_affine_eval` (any field, l ≠ 0): `eval x (charpoly (c • 1 + l • M)) = l^m * eval ((x - c)/l) (charpoly M)`;
  `charpoly_affine` the same as a polynomial identity (infinite field); `charpoly_ray` its specialisation to K.
- `eig_ray_decomp` (l ≠ 0): every complex eigenvalue mu of K(l beta, l t) is `(1 - l) + l nu` for an eigenvalue nu of K(beta, t).
- `posStable_ray`: K(beta, t) positive stable and 0 < l ≤ 1 implies K(l beta, l t) positive stable.
- `ivt_step`: t > 0, t(beta - 2) > 1 implies some l in (2/beta, 1) with `l t (l beta - 2) - 1 = 0`;
  `curve_point_of_root`: then l t > 0, l beta > 2, (l t)(l beta - 2) = 1, l beta = 2 + 1/(l t).
- `explicit_root`: `lamStar = (t + sqrt(t^2 + beta t))/(beta t)` is such a root and lies in (2/beta, 1);
  `explicit_root_abstract` the same with any s > 0, s^2 = t^2 + beta t (no sqrt, no IVT);
  `root_unique`: f has only one positive root.
- `closed_form`: `1 - W = (-delta) • 1 + (1 + delta) • J - (eps + delta) • Aᵀ` (with the transpose).
- `IminusW_eq`: `1 - W = eps • K (J - 1 - Aᵀ) (1 + delta/eps) ((1 - eps)/eps)` (eps ≠ 0).
- `legal_iff` (0 < eps < 1): `delta(1-eps)/eps > 1 ↔ t (beta - 2) > 1`.
- `hurwitz_iff_posStable` (eps > 0): `-1 + W` Hurwitz stable ↔ `K` positive stable.
- `lemmaR_forward`: K(beta, t) positive stable, t > 0, t(beta - 2) > 1 implies some t' in (0, t) with
  2 + 1/t' < beta and K(2 + 1/t', t') positive stable.
- `lemmaR_forward_ctln`: if `-1 + W A eps delta` is Hurwitz stable at a legal point (0 < eps < 1, q > 1), then
  `-1 + W A eps' delta'` is Hurwitz stable for some eps < eps' < 1, delta' > 0 with q(eps', delta') = 1 exactly
  (eps' = 1/(1 + t'), delta' = 1/t').

## Not formalised

- The converse (<=) of Lemma R needs openness of positive stability (continuity of eigenvalues in the
  entries). Mathlib has no ready-made statement for that, so it is not formalised. Only the elementary half,
  `curve_shift_legal` (the point (2 + 1/t' + eta, t') is legal for eta > 0), is proved.
- Conditions 1 and 2 of the stable-fixed-point certificate (positivity, off-support) are not touched here.
