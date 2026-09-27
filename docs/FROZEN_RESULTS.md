# Frozen results

This page records the submission-aligned numerical evidence for this study. Values below are frozen from the audited V3/V3.1 result packages; they should not be silently replaced by later exploratory reruns.

## Representation benchmark

| Representation | Features | Random R² mean ± SD | Scaffold R² mean ± SD | Random RMSE | Scaffold RMSE |
|---|---:|---:|---:|---:|---:|
| Graph19 | 19 | 0.392 ± 0.038 | 0.211 ± 0.136 | 1.068 | 1.197 |
| RDKit2D | 196 | 0.657 ± 0.017 | 0.518 ± 0.100 | 0.802 | 0.933 |
| ECFP4 | 2,048 | 0.735 ± 0.014 | 0.615 ± 0.062 | 0.705 | 0.837 |

Matched ECFP4 gains over Graph19 were **+0.343 R²** (95% CI 0.309–0.377) under random splitting and **+0.404 R²** (95% CI 0.310–0.498) under scaffold splitting; all five matched partitions favored ECFP4.

## Exact Graph19 representation degeneracy

- 10,056 molecules yield **7,382 unique Graph19 descriptor vectors**.
- **1,507 exact collision groups** contain **4,181 molecules (41.6%)**.
- **1,285 groups (85.3%)** contain more than one molecular formula.
- **117 groups** span more than **2 pIC50 units**.
- **108 molecular pairs** combine an identical Graph19 vector, ECFP4 Tanimoto < 0.5, and |ΔpIC50| > 2.
- The whole-dataset empirical RMSE floor induced by exact Graph19 collisions is **0.332 pIC50 units**; this is an in-sample ambiguity diagnostic, not a generalization bound.

## RF feature-sampling sensitivity

Changing the matched Random Forest benchmark from `sqrt(p)` feature sampling to a common **50% feature fraction** preserved the representation ordering:

| Representation | Random R² | Scaffold R² |
|---|---:|---:|
| Graph19 | 0.393 | 0.214 |
| RDKit2D | 0.670 | 0.530 |
| ECFP4 | 0.734 | 0.610 |

Selected full-feature spot checks preserved the same ordering.

## Descriptor-family ablation

| Graph19 family | Features | Random R² | Scaffold R² |
|---|---:|---:|---:|
| Degree only | 10 | 0.265 | 0.065 |
| Distance only | 7 | 0.359 | 0.189 |
| Spectral only | 2 | 0.327 | 0.135 |
| Degree + distance | 17 | 0.385 | 0.204 |
| Full Graph19 | 19 | 0.392 | 0.211 |

## Uncertainty under shift

Nominal 90% split-conformal coverage is close to nominal under the random split but drops under scaffold shift. Scaffold-fold coverage is reported as an empirical stress test only; ordinary split-conformal exchangeability guarantees do not apply under that distribution shift.

## Integrity checkpoints

- V3 final cohort: **10,056 molecules**.
- Graph descriptors: **19/19 regenerated**, maximum discrepancy ≈ `5e-8` relative to the archived matrix.
- Murcko scaffolds: **3,579 unique scaffolds**.
- V3.1 return-package manifest: **59/59 recorded outputs matched SHA-256 hashes**.

Machine-readable sources are placed under `results/frozen/` and `results/v31/`.
