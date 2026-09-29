# Frozen results — V3.3 submission science

This page records the submission-aligned numerical evidence frozen for V3.3. Earlier V3/V3.1 result tables are retained for provenance, but the statements below govern the final manuscript.

## Principal representation benchmark

| Representation | Features | Random R² mean ± SD | Scaffold R² mean ± SD | Random RMSE | Scaffold RMSE |
|---|---:|---:|---:|---:|---:|
| Graph19 | 19 | 0.392 ± 0.038 | 0.211 ± 0.136 | 1.068 | 1.197 |
| RDKit2D | 196 | 0.657 ± 0.017 | 0.518 ± 0.100 | 0.802 | 0.933 |
| ECFP4 | 2,048 | 0.735 ± 0.014 | 0.615 ± 0.062 | 0.705 | 0.837 |

Matched ECFP4 gains over Graph19 were +0.343 R² (95% CI 0.309–0.377) under random splitting and +0.404 R² (95% CI 0.310–0.498) under scaffold splitting.

## Graph19 representation degeneracy at the operational eight-decimal encoding

- cohort: **10,056 molecules**;
- unique Graph19 vectors: **7,382**;
- collision groups: **1,507**;
- molecules in collision groups: **4,181 (41.58%)**;
- maximum within-group pIC50 range: **5.17185**.

The V3.3 canonical-topology/isomorphism audit found that every eight-decimal collision group contained only one hydrogen-suppressed unlabeled simple-graph topology. No observed group required two genuinely different unlabeled simple graphs to share the same Graph19 vector.

Primary earliest information-loss mechanisms:

| Mechanism | Groups | Molecules |
|---|---:|---:|
| Atom/element-label loss | 1,278 | 3,691 |
| Bond-order/aromaticity loss | 104 | 222 |
| Stereo/other full-identity loss | 125 | 268 |
| Different unlabeled simple topology | 0 | 0 |

This is an empirical cohort result, not a proof that the 19-index vector is globally injective on all unlabeled molecular graphs.

## Precision sensitivity

Rounded equality is an operational encoding rather than mathematical real-number equality.

| Encoding | Collision groups | Molecules in collision groups | Cohort fraction |
|---|---:|---:|---:|
| 6 decimals | 1,507 | 4,181 | 41.58% |
| 8 decimals | 1,507 | 4,181 | 41.58% |
| 10 decimals | 1,507 | 4,180 | 41.57% |
| 12 decimals | 1,507 | 4,169 | 41.46% |
| exact binary64 | 1,408 | 3,719 | 36.98% |

Thus the empirical conclusion is highly stable across 6–12 decimal encodings, while exact binary64 equality is somewhat more discriminating but still leaves substantial degeneracy.

## Whole-vector ECFP4 collisions

This analysis concerns different molecules sharing the same complete 2,048-bit fingerprint and is distinct from within-molecule Morgan folding collisions.

- unique ECFP4 vectors: **9,661**;
- identical complete-vector collision groups: **307**;
- molecules in those groups: **702 (6.98%)**;
- groups with pIC50 range >2: **20**;
- maximum pIC50 range: **3.0**.

ECFP4 is therefore not collision-free; it is substantially more discriminating than Graph19 on this cohort.

## Collision-conditioned held-out prediction error

Graph19 ambiguity is not uniformly associated with high error. Under scaffold evaluation, mean Graph19 MAE across folds was approximately **0.636** when the matching collision vector was present in training and **1.018** when the collision-class vector was unseen. The post-hoc high-response-disagreement collision subgroup had Graph19 MAE approximately **1.186**, compared with **0.938** for RDKit2D and **0.863** for ECFP4. These are associative diagnostics, not causal claims.

## Second-learner sensitivity

ExtraTreesRegressor preserved the representation ordering:

| Split | Graph19 R² | RDKit2D R² | ECFP4 R² |
|---|---:|---:|---:|
| Random | 0.366 ± 0.040 | 0.692 ± 0.020 | 0.729 ± 0.018 |
| Scaffold | 0.163 ± 0.165 | 0.556 ± 0.099 | 0.602 ± 0.070 |

The qualitative ordering ECFP4 > RDKit2D > Graph19 is therefore not specific to the Random Forest benchmark.

## Earlier frozen strengthening retained in V3.3

- Graph19 descriptor-family ablation;
- held-out permutation interpretation;
- nearest-training ECFP4/Tanimoto chemical-space diagnostics;
- corrected Williams leverage applicability-domain analysis;
- split-conformal prediction;
- RF feature-sampling sensitivity;
- Y-randomization and historical robustness analyses.

## Integrity checkpoints

- authoritative V3.2 package was verified against its internal SHA-256 manifest: **300/300 files matched**;
- V3.3 Graph19 regeneration covered **10,056/10,056 molecules** and passed parity against the frozen descriptor cache;
- maximum V3.3 regeneration discrepancy was approximately **5.68×10⁻14**;
- V3.3 machine-readable result artifacts were verified against their recorded hashes before manuscript rebuild.

V3.3 machine-readable sources are under `results/v33/`, with source-freeze and environment metadata under `v33/`.
