<p align="center">
  <img src="docs/assets/repository-banner.svg" alt="EGFR Graph QSAR Reproducibility" width="100%" />
</p>

<h1 align="center">EGFR Graph QSAR Reproducibility</h1>

<p align="center">
  <strong>Reproducibility repository for</strong><br/>
  <strong>Mechanisms of Representation Degeneracy and Generalization Limits of Classical Molecular Graph Descriptors in EGFR QSAR</strong>
</p>

<p align="center">
  <a href="docs/FROZEN_RESULTS.md"><img src="https://img.shields.io/badge/reproducibility-V3.3%20frozen-2ea44f.svg" alt="V3.3 frozen results"/></a>
  <a href="environment.yml"><img src="https://img.shields.io/badge/Python-3.11-3776AB.svg?logo=python&logoColor=white" alt="Python 3.11"/></a>
  <a href="docs/DATA_PROVENANCE.md"><img src="https://img.shields.io/badge/data-ChEMBL%20provenance-5b8c85.svg" alt="Data provenance"/></a>
  <a href="CITATION.cff"><img src="https://img.shields.io/badge/citation-CITATION.cff-blue.svg" alt="Citation metadata"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/code%20license-MIT-lightgrey.svg" alt="MIT license"/></a>
</p>

## Why this study?

Classical graph invariants provide an unusually compact and mathematically interpretable encoding of molecular topology. Their compression is useful, but it can also erase chemical distinctions before a predictive model ever sees them.

This study asks two linked questions on a frozen **10,056-compound human EGFR pIC50 cohort**:

1. how much predictive signal remains accessible after compressing each molecule to **19 classical graph invariants (Graph19)**; and
2. **which molecular information is lost, at which abstraction layer, and how that loss relates to generalization and reliability**.

The study is deliberately not framed as a contest to make Graph19 outperform contemporary fingerprints. Graph19 is treated as an interpretable compression whose limitations can be diagnosed explicitly.

## Frozen study design

| Component | Frozen design |
|---|---|
| Target | Human EGFR (`CHEMBL203`) |
| Endpoint | IC50 transformed to pIC50 inhibitory potency |
| Source records | 17,705 |
| Curated activity rows | 17,412 |
| Final unique compounds | **10,056** |
| Compact graph representation | Graph19 - 19 classical invariants |
| Descriptor baseline | RDKit2D - 196 retained descriptors |
| Fingerprint baseline | ECFP4 - Morgan radius 2, 2,048-bit binary, chirality disabled |
| Random validation | 5 fixed 80:20 partitions |
| Scaffold validation | 5 Bemis-Murcko GroupKFold partitions |
| Principal learner | Random Forest |
| Second learner sensitivity | ExtraTreesRegressor |
| Reliability analyses | chemical similarity, leverage, conformal intervals, collision-conditioned held-out error |

## Representation benchmark

| Representation | Features | Random R² mean ± SD | Scaffold R² mean ± SD | Random RMSE | Scaffold RMSE |
|---|---:|---:|---:|---:|---:|
| Graph19 | 19 | 0.392 ± 0.038 | 0.211 ± 0.136 | 1.068 | 1.197 |
| RDKit2D | 196 | 0.657 ± 0.017 | 0.518 ± 0.100 | 0.802 | 0.933 |
| ECFP4 | 2,048 | **0.735 ± 0.014** | **0.615 ± 0.062** | **0.705** | **0.837** |

Graph19 therefore retains measurable EGFR pIC50 predictive signal at extreme dimensional compression, but richer chemical representations recover substantially more signal, especially under scaffold shift.

## V3.3 mechanism-resolved representation degeneracy

At the operational **8-decimal Graph19 encoding**:

- **7,382** unique Graph19 vectors;
- **1,507** collision groups;
- **4,181 molecules (41.58%)** in collision groups;
- maximum within-class pIC50 range **5.17185**.

The new V3.3 graph-isomorphism/canonical-topology audit localizes where those collisions arise. Every 8-decimal collision group contained only **one hydrogen-suppressed unlabeled simple-graph topology**. No observed collision group required two genuinely different unlabeled simple graphs to share the same Graph19 vector.

Primary earliest information-loss mechanisms were:

| Mechanism | Collision groups | Molecules |
|---|---:|---:|
| Atom/element-label loss | **1,278** | **3,691** |
| Bond-order/aromaticity loss | 104 | 222 |
| Stereo/other full-identity loss | 125 | 268 |
| Different unlabeled simple topology | **0** | **0** |

This supports the interpretation that the dominant empirical degeneracy in this EGFR cohort occurs at the **chemical structure -> hydrogen-suppressed unlabeled simple graph** abstraction step, rather than from observed non-isomorphic unlabeled simple graphs colliding only after the 19-index compression.

## Precision sensitivity

Rounded equality is treated as an operational encoding, not as mathematical real-number equality.

| Encoding | Collision groups | Molecules in collisions | Cohort fraction |
|---|---:|---:|---:|
| 6 decimals | 1,507 | 4,181 | 41.58% |
| 8 decimals | 1,507 | 4,181 | 41.58% |
| 10 decimals | 1,507 | 4,180 | 41.57% |
| 12 decimals | 1,507 | 4,169 | 41.46% |
| exact binary64 | 1,408 | 3,719 | 36.98% |

The empirical collision burden is therefore highly stable across 6-12 decimal encodings, while exact binary64 equality is somewhat more discriminating but still leaves a substantial collision burden.

## Whole-vector ECFP4 comparison

This analysis is distinct from within-molecule Morgan hash-folding collisions. It asks whether **different molecules share the same complete 2,048-bit vector**.

- **9,661** unique ECFP4 vectors;
- **307** identical full-vector collision groups;
- **702 molecules (6.98%)** in those groups;
- **20 groups** span more than 2 pIC50 units;
- maximum pIC50 range **3.0**.

ECFP4 is therefore not described as collision-free; it is simply much more discriminating than Graph19 on this cohort.

## Collision-conditioned held-out prediction error

Graph19 ambiguity is not uniformly associated with larger prediction error. The important distinction is whether the held-out encoded class is represented in training and whether that class contains strong response disagreement.

Under scaffold evaluation, Graph19 MAE was approximately:

- **0.636** when the matching Graph19 vector was seen in training;
- **1.018** when the collision-class vector was not seen in training;
- **1.186** for high-response-disagreement collision classes.

For the high-disagreement scaffold subgroup, RDKit2D and ECFP4 reduced MAE to approximately **0.938** and **0.863**, respectively. These are associative diagnostics, not causal claims.

## Second-learner representation sensitivity

ExtraTreesRegressor preserved the same representation ordering:

| Split | Graph19 R² | RDKit2D R² | ECFP4 R² |
|---|---:|---:|---:|
| Random | 0.366 ± 0.040 | 0.692 ± 0.020 | **0.729 ± 0.018** |
| Scaffold | 0.163 ± 0.165 | 0.556 ± 0.099 | **0.602 ± 0.070** |

Thus the qualitative ordering **ECFP4 > RDKit2D > Graph19** is not specific to the Random Forest benchmark.

## Data and code availability

The public repository contains the reproducibility code, V3.3 strengthening code, environment locks, data provenance, checksums, frozen result summaries, and the final modeled EGFR cohort as a compressed CSV archive under `data/egfr_qspr_final_v33.zip`.

ChEMBL-derived data retain their source licensing terms and attribution requirements. Original project code is released under the MIT License. The missing historical ChEMBL release identifier is disclosed rather than inferred.

## Reproducing the study

```bash
conda env create -f environment.yml
conda activate egfr-graph-qsar
```

For the original V3 pipeline, use `RUN_EGFR_QSAR_REPRODUCIBILITY.cmd`. The V3.3 incremental strengthening code and its pinned environment are under `v33/`.

The V3.3 code reuses the frozen cohort and validation design rather than substituting a newer ChEMBL release.

## Result-to-source map

| Evidence | Repository source |
|---|---|
| Random/scaffold representation benchmark | `results/frozen/representation_benchmark_summary.csv` |
| Descriptor-family ablation | `results/frozen/descriptor_family_ablation_summary.csv` |
| Chemical-space shift | `results/frozen/chemical_space_similarity_bins.csv` |
| Conformal uncertainty | `results/frozen/conformal_metrics.csv` |
| V3.1 paired representation effects | `results/v31/paired_representation_effect_sizes.csv` |
| V3.3 mechanism decomposition | `results/v33/action_a_mechanism_summary.csv` |
| V3.3 precision sensitivity | `results/v33/action_b_precision_sensitivity_summary.csv` |
| V3.3 Graph19 vs ECFP4 collision comparison | `results/v33/action_c_graph19_vs_ecfp4_collision_comparison.csv` |
| V3.3 collision-conditioned error | `results/v33/action_d_summary_across_partitions.csv` |
| V3.3 ExtraTrees sensitivity | `results/v33/action_e_extratrees_summary.csv` |
| Complete V3.3 scientific freeze | `results/v33/V33_SCIENCE_FREEZE_SUMMARY.json` |

## Scientific scope and limitations

- pIC50 is an inhibitory-potency endpoint, not a thermodynamic binding-affinity constant.
- Graph19 intentionally omits atom labels, bond orders, stereochemistry, 3D geometry, protein state and assay context.
- Eight-decimal Graph19 equality is an operational encoding; exact binary64 results are reported separately.
- The observed absence of non-isomorphic unlabeled-topology collisions is cohort-specific and is not a theorem that the 19-index vector is globally injective.
- ECFP4 is not collision-free.
- Representation comparisons are learner- and protocol-mediated, not information-theoretic proofs.
- Scaffold conformal coverage is an empirical distribution-shift stress test, not an exchangeability-guaranteed coverage result.
- The historical ChEMBL release identifier was not preserved.

## Authors

- **Sunilgar L. Gusai** - corresponding author, Marwadi University  
  ORCID: https://orcid.org/0009-0004-0739-4812
- **Manoharsinh R. Jadeja** - Marwadi University

## Citation

Citation metadata are provided in `CITATION.cff`. Publication DOI/journal metadata will be added when available.
