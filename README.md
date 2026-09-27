# EGFR Graph QSAR Reproducibility

**Reproducibility repository for:**  
**Classical Molecular Graph Descriptors versus Contemporary Representations for EGFR pIC50 Prediction under Scaffold and Chemical-Space Shift**

[![License: MIT](https://img.shields.io/badge/code%20license-MIT-blue.svg)](LICENSE)
[![Data: ChEMBL](https://img.shields.io/badge/data-ChEMBL-5b8c85.svg)](DATA_LICENSE_NOTICE.md)
[![Status](https://img.shields.io/badge/status-pre--submission-orange.svg)](#release-status)

## Why this study?

The paper treats **19 classical mathematical molecular-graph invariants (Graph19)** as an intentionally compressed representation rather than claiming that they are an accuracy-optimized replacement for contemporary cheminformatics features. The central question is:

> **What EGFR pIC50 predictive signal survives compression to only 19 classical graph invariants, and where does that compression fail?**

The same frozen 10,056-compound EGFR cohort is evaluated using Graph19, 196 RDKit2D descriptors and 2,048-bit ECFP4 fingerprints under matched Random Forest principles, random splits, Bemis-Murcko scaffold folds, chemical-space diagnostics and uncertainty analysis.

## Key findings

| Result | Graph19 | RDKit2D | ECFP4 |
|---|---:|---:|---:|
| Features | 19 | 196 | 2,048 |
| Random-split mean R2 | 0.392 | 0.657 | **0.735** |
| Scaffold mean R2 | 0.211 | 0.518 | **0.615** |
| Random RMSE | 1.068 | 0.802 | **0.705** |
| Scaffold RMSE | 1.197 | 0.933 | **0.837** |

Matched ECFP4 gains over Graph19 were **+0.343 R2** (95% CI 0.309-0.377) under random splitting and **+0.404 R2** (0.310-0.498) under scaffold splitting; the direction was positive in every matched split/fold.

### Representation ambiguity is measurable

- 10,056 molecules collapse to **7,382 distinct Graph19 vectors**.
- **1,507 exact Graph19 descriptor-vector collision groups** contain **4,181 molecules (41.6%)**.
- **1,285 (85.3%)** collision groups contain more than one molecular formula.
- 117 groups span more than **2 pIC50 units**.
- **108 molecular pairs** have identical Graph19 vectors, ECFP4 Tanimoto <0.5, and |delta pIC50| >2.
- The representation ranking remains **ECFP4 > RDKit2D > Graph19** when RF feature sampling is changed from `sqrt(p)` to a common 50% feature fraction.

## Frozen data lineage

```text
17,705 raw IC50 records
      ↓
17,412 curated activity records
      ↓
10,056 unique compounds (InChIKey deduplication, median pIC50)
```

Of the 17,412 curated activity labels, 17,156 use ChEMBL `pchembl_value`; 256 use the compatible unit-derived fallback. All 19 graph descriptors were independently regenerated and matched the archived matrix to within `5e-8`. The study recovered **3,579 unique Bemis-Murcko scaffolds**.

See `DATA_LICENSE_NOTICE.md` for ChEMBL attribution/licensing and `INPUT_MANIFEST.json` for frozen input hashes.

## Repository map

```text
.
├── README.md
├── LICENSE                         # MIT - original study code
├── DATA_LICENSE_NOTICE.md          # ChEMBL-derived data terms
├── CITATION.cff
├── requirements.txt                # supported dependency ranges
├── requirements-lock.txt           # exact definitive environment
├── MASTER_RUN_PAPER003_V3.cmd      # one-click full Windows run
├── SMOKE_TEST_PAPER003_V3.cmd
├── config.json
├── INPUT_MANIFEST.json
├── data/README.md                  # frozen-data access + archive route
├── scripts/                        # V3 assurance/benchmark pipeline
├── v31/                            # representation-ambiguity + RF-sensitivity strengthening
├── results_key/                    # principal V3 machine-readable outputs
├── results_v31/                    # compact V3.1 machine-readable outputs
├── figures/                        # selected publication figures used by the README
├── docs/                           # design freeze, QA and review audit
└── .github/workflows/qa.yml        # lightweight source/manifest QA
```

## Reproduce the definitive V3 benchmark on Windows

1. Clone or download this repository.
2. Download/extract the frozen data payload from the Zenodo submission archive into the paths documented by `INPUT_MANIFEST.json`.
3. Keep the frozen inputs unchanged.
4. Double-click `MASTER_RUN_PAPER003_V3.cmd`.
5. The script bootstraps a virtual environment, validates input hashes, executes the phased pipeline, and writes a return package under `results/`.

A reduced software-path check is available as `SMOKE_TEST_PAPER003_V3.cmd`. The smoke test verifies software plumbing only; its reduced results are **not** scientific results.

## Reproduce the V3.1 strengthening analysis

`v31/v31_strengthening.py` implements:

1. exact Graph19 descriptor-vector collision analysis;
2. within-collision ECFP4/pIC50 discordance analysis and representative structures;
3. paired representation effect sizes across the matched V3 splits/folds;
4. RF feature-sampling sensitivity at a common 50% feature fraction;
5. selected full-feature RF spot checks.

The immutable Zenodo reproducibility archive retains the frozen V3 return archive consumed by the strengthening script so V3.1 can be independently audited without reconstructing an intermediate archive by hand.

## Environment

Definitive V3 environment:

- Python 3.11.9
- RDKit 2026.03.6
- scikit-learn 1.8.0
- NumPy 2.4.6
- SciPy 1.17.1
- pandas 2.3.3
- matplotlib 3.11.2

Use `requirements-lock.txt` for the exact frozen package versions and `requirements.txt` for compatible ranges.

## Result-to-source map

| Manuscript result | Machine-readable source |
|---|---|
| Graph19/RDKit2D/ECFP4 random + scaffold benchmark | `results_key/representation_benchmark_metrics.csv` |
| Benchmark summary | `results_key/representation_benchmark_summary.csv` |
| Paired representation effects | `results_v31/paired_representation_effect_sizes.csv` |
| Graph19 family ablation | `results_key/descriptor_family_ablation_summary.csv` |
| Individual/grouped permutation importance | `results_key/interpretability_*` |
| Chemical-space shift | `results_key/chemical_space_similarity_bins.csv`, `tanimoto_ad_metrics.csv` |
| Conformal coverage | `results_key/conformal_metrics.csv` |
| Chemically discordant collisions | `results_v31/chemically_discordant_collision_pairs.csv` |
| RF feature-sampling sensitivity | `results_v31/rf_sensitivity_metrics.csv`, `rf_sensitivity_summary.csv` |

## Scientific scope and limitations

This repository supports an **EGFR-specific representation study**, not a claim that Graph19 or ECFP4 universally dominates for all targets or learners. Graph19 intentionally omits atom labels, bond orders, 3D geometry, protein state and assay context. Descriptor-vector collisions therefore quantify ambiguity in this specific compact representation and do not imply that the underlying unlabeled graphs are non-isomorphic. Scaffold conformal coverage is reported as a distribution-shift stress test rather than a formal exchangeability-guaranteed coverage result.

## Release status

**Current status: pre-submission freeze.**  
The repository is being frozen for the `v1.0.0-submission` release. That release will be archived through Zenodo together with the complete frozen data/intermediate payload. The resulting Zenodo DOI will then be inserted into the manuscript, `CITATION.cff`, and this README before submission.

## Authors

- **Sunilgar L. Gusai** - corresponding author, Marwadi University
- **Manoharsinh R. Jadeja** - Marwadi University

## License

Original code in this repository is released under the **MIT License**. ChEMBL-derived data remain subject to the ChEMBL licensing terms described in `DATA_LICENSE_NOTICE.md`.

## Full immutable archive

The GitHub clone is intentionally lightweight. The submission-linked Zenodo archive contains the complete frozen ChEMBL-derived raw/intermediate/final data, the frozen V3 return archive, all publication figure formats, and the complete machine-readable supporting output set. This separation keeps normal cloning practical while preserving third-party reproducibility in the immutable archive.
