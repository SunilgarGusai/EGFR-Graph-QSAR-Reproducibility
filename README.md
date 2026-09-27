# EGFR Graph QSAR Reproducibility

Reproducibility repository for the pre-submission study on **classical molecular graph descriptors versus contemporary molecular representations for EGFR pIC50 prediction under scaffold and chemical-space shift**.

> **Status:** final pre-submission validation (V3.1 strengthening in progress). The frozen submission release will be tagged only after the final scientific and visual audit.

## Scientific question

How much EGFR pIC50 predictive signal remains accessible when molecular structure is compressed to only **19 classical mathematical graph invariants**, and what accuracy, generalization and reliability costs arise relative to modern RDKit2D descriptors and ECFP4 fingerprints?

## Core design

- Frozen ChEMBL-derived human EGFR IC50/pIC50 cohort
- 10,056 unique compounds after reproducible curation and InChIKey aggregation
- 19 classical degree-, distance- and spectral graph invariants (Graph19)
- Matched comparison with RDKit2D and 2,048-bit ECFP4
- Repeated random splits and Bemis–Murcko scaffold evaluation
- Chemical-space applicability diagnostics
- Descriptor-family ablation and held-out permutation importance
- Split-conformal reliability analysis
- V3.1 representation-ambiguity / descriptor-collision analysis
- Reproducible checkpoints, machine-readable outputs and figure regeneration

## Frozen primary benchmark

| Representation | Features | Random-split R² | Scaffold R² |
|---|---:|---:|---:|
| Graph19 | 19 | 0.392 ± 0.038 | 0.211 ± 0.136 |
| RDKit2D | 196 | 0.657 ± 0.017 | 0.518 ± 0.100 |
| ECFP4 | 2,048 | 0.735 ± 0.014 | 0.615 ± 0.062 |

These values are the frozen V3 matched-Random-Forest benchmark. The V3.1 sensitivity analysis does **not** replace this benchmark; it tests whether the representation hierarchy is robust to RF feature-sampling choices.

## Repository roadmap

The submission release will contain:

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── environment/
├── data/
│   ├── README.md
│   ├── manifests/
│   └── frozen_processed/
├── src/
│   ├── curation/
│   ├── descriptors/
│   ├── representations/
│   ├── modeling/
│   ├── reliability/
│   └── v31_ambiguity/
├── configs/
├── splits/
├── results/
│   ├── tables/
│   └── machine_readable/
├── figures/
├── supplementary/
└── reproduce/
    ├── Windows/
    └── README.md
```

## Reproducibility principle

The repository is being prepared so that the manuscript's principal tables and figures can be traced to frozen inputs, explicit split assignments, fixed software/configuration metadata and deterministic scripts. Dataset redistribution will follow the applicable ChEMBL terms and attribution requirements; where redistribution is not appropriate, acquisition/curation instructions and checksums will be supplied.

## Authors

- **Sunilgar L. Gusai** — corresponding author
- **Manoharsinh R. Jadeja**

Marwadi University, Rajkot, Gujarat, India.

## License

Original code in this repository is released under the MIT License unless a file states otherwise. Third-party data and software retain their respective licenses.

## Citation

Citation metadata will be finalized in `CITATION.cff` when the manuscript DOI is available.
