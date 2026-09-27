# Frozen experiment protocol

## Scientific question

The study asks how much EGFR pIC50 predictive signal remains accessible when molecular structure is compressed to 19 classical graph invariants, where that representation fails, and how those failures compare with contemporary molecular descriptors and fingerprints under chemical-space shift.

## Frozen endpoint and cohort

- Target: human EGFR (`CHEMBL203`).
- Activity endpoint: IC50 transformed/reported as pIC50.
- Assay class: ChEMBL binding assays (`assay_type = B`).
- Source records: 17,705.
- Stage-1 curated records: 17,412.
- Final deduplicated cohort: 10,056 unique compounds.
- Duplicate handling: InChIKey grouping, median pIC50 retained.
- Primary label source: ChEMBL `pchembl_value` when present; compatible unit-derived pIC50 fallback otherwise.

The historical ChEMBL release identifier was not preserved in the archived source file. The frozen raw file checksum, exact curation logic, source target identifier, and audit counts are therefore used as the reproducibility anchors rather than reconstructing a release number by guesswork.

## Molecular representations

1. **Graph19** — 19 classical graph invariants spanning degree-, distance-, and spectral-based families.
2. **RDKit2D** — 196 retained scalar descriptors after finite/constant filtering in the frozen environment.
3. **ECFP4** — binary Morgan fingerprint, radius 2, 2,048 bits, no chirality.

The principal representation comparison uses the same frozen compounds, target values, split assignments, metrics, and Random Forest principles across all representations.

## Validation families

### Random split benchmark

Five random seeds: `42, 101, 202, 303, 404`.

### Scaffold benchmark

Five Bemis-Murcko GroupKFold partitions derived from 3,579 unique scaffolds.

### Historical assurance model

A provenance-backed 600-tree Random Forest reproduces the V2-style external holdout result and is used for validation continuity, 10-fold CV, Y-randomization, and scaffold assurance.

## Reliability and interpretation analyses

- nearest-training ECFP4 Tanimoto similarity and error-by-similarity analysis;
- Williams leverage diagnostics with training-set leverage threshold;
- 90% split-conformal intervals, with scaffold coverage reported only as a distribution-shift stress test;
- Graph19 descriptor-family ablation;
- held-out individual and grouped permutation importance;
- exact Graph19 descriptor-vector collision analysis;
- paired representation effect sizes across matched splits/folds;
- RF feature-sampling sensitivity using a common 50% feature fraction and selected full-feature spot checks.

## Reporting boundary

The study does **not** claim that pIC50 is a thermodynamic binding-affinity constant, that Graph19 collisions imply graph isomorphism, or that the observed representation ranking is universal across targets and learners.
