# Reproducibility

This repository is designed as an auditable research package rather than a code-only supplement.

## What is reproducible

The frozen workflow supports:

- input checksum validation;
- exact data-lineage reconstruction;
- regeneration and verification of all 19 graph descriptors;
- RDKit2D and ECFP4 representation generation;
- fixed random and Murcko-scaffold partitioning;
- matched Random Forest representation benchmarking;
- 10-fold CV and Y-randomization assurance;
- descriptor-family ablation and permutation importance;
- chemical-space and leverage diagnostics;
- split-conformal uncertainty analysis;
- exact Graph19 representation-degeneracy analysis;
- paired effect-size reporting;
- RF feature-sampling sensitivity analysis.

## Environment

The definitive V3 environment is documented by `requirements.txt`, `requirements-lock.txt`, configuration files, and machine-readable environment outputs. The V3.1 strengthening workflow also records its own environment.

## Frozen-result principle

Submission-supporting values are copied into `results/frozen/` and `results/v31/`. Exploratory reruns should not overwrite those files. Any future update should either reproduce them or be versioned separately with an explicit explanation.

## Data redistribution

ChEMBL-derived material remains subject to ChEMBL source terms. The repository includes source/provenance and licensing notes. If the frozen CSV payload is publicly redistributed, it should remain accompanied by `DATA_LICENSE_NOTICE.md` and the source attribution.

## Validation

The public QA workflow performs lightweight repository checks. Full scientific retraining is intentionally not run on every GitHub push because it is substantially heavier than repository-integrity validation.
