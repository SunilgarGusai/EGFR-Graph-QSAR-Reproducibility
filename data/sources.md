# Data sources

## Primary source

The frozen EGFR Graph QSAR cohort is derived from **ChEMBL** activity records for human EGFR (`CHEMBL203`) using IC50 activity records from binding assays (`assay_type = B`).

The historical source snapshot was archived locally during the original EGFR Graph QSAR work. The exact ChEMBL release identifier was not retained in that archive, so this repository does **not** infer or invent a version number. Reproducibility is anchored to the frozen-file checksums, source target identifier, curation rules, and row-count audit.

## Frozen lineage

- archived source records: 17,705
- curated activity rows: 17,412
- final unique compounds: 10,056
- deduplication key: RDKit InChIKey
- repeated measurements: median pIC50

## Attribution and terms

ChEMBL is an EMBL-EBI resource. ChEMBL-derived files remain subject to applicable ChEMBL/EMBL-EBI terms and attribution requirements. See `../DATA_LICENSE_NOTICE.md`.
