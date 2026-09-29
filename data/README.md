# Frozen EGFR data provenance

The definitive study uses a frozen ChEMBL-derived human EGFR pIC50 cohort containing 10,056 modeled compounds.

For the public submission repository, this directory records the exact data lineage and cryptographic checksums used to identify the frozen cohort. ChEMBL-derived material remains subject to its source licensing and attribution requirements.

The public repository focuses on code, provenance, environment locks and machine-readable result tables. The exact frozen modeled cohort and fixed random/scaffold split assignments are included in the submission machine-readable supplementary archive; they can be deposited publicly by the corresponding author or journal without changing the scientific population.

Key provenance anchors:

- final modeled cohort SHA-256: `d2cbd3eec417be9bca5debac435c8639dd35261042ef7034757cac190b76423f`
- source records: 17,705
- curated activity rows: 17,412
- final unique compounds: 10,056
- random seeds: 42, 101, 202, 303, 404
- scaffold evaluation: five fixed Bemis-Murcko folds

See `checksums.csv`, `sources.md`, `../INPUT_MANIFEST.json`, `../docs/DATA_PROVENANCE.md`, and `../DATA_LICENSE_NOTICE.md` for the complete provenance boundary.
