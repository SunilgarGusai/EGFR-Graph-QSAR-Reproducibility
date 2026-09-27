# Data provenance and frozen-input policy

The definitive study uses a frozen ChEMBL-derived EGFR cohort.

For the public submission repository, the data directory documents provenance, checksums and source terms rather than assuming that the software license governs third-party molecular data.

- `sources.md` — source target, endpoint and frozen lineage
- `checksums.csv` — SHA-256 hashes of the archived PAPER003 source/intermediate files
- `../INPUT_MANIFEST.json` — machine-readable input manifest
- `../DATA_LICENSE_NOTICE.md` — ChEMBL attribution and data-usage notice
- `../docs/DATA_PROVENANCE.md` — complete provenance narrative

The full frozen CSV payload is retained by the authors and can be made public with the source notice if that redistribution route is chosen for the journal submission. The repository does not infer a missing historical ChEMBL release identifier.
