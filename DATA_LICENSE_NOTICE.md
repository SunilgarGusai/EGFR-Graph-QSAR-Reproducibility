# Data license and provenance notice

The activity records used by this study are derived from **ChEMBL** (human EGFR / CHEMBL203) and are redistributed for reproducibility with attribution under the license stated by ChEMBL: **Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)**.

ChEMBL: https://www.ebi.ac.uk/chembl/

The repository-level MIT license applies only to original software authored for this study. It does **not** replace the licensing terms applicable to ChEMBL-derived data.

Frozen lineage used by the study:

- raw activity archive: 17,705 records;
- curated Stage-1 archive: 17,412 records;
- final deduplicated cohort: 10,056 unique compounds;
- target: CHEMBL203 (human EGFR);
- endpoint: IC50 transformed/reported as pIC50;
- replicate aggregation: median pIC50 by InChIKey.

Input SHA-256 hashes are recorded in `INPUT_MANIFEST.json` and the execution packages. The study intentionally preserves this frozen snapshot rather than silently refreshing the cohort to a later ChEMBL release.
