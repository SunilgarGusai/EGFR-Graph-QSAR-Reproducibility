# Data provenance

## Frozen lineage

```text
17,705 archived ChEMBL-derived EGFR IC50 records
        ↓
17,412 curated binding-assay activity records
        ↓
10,056 unique compounds after InChIKey grouping
        ↓
median pIC50 per compound
```

The frozen source is EGFR (`CHEMBL203`) and the retained assay class is ChEMBL binding assays (`B`). The archived raw file contains the source activity, assay, document, molecule, relation, units, canonical SMILES, normalized IC50 and pIC50 fields needed to reproduce the historical curation.

## Label provenance

Of the 17,412 curated activity rows:

- **17,156** use ChEMBL `pchembl_value` as the primary pIC50 value;
- **256** use a compatible unit-derived fallback.

Where both sources were present, the archived audit found only small rounding-scale differences.

## Deduplication

RDKit InChIKeys are used to group repeated molecule records. Median pIC50 is retained for each final compound. The reconstructed pipeline reproduces all 10,056 archived final InChIKeys and canonical SMILES.

## Descriptor provenance

All 19 Graph19 descriptors were independently regenerated from the frozen structures and matched the archived descriptor matrix to approximately `5e-8` maximum absolute discrepancy.

## Historical release limitation

The original ChEMBL release/version identifier was not preserved in the historical PAPER003 archive. This repository does not invent one. Instead it exposes the frozen-file checksums, target identifier, curation rules, counts, and source-field structure needed to identify exactly what was analyzed.

See `INPUT_MANIFEST.json` and `DATA_LICENSE_NOTICE.md` for file-level and licensing notes.
