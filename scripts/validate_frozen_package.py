from __future__ import annotations
import csv, json, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REQUIRED = [
    "README.md",
    "CITATION.cff",
    "INPUT_MANIFEST.json",
    "docs/FROZEN_RESULTS.md",
    "docs/EXPERIMENT_PROTOCOL.md",
    "docs/DATA_PROVENANCE.md",
    "docs/REPRODUCIBILITY.md",
    "results/frozen/representation_benchmark_summary.csv",
    "results/frozen/descriptor_family_ablation_summary.csv",
    "results/frozen/chemical_space_similarity_bins.csv",
    "results/frozen/conformal_metrics.csv",
    "results/v31/graph19_ambiguity_summary.json",
    "results/v31/paired_representation_effect_sizes.csv",
    "results/v31/representative_collision_pairs.csv",
    "results/v31/rf_sensitivity_summary.csv",
    "v31/v31_strengthening.py",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    missing = [p for p in REQUIRED if not (ROOT / p).exists()]
    if missing:
        raise SystemExit("Missing required files:\n- " + "\n- ".join(missing))

    manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text(encoding="utf-8"))
    ambiguity = json.loads((ROOT / "results/v31/graph19_ambiguity_summary.json").read_text(encoding="utf-8"))
    if ambiguity.get("n_molecules") != 10056:
        raise SystemExit("Unexpected frozen molecule count in ambiguity summary")
    if ambiguity.get("n_collision_groups") != 1507:
        raise SystemExit("Unexpected collision-group count")

    with (ROOT / "results/frozen/representation_benchmark_summary.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 6:
        raise SystemExit("Representation benchmark should contain 6 aggregate rows")

    print("Frozen repository validation PASSED")
    print("Required files:", len(REQUIRED))
    print("Input manifest entries:", len(manifest))
    print("Frozen molecules:", ambiguity["n_molecules"])
    print("Exact Graph19 collision groups:", ambiguity["n_collision_groups"])
    print("Benchmark rows:", len(rows))


if __name__ == "__main__":
    main()
