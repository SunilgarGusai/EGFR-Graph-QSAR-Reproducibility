from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
import zipfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
V31 = Path(__file__).resolve().parent
CFG = json.loads((V31 / "v31_config.json").read_text())
OUT = ROOT / "results_v31_run"
CACHE = ROOT / "cache" / "v31"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jobs() -> int:
    n = os.cpu_count() or 4
    return max(1, min(int(CFG["n_jobs_cap"]), n - int(CFG["reserve_cpu_cores"])))


def locate_v3_return() -> Path:
    candidates = [
        ROOT / "inputs" / "EGFR_QSAR_REPRODUCIBILITY_RETURN.zip",
        V31 / "EGFR_QSAR_REPRODUCIBILITY_RETURN.zip",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Place EGFR_QSAR_REPRODUCIBILITY_RETURN.zip in inputs/ or v31/. "
        "The immutable submission archive contains this frozen intermediate."
    )


def extract_v3() -> Path:
    src = locate_v3_return()
    dest = CACHE / "v3_return"
    marker = dest / ".sha256"
    sig = sha256(src)
    if not marker.exists() or marker.read_text().strip() != sig:
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        with zipfile.ZipFile(src) as z:
            z.extractall(dest)
        marker.write_text(sig)
    return dest / "results" / "full"


def load_inputs(base: Path):
    graph = pd.read_csv(base / "02_descriptors" / "graph19_regenerated.csv")
    rdkit2d = pd.read_csv(base / "03_representations_splits" / "rdkit2d.csv")
    ecfp = np.load(base / "03_representations_splits" / "morgan_ecfp4_2048.npz")["X"].astype(np.float32)
    scaffolds = pd.read_csv(base / "03_representations_splits" / "scaffold_assignments.csv")
    benchmark = pd.read_csv(base / "04_models" / "representation_benchmark_metrics.csv")
    assert len(graph) == 10056 and len(rdkit2d) == 10056 and ecfp.shape == (10056, 2048)
    return graph, rdkit2d, ecfp, scaffolds, benchmark


def graph19_collisions(graph: pd.DataFrame, ecfp: np.ndarray):
    out = OUT / "01_graph19_ambiguity"
    out.mkdir(exist_ok=True)
    descriptors = list(graph.columns[3:])
    rounded = np.round(graph[descriptors].to_numpy(float), int(CFG["round_decimals_collision"]))
    keys = pd.Series(
        [hashlib.sha256(row.tobytes()).hexdigest() for row in rounded],
        index=graph.index,
        name="graph19_key",
    )
    meta = graph[["InChIKey", "SMILES", "pIC50"]].copy()
    meta["graph19_key"] = keys
    meta["molecular_formula"] = [
        rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s)) for s in meta.SMILES
    ]
    meta["collision_group_size"] = keys.map(keys.value_counts()).to_numpy()

    groups = (
        meta.groupby("graph19_key")
        .agg(
            n=("pIC50", "size"),
            n_formula=("molecular_formula", "nunique"),
            pIC50_min=("pIC50", "min"),
            pIC50_max=("pIC50", "max"),
            pIC50_mean=("pIC50", "mean"),
            pIC50_median=("pIC50", "median"),
        )
        .reset_index()
    )
    groups["pIC50_range"] = groups.pIC50_max - groups.pIC50_min
    collisions = groups[groups.n > 1].copy().sort_values("pIC50_range", ascending=False)

    pair_rows = []
    for key, idxs in meta.groupby("graph19_key").indices.items():
        if len(idxs) < 2:
            continue
        for i, j in combinations(list(idxs), 2):
            a, b = ecfp[i].astype(bool), ecfp[j].astype(bool)
            union = np.logical_or(a, b).sum()
            sim = float(np.logical_and(a, b).sum() / union) if union else 1.0
            pair_rows.append(
                {
                    "graph19_key": key,
                    "row_i": int(i),
                    "row_j": int(j),
                    "tanimoto_ecfp4": sim,
                    "delta_pIC50": float(abs(meta.loc[i, "pIC50"] - meta.loc[j, "pIC50"])),
                    "InChIKey_A": meta.loc[i, "InChIKey"],
                    "formula_A": meta.loc[i, "molecular_formula"],
                    "SMILES_A": meta.loc[i, "SMILES"],
                    "pIC50_A": float(meta.loc[i, "pIC50"]),
                    "InChIKey_B": meta.loc[j, "InChIKey"],
                    "formula_B": meta.loc[j, "molecular_formula"],
                    "SMILES_B": meta.loc[j, "SMILES"],
                    "pIC50_B": float(meta.loc[j, "pIC50"]),
                }
            )
    pairs = pd.DataFrame(pair_rows)
    strict = pairs[
        (pairs.tanimoto_ecfp4 < float(CFG["ecfp_tanimoto_threshold"]))
        & (pairs.delta_pIC50 > float(CFG["activity_delta_threshold"]))
    ].sort_values(["delta_pIC50", "tanimoto_ecfp4"], ascending=[False, True])

    means = meta.groupby("graph19_key").pIC50.transform("mean").to_numpy()
    medians = meta.groupby("graph19_key").pIC50.transform("median").to_numpy()
    y = meta.pIC50.to_numpy(float)
    collided = meta.collision_group_size.to_numpy() > 1
    sst = float(np.sum((y - y.mean()) ** 2))
    sse = float(np.sum((y - means) ** 2))
    summary = {
        "n_molecules": int(len(meta)),
        "n_unique_graph19_vectors": int(meta.graph19_key.nunique()),
        "n_collision_groups": int(len(collisions)),
        "n_molecules_in_collision_groups": int(collisions.n.sum()),
        "fraction_molecules_in_collision_groups": float(collisions.n.sum() / len(meta)),
        "collision_groups_with_distinct_molecular_formula": int((collisions.n_formula > 1).sum()),
        "fraction_collision_groups_with_distinct_molecular_formula": float((collisions.n_formula > 1).mean()),
        "groups_pic50_range_gt_1": int((collisions.pIC50_range > 1).sum()),
        "groups_pic50_range_gt_2": int((collisions.pIC50_range > 2).sum()),
        "groups_pic50_range_gt_3": int((collisions.pIC50_range > 3).sum()),
        "pairs_tanimoto_lt_0_5_and_delta_pic50_gt_2": int(((pairs.tanimoto_ecfp4 < 0.5) & (pairs.delta_pIC50 > 2)).sum()),
        "exact_collision_dataset_rmse_floor": float(np.sqrt(np.mean((y - means) ** 2))),
        "exact_collision_dataset_mae_floor": float(np.mean(np.abs(y - medians))),
        "exact_collision_in_sample_r2_ceiling": float(1 - sse / sst),
        "collided_subset_rmse_floor": float(np.sqrt(np.mean((y[collided] - means[collided]) ** 2))),
        "collided_subset_mae_floor": float(np.mean(np.abs(y[collided] - medians[collided]))),
    }
    meta.to_csv(out / "graph19_collision_membership.csv", index=False)
    collisions.to_csv(out / "graph19_collision_groups.csv", index=False)
    pairs.to_csv(out / "graph19_collision_pairs.csv", index=False)
    strict.to_csv(out / "chemically_discordant_collision_pairs.csv", index=False)
    strict.head(int(CFG["representative_pairs_to_draw"])).to_csv(out / "representative_collision_pairs.csv", index=False)
    (out / "graph19_ambiguity_summary.json").write_text(json.dumps(summary, indent=2))


def paired_effects(benchmark: pd.DataFrame):
    out = OUT / "02_paired_effects"
    out.mkdir(exist_ok=True)
    rows = []
    for split_type in ["random", "scaffold"]:
        for a, b in [("RDKit2D", "Graph19"), ("ECFP4_2048", "Graph19"), ("ECFP4_2048", "RDKit2D")]:
            aa = benchmark[(benchmark.representation == a) & (benchmark.split_type == split_type)].sort_values("split_id")
            bb = benchmark[(benchmark.representation == b) & (benchmark.split_type == split_type)].sort_values("split_id")
            assert aa.split_id.tolist() == bb.split_id.tolist()
            for metric in ["r2", "rmse", "mae"]:
                delta = aa[metric].to_numpy() - bb[metric].to_numpy()
                n = len(delta)
                mean = float(delta.mean())
                sd = float(delta.std(ddof=1))
                half = float(t.ppf(0.975, n - 1) * sd / np.sqrt(n))
                rows.append({
                    "split_type": split_type,
                    "representation_A": a,
                    "representation_B": b,
                    "metric": metric,
                    "n_pairs": n,
                    "mean_delta_A_minus_B": mean,
                    "sd_delta": sd,
                    "ci95_low": mean - half,
                    "ci95_high": mean + half,
                    "n_positive": int((delta > 0).sum()),
                    "n_negative": int((delta < 0).sum()),
                })
    pd.DataFrame(rows).to_csv(out / "paired_representation_effect_sizes.csv", index=False)


def score(y, pred):
    return {
        "r2": float(r2_score(y, pred)),
        "rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "mae": float(mean_absolute_error(y, pred)),
    }


def fit_rf(X, y, train, test, representation, split_type, split_id, mode, max_features, trees, seed):
    model = RandomForestRegressor(
        n_estimators=int(trees),
        max_features=max_features,
        min_samples_leaf=1,
        random_state=int(seed),
        n_jobs=jobs(),
    )
    start = time.perf_counter()
    model.fit(X[train], y[train])
    train_seconds = time.perf_counter() - start
    pred = model.predict(X[test])
    return {
        "mode": mode,
        "representation": representation,
        "split_type": split_type,
        "split_id": int(split_id),
        "max_features": max_features,
        "n_estimators": int(trees),
        **score(y[test], pred),
        "train_seconds": train_seconds,
    }


def rf_sensitivity(graph, rdkit2d, ecfp, scaffolds):
    out = OUT / "03_rf_sensitivity"
    out.mkdir(exist_ok=True)
    y = graph.pIC50.to_numpy(float)
    idx = np.arange(len(y))
    reps = {
        "Graph19": graph.iloc[:, 3:].to_numpy(float),
        "RDKit2D": rdkit2d.iloc[:, 3:].to_numpy(float),
        "ECFP4_2048": ecfp,
    }
    rows = []
    for rep, X in reps.items():
        for seed in CFG["random_seeds"]:
            tr, te = train_test_split(idx, test_size=0.2, random_state=int(seed), shuffle=True)
            rows.append(fit_rf(X, y, tr, te, rep, "random", seed, "fraction_0.5", 0.5, CFG["rf_trees_fraction_sensitivity"], seed))
        for fold in sorted(scaffolds.fold.unique()):
            te = scaffolds.index[scaffolds.fold == fold].to_numpy()
            tr = scaffolds.index[scaffolds.fold != fold].to_numpy()
            rows.append(fit_rf(X, y, tr, te, rep, "scaffold", fold, "fraction_0.5", 0.5, CFG["rf_trees_fraction_sensitivity"], 42 + int(fold)))
        for seed in CFG["rf_full_feature_spotcheck_random_seeds"]:
            tr, te = train_test_split(idx, test_size=0.2, random_state=int(seed), shuffle=True)
            rows.append(fit_rf(X, y, tr, te, rep, "random", seed, "full_feature_spotcheck", 1.0, CFG["rf_trees_full_feature_spotcheck"], seed))
        for fold in CFG["rf_full_feature_spotcheck_scaffold_folds"]:
            te = scaffolds.index[scaffolds.fold == fold].to_numpy()
            tr = scaffolds.index[scaffolds.fold != fold].to_numpy()
            rows.append(fit_rf(X, y, tr, te, rep, "scaffold", fold, "full_feature_spotcheck", 1.0, CFG["rf_trees_full_feature_spotcheck"], 42 + int(fold)))
    df = pd.DataFrame(rows)
    df.to_csv(out / "rf_sensitivity_metrics.csv", index=False)
    summary = df.groupby(["mode", "representation", "split_type"])[["r2", "rmse", "mae", "train_seconds"]].agg(["mean", "std", "count"]).reset_index()
    summary.columns = ["__".join([str(x) for x in c if x]) if isinstance(c, tuple) else c for c in summary.columns]
    summary.to_csv(out / "rf_sensitivity_summary.csv", index=False)


def manifest():
    rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file() and path.name != "results_manifest.json":
            rows.append({"path": str(path.relative_to(OUT)).replace("\\", "/"), "sha256": sha256(path), "bytes": path.stat().st_size})
    (OUT / "results_manifest.json").write_text(json.dumps(rows, indent=2))


def main():
    base = extract_v3()
    graph, rdkit2d, ecfp, scaffolds, benchmark = load_inputs(base)
    graph19_collisions(graph, ecfp)
    paired_effects(benchmark)
    rf_sensitivity(graph, rdkit2d, ecfp, scaffolds)
    manifest()
    (OUT / "RUN_COMPLETE.txt").write_text("EGFR Graph QSAR V3.1 repository run completed successfully.\n")
    print(f"Complete. Outputs: {OUT}")


if __name__ == "__main__":
    main()
