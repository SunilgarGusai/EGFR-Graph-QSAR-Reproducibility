from __future__ import annotations

# Public V3.3 strengthening source.
# The complete audited implementation used for the manuscript is bundled in
# v33/EGFR_QSAR_V33_PUBLIC_PAYLOAD.zip. This public entry point documents the
# frozen protocol and points to the full executable source in that bundle.

V33_PROTOCOL = {
    "cohort_n": 10056,
    "graph19_primary_precision_decimals": 8,
    "precision_decimals": [6, 8, 10, 12],
    "ecfp4": {"radius": 2, "n_bits": 2048, "use_chirality": False},
    "random_seeds": [42, 101, 202, 303, 404],
    "scaffold_folds": 5,
    "second_learner": "ExtraTreesRegressor",
    "actions": [
        "A_collision_mechanism_decomposition",
        "B_graph19_precision_sensitivity",
        "C_ecfp4_complete_vector_collisions",
        "D_collision_conditioned_heldout_error",
        "E_second_learner_representation_sensitivity",
    ],
}

if __name__ == "__main__":
    import json
    print(json.dumps(V33_PROTOCOL, indent=2))
    print("Full audited implementation: v33/EGFR_QSAR_V33_PUBLIC_PAYLOAD.zip")
