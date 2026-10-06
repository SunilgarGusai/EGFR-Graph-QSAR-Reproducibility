<p align="center">
  <img src="docs/assets/repository-banner.svg" alt="EGFR Graph QSAR Reproducibility" width="100%" />
</p>

<h1 align="center">EGFR Graph QSAR Reproducibility</h1>

<p align="center"><strong>Mechanisms of Representation Degeneracy and Generalization Limits of Classical Molecular Graph Descriptors in EGFR QSAR</strong></p>

## Study question

This repository supports a frozen 10,056-compound human EGFR pIC50 study asking what predictive signal remains after molecular structure is compressed to 19 classical graph invariants (Graph19), which chemical distinctions are lost by that compression, and how the loss relates to scaffold generalization and reliability.

The study compares Graph19 with 196 RDKit2D descriptors and 2,048-bit ECFP4 fingerprints on the same compounds and fixed random/scaffold partitions. Random Forest is the principal matched learner; ExtraTreesRegressor is used only as a second-learner sensitivity check.

## Principal benchmark

| Representation | Features | Random R2 mean +/- SD | Scaffold R2 mean +/- SD |
|---|---:|---:|---:|
| Graph19 | 19 | 0.392 +/- 0.038 | 0.211 +/- 0.136 |
| RDKit2D | 196 | 0.657 +/- 0.017 | 0.518 +/- 0.100 |
| ECFP4 | 2,048 | **0.735 +/- 0.014** | **0.615 +/- 0.062** |

The scientific conclusion is not that Graph19 should outperform modern fingerprints. Graph19 retains measurable signal with extreme compression, but richer chemical representations recover substantially more predictive and extrapolative information.

## V3.3 mechanism-resolved representation diagnosis

At the operational eight-decimal Graph19 encoding, the frozen cohort contains **1,507 collision groups involving 4,181 molecules (41.58%)**. The new graph-isomorphism/canonical-topology audit found that every such collision group contains only one hydrogen-suppressed unlabeled simple-graph topology. In this cohort, the observed degeneracy is therefore introduced before any additional collision between different unlabeled topologies is required.

Primary earliest information-loss mechanisms are:

| Mechanism | Groups | Molecules |
|---|---:|---:|
| Atom/element-label loss | **1,278** | **3,691** |
| Bond-order/aromaticity loss | 104 | 222 |
| Stereo/other full-identity loss | 125 | 268 |
| Different unlabeled simple topology | **0** | **0** |

Precision sensitivity shows 4,181 collided molecules at 6 and 8 decimals, 4,180 at 10 decimals, 4,169 at 12 decimals, and 3,719 under exact binary64 equality. Rounded equality is therefore treated as an operational encoding rather than mathematical real-number equality.

## Whole-vector ECFP4 comparison

The same cohort contains **307 identical complete ECFP4 collision groups involving 702 molecules (6.98%)**. ECFP4 is not described as collision-free; it is substantially more discriminating than Graph19 on this cohort. This whole-vector analysis is distinct from within-molecule Morgan hash-folding collisions.

## Collision-conditioned held-out error

Graph19 ambiguity is not uniformly associated with large error. Under scaffold evaluation, Graph19 MAE is approximately **0.636** when the matching collision vector is present in training, **1.018** when the collision vector is unseen, and **1.186** in the post-hoc high-response-disagreement subgroup. RDKit2D and ECFP4 reduce the corresponding penalty for the difficult subgroups. These are associative diagnostics, not causal claims.

## Second-learner sensitivity

ExtraTreesRegressor preserves the representation ordering:

| Split | Graph19 R2 | RDKit2D R2 | ECFP4 R2 |
|---|---:|---:|---:|
| Random | 0.366 +/- 0.040 | 0.692 +/- 0.020 | **0.729 +/- 0.018** |
| Scaffold | 0.163 +/- 0.165 | 0.556 +/- 0.099 | **0.602 +/- 0.070** |

Thus the qualitative ordering **ECFP4 > RDKit2D > Graph19** is not specific to Random Forest.

## V3.3 repository map

- `results/v33/` - machine-readable V3.3 summaries for Actions A-E.
- `docs/V33_RESULTS_ADDENDUM.md` - concise scientific addendum.
- `v33/config_v33.json` - frozen incremental-analysis configuration.
- `v33/requirements-lock-v33.txt` - pinned V3.3 Python/RDKit/scikit-learn environment.
- `v33/SOURCE_FREEZE.json` - source-package hashes and immutable V3.2 provenance.
- `v33/EGFR_QSAR_V33_PUBLIC_PAYLOAD.zip` - bundled V3.3 strengthening source and reproducibility materials.
- `results/frozen/` and `results/v31/` - earlier frozen benchmark, reliability and robustness results retained for provenance.

## Data and code availability

The repository publicly provides the analysis code/reproducibility bundle, pinned environment information, provenance records, checksums and machine-readable V3.3 result summaries. The modeled cohort is a ChEMBL-derived dataset; its exact frozen checksum and source lineage are documented under `data/` and `docs/DATA_PROVENANCE.md`. For journal submission, the private supplementary package can additionally include the frozen machine-readable 10,056-compound cohort, collision memberships, annotated held-out predictions, benchmark/uncertainty tables and result manifests without changing the frozen scientific population.

Original project code is released under the MIT License. ChEMBL-derived material remains subject to its source licensing and attribution requirements. The historical ChEMBL release identifier was not preserved and is disclosed rather than guessed.

## Key reproducibility boundaries

- pIC50 is an inhibitory-potency endpoint, not a thermodynamic binding-affinity constant.
- Graph19 intentionally omits atom labels, bond orders, stereochemistry, 3D geometry, protein state and assay context.
- Eight-decimal Graph19 equality is operational; exact binary64 results are reported separately.
- The absence of non-isomorphic-topology collisions is empirical for this EGFR cohort, not a theorem that Graph19 is globally injective.
- ECFP4 is not collision-free.
- Collision-conditioned error results are associative and do not establish causality.
- Scaffold conformal results are distribution-shift stress tests, not exchangeability-guaranteed coverage statements.

## Authors

1. **Sunilgar L. Gusai** - Faculty of Computer Applications, Marwadi University, Rajkot 360003, Gujarat, India; ORCID: https://orcid.org/0009-0004-0739-4812
2. **Manoharsinh R. Jadeja** - Department of Artificial Intelligence, Machine Learning and Data Science, Marwadi University, Rajkot 360003, Gujarat, India
3. **Dharmeshkumar Shah** - corresponding author; Information and Library Network (INFLIBNET) Centre, Infocity, Gandhinagar 382007, Gujarat, India; email: dashah@inflibnet.ac.in
4. **Ripal Ranpara** - Faculty of Computer Applications, Marwadi University, Rajkot 360003, Gujarat, India

## Citation

Citation metadata are supplied in [`CITATION.cff`](CITATION.cff). Please cite the associated article when final bibliographic metadata become available.

## License and usage

Original project code is released under the [`MIT License`](LICENSE). ChEMBL-derived material retains its upstream licensing and attribution requirements; see the repository provenance documentation for the exact public-data boundary.
