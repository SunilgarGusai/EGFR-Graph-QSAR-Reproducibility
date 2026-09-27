# PAPER003 V3.1 - Five-cycle pre-submission review

Target: Journal of Cheminformatics  
Review date: 2026-09-27

## Cycle 1 - Journal scope, article structure and formal requirements
Status: PASS after revision.

- Scientific scope matches molecular representation, QSAR, data mining, applicability/reliability and reproducible cheminformatics.
- Abstract remains below 350 words and contains a two-sentence Scientific Contribution section.
- Manuscript uses double-line spacing plus line/page numbering.
- Declarations include availability, competing interests, funding, contributions and acknowledgements.
- LLM use is documented in Methods, as required by the current journal Research instructions.
- Graphical abstract is 920 x 300 px, white background, and <150 KB.
- Additional files are explicitly described in the manuscript.
- Remaining submission gate: insert immutable Zenodo DOI and confirm final co-author approval.

## Cycle 2 - Novelty, positioning and claim defensibility
Status: PASS after major reframing.

- Dataset size and EGFR modeling are not presented as novelty because a 2026 multi-target benchmark already includes approximately 10k EGFR compounds.
- Endpoint wording changed from binding affinity to pIC50 / inhibitory potency; IC50 is not presented as a thermodynamic affinity constant.
- QSPR framing changed to QSAR.
- Broad intrinsic "information content" wording was replaced by predictive signal accessible under matched RF evaluation.
- No first/unprecedented claim is made.
- The differentiated contribution is extreme compact representation + exact descriptor-vector degeneracy + matched random/scaffold comparison + chemical-space reliability + uncertainty.

## Cycle 3 - Methods, statistics, leakage and reproducibility
Status: PASS.

- Frozen lineage reproduced exactly: 17,705 -> 17,412 -> 10,056.
- All 19 graph descriptors independently regenerated within 5e-8 of the archive.
- 3,579 Murcko scaffolds recovered.
- Same compounds/splits/folds/model principles used across representations.
- Paired representation effects use matched split/fold differences with 95% Student-t CIs (n=5).
- RF feature-sampling sensitivity at max_features=0.5 preserves ECFP4 > RDKit2D > Graph19 under random and scaffold evaluation.
- Selected max_features=1.0 checks preserve the same ordering but are explicitly treated as spot checks.
- Exact Graph19 collisions are operationally defined at 8-decimal regenerated-vector precision; they are not mislabelled as graph-isomorphism collisions.
- Conformal scaffold results are distribution-shift stress tests, not formal exchangeability-guaranteed coverage.
- Williams threshold corrected to training-set h*=0.00745898.

## Cycle 4 - Figures, tables and explanatory power
Status: PASS after redesign.

- Fig. 1 rebuilt as a clean six-stage workflow.
- Fig. 2 is the signature compactness/performance/generalization figure; paired effect sizes are reported in the text and supplement.
- Fig. 3 combines family ablation and held-out individual permutation importance.
- Fig. 4 makes representation ambiguity concrete using collision-group distributions and representative chemically discordant molecular pairs.
- Fig. 5 combines chemical-space error behavior with conformal reliability under scaffold shift.
- Williams diagnostics and RF-sampling detail moved to Supplementary Information.
- Main figures are supplied as vector PDF plus high-resolution PNG; captions are outside the graphic in the manuscript.

## Cycle 5 - Repository, data/software availability and submission ecosystem
Status: PASS for pre-submission; archival DOI pending.

- Public GitHub repository: https://github.com/SunilgarGusai/EGFR-Graph-QSAR-Reproducibility
- Code license: MIT.
- ChEMBL-derived data licensing/attribution is separated from the code license.
- Frozen inputs, manifests/checksums, pipeline code, V3.1 code, key result tables, figures and environment metadata are prepared for archival release.
- README exposes the scientific question, key results, reproduction route, provenance and limitations.
- Final archival step: freeze release v1.0.0-submission, archive with Zenodo, insert DOI in manuscript/CITATION.cff/README, then compile the immutable submission version.

## Editorial-style verdict
The strengthened study is materially above the former V2/V3 baseline. Its most distinctive evidence is not the absolute predictive score but the direct demonstration that a 19-variable unlabeled topological encoding creates widespread exact representation aliasing, including chemically dissimilar molecules with multi-log-unit potency differences, while the representation hierarchy remains robust under scaffold holdout and RF feature-sampling sensitivity. No further scientific model run is currently required; the remaining gate is archival freeze and final author approval.
