# V3.3 results addendum

This update adds mechanism-resolved Graph19 collision analysis, precision sensitivity, complete-vector ECFP4 collision comparison, collision-conditioned held-out prediction error, and ExtraTrees representation sensitivity.

Headline frozen findings:
- 1,507 eight-decimal Graph19 collision groups / 4,181 molecules.
- No collision group spanned different hydrogen-suppressed unlabeled simple-graph topologies.
- Primary earliest loss: atom/element labels 1,278 groups; bond-order/aromaticity 104; stereo/other full identity 125.
- Exact binary64 Graph19 equality: 1,408 groups / 3,719 molecules (36.98%).
- ECFP4 complete-vector equality: 307 groups / 702 molecules (6.98%).
- Scaffold Graph19 MAE: seen collision vector 0.636; unseen collision vector 1.018; high-disagreement class 1.186.
- ExtraTrees preserves ECFP4 > RDKit2D > Graph19.

