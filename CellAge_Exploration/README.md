# CellAge — Score provenance and current-model state analysis

The historical 20261030 snapshot contains hUSI scores for 1,665,482/536,017 fibroblast/endothelial-lineage cells and recSID recommendations for 1,372,858/352,228 cells. Its K8000 partial clock branch is not the current compact clock. hUSI is a senescence-associated transcriptomic score, not age in years; recSID dimensions are identity recommendations, not verified senescence intensities.

The current analysis of the 27 selected states per lineage is published under [F3 ModelFeature_CellAge](../F3_ds_analysis/F3_features/ModelFeature_cellAge/). It covers 731/529 adult sample proxies before state-specific primary coverage filtering. Relative hUSI state contrasts and adjusted age slopes are different tests; no primary adjusted age-dependent hUSI/SSE/composition comparison survives its metric-specific correction.

![Current selected-state CellAge](../F3_ds_analysis/F3_features/ModelFeature_cellAge/figures/Fig_ModelFeature_CellAge_Integration.png)

## Historical snapshot

![Historical CellAge figure](figures/Fig_CellAge_Stromal_Integration.png)

The historical figure is retained as provenance. Its equal-state relative comparator and older clock cannot be substituted for the current cell-weighted selected-state analysis. [Historical aggregate tables](tables/) · [Current Results Part 3](../docs/results/Results_Part3_StableFeatures_Mechanisms.md) · [Evidence](evidence/evidence_registry.csv)
