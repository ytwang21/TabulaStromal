# Source differences and interpretation controls

| Issue | Resolution used in this release |
| --- | --- |
| Four named paths versus the wording “three result directories” | All four named directories were audited and given illustrated result pages. |
| Historical StromalAtlas name | Current source directory is F1_StromalAtlas. Exported text paths are generalized; original source hashes remain recorded. |
| Full atlas versus age-clock cohort | 5,916,760/2,239,277 cells belong to the full atlas. Clock development has 640/437 sample units; the counts are not interchangeable. |
| Sample versus donor | No independently harmonized donor identity was found. Historical donor-named fields retain their source labels and are explained as sample-designated units. |
| Tissue versus annotation suffix | Raw metadata have 39 labels per lineage. Figure 1C–E use 12/11 terminal annotation suffixes including pan, not 12/11 anatomical organs. |
| Figure 2 versus fivefold benchmark | Figure 2 is source fold 4 with internally selected K. Fivefold means/SDs come from the complete benchmark tables. |
| Global winner versus final EN | XGBoost/ExtraTrees win the internal descriptive benchmark; the final interpretable branch is EN at K=20000/10000. |
| External scopes | No adaptation, matched held-out adaptation and adapted ALL100 are separate rows. ALL100 includes fitted observations. |
| MAD | Source MAD is mean absolute deviation of chronological-age labels from their mean. MedAE is a separate prediction-error metric. |
| Source “biological age” wording | The actual trained target is recorded chronological age. The source figure is copied unchanged and the caption explains the boundary. |
| CellAge version | The older K8000 prediction/stability branch is not substituted for the completed final clocks. |
| SenCID recSID | Available recommendation columns do not establish genuine senescence intensity. |
| Mechanisms | Coefficients, enrichment and expression–age associations identify candidates; no causal driver or therapeutic target is claimed. |
| Conceptual overview | Illustrative r/MAE/RMSE, healthy-only wording, donor icons and GWAS/TF sketches are not treated as corresponding empirical results. |

No prettier historical metric was substituted for a poorer final-protocol result. Negative external R² and non-significant adjusted CellAge associations are retained.
