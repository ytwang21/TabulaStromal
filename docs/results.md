# Results, provenance and interpretation

This release reports aggregate values from the existing project results. It does not rerun the full study or infer missing statistics from illustrative figures. Every exported source file is indexed in `provenance/source_manifest.csv`; exported-file hashes are recorded separately.

## Result families

| Directory/table | Source and scope |
| --- | --- |
| `results/atlas_summary.csv` | Full-atlas metadata audit, original `meta_full_20261031_232542_2002510` run |
| `results/{Fibro,Endo}/benchmark/` | Original fivefold benchmark; fold metrics and mean/std summaries, all model/K combinations |
| `results/final_clock_results_summary.csv` | Final full-development ElasticNet branch; selected-K conditional CV and external scopes |
| `results/{Fibro,Endo}/stability/` | Original final-branch stability pairs, subtype/gene summaries and conditional internal metrics |
| `results/{Fibro,Endo}/enrichment/` | Aggregate test counts and versioned-library provenance; no causal inference |
| `results/historical_figure2/` | Historical fold-4 figure data; external FT30 ALL100 scope; not final-clock performance |
| `results/cellage/` | Run `20261030_220917_407544`; historical partial K=8,000 ElasticNet integration |

Paths in copied tables are generalized, but numerical results are retained. The released `final_model_path` names indicate expected local output locations, not downloadable weights. Per-sample prediction files are not included, so full independent recalculation of study metrics requires the authorized original archive. Synthetic metric-formula tests check implementation independently of that archive.

## Atlas scope

Full-atlas counts are 5,916,760 fibroblast cells (4,818 sample labels; 848 GSE labels) and 2,239,277 endothelial cells (3,553 sample labels; 748 GSE labels). No independent donor column was verified. Each lineage has 39 raw tissue labels, including nonstandard categories; this is not a count of 39 harmonized anatomical organs. Full-atlas disease metadata are heterogeneous.

The age-clock development cohorts have 640/437 sample units and 634/431 GSE/GSM proxy groups (fibroblast/endothelial). These are separate from the full-atlas totals. Existing metacell trajectories use another subset and cannot be treated as all atlas cells or temporal lineage tracing.

## Benchmark and final branch

The global internal benchmark winners are fibroblast XGBoost (K=20,000, mean Valid R²=0.802880, MAE=9.714395 years) and endothelial ExtraTrees (K=20,000, mean Valid R²=0.716193, MAE=11.544655 years). These descriptive selections are based on the same internal comparison and are not estimates from an additional independent model-selection cohort.

The final interpretable EN branch uses K=20,000/10,000, alpha=1.0, l1_ratio=0.1. Its original-feature CV results exceed those of reduced stable-feature refits, so both final clocks retain the original feature sets. Internal results are conditional on K selected using CV; they should not be described as fully nested, unbiased estimates of the entire K-selection procedure.

Full-development stability identifies 438 pairs / 416 unique genes / 5 subtypes in fibroblasts and 23 pairs / 23 genes / 1 subtype in endothelial cells at frequency ≥0.70 over 300 group bootstrap fits per scope. Five outer-training scopes plus full development give 1,800 bootstrap records per lineage. Full-development stability is explanatory, while fold-specific stability must be learned on training folds only. Stable features are associations, not perturbation-validated mechanisms.

## External scope

There are 93/77 combined external sample units. Primary no-adaptation final-EN performance is poor (R² −2.943346/−1.719572, MAE 29.322619/24.889889 years). This limitation is not hidden by plotting only correlations or an adapted pooled result.

FT30 uses 28/23 external samples for adaptation and 65/54 for held-out evaluation. Held-out R² is −0.000803/−0.355121 and MAE 13.140648/13.175828 years. `FT30_all100_*` pools adaptation and held-out samples (R² 0.357544/0.172103; MAE 9.340692/9.558580); it is descriptive and **not independent validation**. Per-GSE and alternate integrations may reuse overlapping observations and are not additional independent cohorts.

## CellAge and mechanisms

CellAge integrates existing scores with a historical partial K=8,000 EN snapshot, not the final clocks. At the primary ≥20-cell threshold, adjusted identifiable hUSI–age tests comprise 15 fibroblast and 6 endothelial subtype models; none has BH q<0.05 in the archived tables. Unsupported or rank-deficient fits remain explicit. Historical “donor” column labels refer to configured sample aggregation and do not independently establish donor identity.

The inspected SenCID recommendation columns cannot support a claim about actual SenCID senescence intensity. hUSI is not measured in age years. Pathway enrichment and coefficient signs are exploratory associations with a specified feature universe and correction; they do not demonstrate causality or successful biological-age validation. No verified GWAS or transcription-factor network result is added from the conceptual overview image.

## Metric definitions

For paired chronological ages y and predictions p:

- R² = 1 − sum((y−p)²) / sum((y−mean(y))²).
- MAE = mean(abs(y−p)); RMSE = sqrt(mean((y−p)²)).
- MedAE = median(abs(y−p)).
- Source `MAD` and `Age_MAD` = mean(abs(y−mean(y))); they describe label dispersion.
- PCC = Pearson correlation; it does not establish calibrated prediction accuracy.

Undefined small-sample/constant-age metrics remain missing. `MAD` values should agree across models evaluated on identical age labels. In the historical Figure 2 table the column name is retained for fidelity; its interpretation is age dispersion.
