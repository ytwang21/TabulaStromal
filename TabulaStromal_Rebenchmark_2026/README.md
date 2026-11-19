# Rebenchmark — Current compact model and earlier EN protocols

The current performance-first representation retains **750/350 Celltype_Gene pairs**, **678/337 unique genes**, and **27 states per lineage**. Feature masks and optimization settings were selected using ALL100 outcomes. Development-set R²≈0.950 is in-sample; adapted ALL100 R²=0.377/0.391 is outcome-selected, not independent validation. The reused 70% subset also participated in selection.

[Current selection summary](tables/current_performance_first_summary.csv) · [Fibroblast features/protocol](tables/Fibro/current_performance_first/) · [Endothelial features/protocol](tables/Endo/current_performance_first/) · [Current downstream interpretation](../F3_ds_analysis/) · [Results Part 2](../docs/results/Results_Part2_ClockPerformance_CellAge.md)

![Current compact EN source composite](figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.png)

*The source composite now depicts the 750/350-feature models. Its nested-GSE display is retrospective and is not a new independent cohort. The reused filename must not be read using the older original-EN panel-source CSV.*

## Earlier unreduced ElasticNet branch

The preserved `final_clock_*` filenames and figures below describe the earlier original-EN branch (K=20,000/10,000), not the compact 750/350-feature model. Its conditional internal R² was 0.742/0.633; unadapted external R² was −2.943/−1.720, with MAE 29.32/24.89 years. Adaptation on 28/23 samples reduced matched held-out MAE to 13.14/13.18 in 65/54 observations, while R² remained −0.001/−0.355. Its 300-bootstrap, frequency-0.70 stable sets contained 438/23 pairs and are not the current F3 feature sets.

[Earlier tables](tables/) · [Earlier figure archive](figures/) · [Evidence](evidence/evidence_registry.csv) · [Registered assets](audit/released_assets.csv)
