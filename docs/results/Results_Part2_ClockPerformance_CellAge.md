# Results Part 2 — Celltype-resolved age clocks distinguish internal prediction, external transfer and cellular senescence

## Clock overview and complete grouped benchmark

TabulaStromal predicts recorded chronological age from sample-level Celltype_Gene expression features ([Fig. 2A](../../F2_ModelPerformance/figures/TabulaStromal_Workflow_V6.png)). Expression was aggregated within sample and annotated cell type; explicit cell-type proportions or ratios were excluded. Feature filtering, ranking and scaling were fitted within training partitions. The local DeepSets, Pseudobulk and Ensemble implementations are IDEAL-Age-inspired estimators and are not asserted to reproduce an external method protocol exactly.

The development cohorts contained 640 fibroblast samples in 634 conservative GSE/GSM proxy groups and 437 endothelial samples in 431 groups. Thirteen model families were evaluated at K=8,000, 10,000, 15,000 and 20,000 using grouped five-fold CV. Saved groups were disjoint within every training/validation split. This check cannot resolve unannotated donor reuse across studies.

## Internal performance and source-figure scope

For the fibroblast collection, XGBoost at K=20,000 ranked first by mean validation R², with R²=0.803 ± 0.020, MAE=9.71 ± 0.66 years and PCC=0.899 ± 0.012 (mean ± sample SD across folds).

For the endothelial collection, ExtraTrees at K=20,000 ranked first by mean validation R², with R²=0.716 ± 0.116, MAE=11.54 ± 1.84 years and PCC=0.848 ± 0.074 (mean ± sample SD across folds).

This is a descriptive ranking across the completed benchmark, not a model-superiority significance test. The supplied [Fig. 2B–F](../../F2_ModelPerformance/figures/Fig2_StromalClock_ModelPerformance.pdf) uses source fold 4, with each model's K chosen from that fold's validation result. It is not a five-fold mean. Its external values use FT30 followed by ALL100 evaluation, including fitted adaptation samples. MAD is the mean absolute deviation of chronological-age labels from their mean; it describes label dispersion rather than prediction error. The figure's biological-age axis wording refers to prediction of the recorded chronological-age target.

## A separate final interpretable ElasticNet branch

The fibroblast EN branch selected K=20,000 and achieved conditional validation R²=0.742 ± 0.031 and MAE=11.86 ± 0.95 years. The endothelial EN branch selected K=10,000 and achieved conditional validation R²=0.633 ± 0.124 and MAE=14.18 ± 1.98 years. Both clocks retained the original selected feature representation because stable-feature restriction predicted less well. They were refitted to all development samples ([Fig. 3A,B](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.pdf)). Selection of the EN family and its K is separate from the global benchmark ranking. Because K was selected using the same internal comparison, these CV estimates are conditional on that selection rather than fully nested estimates of the complete selection procedure.

## External transfer depends on adaptation and the evaluation population

Without adaptation, the full-development fibroblast EN yielded R²=−2.943 and MAE=29.32 years in 93 merged external samples. The endothelial EN yielded R²=−1.720 and MAE=24.89 years in 77 samples (Fig. 3E,F). Negative R² preserves an important transfer limitation: the predictions had greater squared error than the external cohort-mean reference.

For the final EN adaptation protocol, 28 fibroblast and 23 endothelial samples were used for regularized refitting with frozen features, scaler and source hyperparameters, leaving 65 and 54 group-disjoint samples for evaluation. On the same held-out fibroblast samples, MAE decreased from 29.02 to 13.14 years and R² changed from −3.168 to −0.001. For endothelial samples, MAE decreased from 25.65 to 13.18 years and R² changed from −2.194 to −0.355. Adaptation reduced error, while held-out R² remained near zero or negative (Fig. 3M,N).

| Lineage | Evaluation population | n | R² | MAE (years) |
| --- | --- | ---: | ---: | ---: |
| Fibro | No adaptation, all external | 93 | -2.943 | 29.32 |
| Fibro | No adaptation, held-out 70% | 65 | -3.168 | 29.02 |
| Fibro | After adaptation, held-out 70% | 65 | -0.001 | 13.14 |
| Fibro | After adaptation, ALL100 (mixed) | 93 | 0.358 | 9.34 |
| Endo | No adaptation, all external | 77 | -1.720 | 24.89 |
| Endo | No adaptation, held-out 70% | 54 | -2.194 | 25.65 |
| Endo | After adaptation, held-out 70% | 54 | -0.355 | 13.18 |
| Endo | After adaptation, ALL100 (mixed) | 77 | 0.172 | 9.56 |

Pooling the fitted and held-out observations gives the more favorable ALL100 metrics shown in the last row for each lineage. Those pooled metrics are mixed apparent/unseen performance and are not independent validation. Primary no-adaptation performance, matched held-out adaptation comparisons and the historical source-fold figure therefore remain distinct results. No independent health-outcome analysis establishes these predictions as biological-age measurements.

## CellAge identifies subtype differences without a confirmed adjusted age association

The audited CellAge snapshot included hUSI outputs for 1,665,482 fibroblast-lineage and 536,017 endothelial-lineage cells; recSID recommendations covered 1,372,858 and 352,228 cells, respectively ([Fig. 4](../../CellAge_Exploration/figures/Fig_CellAge_Stromal_Integration.pdf)). hUSI was analyzed as a senescence-associated score, not an age in years or a cross-lineage calibrated quantity. The inspected SenCID_rec_SID columns contain identity recommendations and do not provide genuine senescence-intensity measurements.

iFibro_CD74_pan had a positive paired hUSI contrast against the equal-weight mean of other observed subtypes within the same sample/tissue (effect 0.0539, 95% bootstrap CI 0.0427–0.0655; BH q=8.9e-12; 62 source-designated sampling units). EndoMT_PTN_pan had a positive paired hUSI contrast against the equal-weight mean of other observed subtypes within the same sample/tissue (effect 0.0826, 95% bootstrap CI 0.0630–0.1012; BH q=3.3e-11; 67 source-designated sampling units). These differences describe subtype-associated scores; they are not slopes with chronological age. The original source calls these units donors under its sample=donor assumption, which the available cross-study metadata cannot independently verify.

At the primary threshold of at least 20 scored cells per sample/tissue/subtype, 15 fibroblast and 6 endothelial adjusted hUSI–age models were identifiable; none met BH q<0.05 (Fig. 4E,F). The historical K=8,000 EN integration did not yield an identifiable fully adjusted hUSI–delta-age association. That partial historical prediction snapshot and its older stable features differ from the final EN clocks. A new joint analysis using complete final-clock OOF predictions and genuine SenCID intensities remains unavailable; it is not reported as a completed result.

The source-to-claim mapping is available in [the evidence registry](../../evidence/evidence_registry.csv) and [the paragraph evidence report](../../audit/claim_to_evidence_report.csv).
