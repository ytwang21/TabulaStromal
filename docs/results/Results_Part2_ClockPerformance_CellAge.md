# Results Part 2. Celltype-resolved clocks distinguish internal age information from external performance

## Clock overview and grouped model comparison

TabulaStromal predicts recorded chronological age from expression aggregated within each sample and annotated cell state, producing Celltype_Gene features (Fig. 2A). Explicit cell-type proportions or ratios were excluded. The development sets contained 640 fibroblast samples assigned to 634 conservative GSE/GSM proxy groups and 437 endothelial samples assigned to 431 groups. Grouped five-fold cross-validation compared thirteen model families at requested feature counts of 8,000, 10,000, 15,000 and 20,000. Each lineage had all 260 model–feature-count–fold combinations. Proxy groups were disjoint within saved training/validation splits, although unresolved cross-study donor reuse could not be excluded.

*Source record: [C0007](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Across the full benchmark, fibroblast XGBoost at 20,000 features ranked first by mean validation R², achieving 0.803 ± 0.020 with MAE 9.71 ± 0.66 years. Endothelial ExtraTrees at 20,000 features ranked first with R² 0.716 ± 0.116 and MAE 11.54 ± 1.84 years (mean ± sample standard deviation across five folds). These are descriptive model rankings; no paired superiority test established that either configuration outperformed every alternative. The three IDEAL-Age-inspired implementations were local variants, and their results do not constitute reproduction of the published immune-aging model [9].

*Source record: [C0008](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The heatmaps in Fig. 2B–D present a different summary: source fold 4, with feature count selected for each model using that fold's validation result. They contain training, validation and FT30→ALL100 external values and are neither five-fold means nor a benchmark of the later compact models. Moreover, MAD in panel D is the mean absolute deviation of observed ages from their mean. It describes age-label dispersion rather than prediction error.

*Source record: [C0009](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Current performance-first ElasticNet fits

The current Fig. 2E,F and the model-fit elements in Fig. 2A use the subsequent performance-first ElasticNet source models with 750 fibroblast and 350 endothelial Celltype_Gene features. Applied to the development data used for fitting, these models achieved R² 0.950 and 0.950, MAE 4.98 and 5.03 years, and median absolute error 3.66 and 3.50 years, respectively (Table 1). Each plotted point represents one development sample; the displayed 640 and 437 observations are in-sample predictions. These panels do not establish out-of-fold performance of the compact feature sets.

*Source record: [C0010](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The adapted performance-first models were evaluated on 93 fibroblast and 77 endothelial external samples, yielding R² 0.377 and 0.391, MAE 8.71 and 9.28 years, and median absolute error 4.69 and 6.42 years. These are the values in the lower external panels of Fig. 2A. The ALL100 designation identifies all eligible external observations rather than a literal sample count of one hundred. External outcomes participated in selecting the feature representation and optimization settings, in addition to supplying labeled adaptation samples. Accordingly, these results are outcome-informed external-cohort performance, not independent external validation.

*Source record: [C0011](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Table 1. Current performance-first models. Values describe distinct evaluation populations and must not be compared as equivalent validation estimates.

| Lineage | Evaluation population | Features | n | R² | MAE, years | MedAE, years |
|---|---|---:|---:|---:|---:|---:|
| Fibroblast | Development, in-sample source-model fit | 750 | 640 | 0.950 | 4.98 | 3.66 |
| Endothelial | Development, in-sample source-model fit | 350 | 437 | 0.950 | 5.03 | 3.50 |
| Fibroblast | Adapted, ALL100 used in model selection | 750 | 93 | 0.377 | 8.71 | 4.69 |
| Endothelial | Adapted, ALL100 used in model selection | 350 | 77 | 0.391 | 9.28 | 6.42 |
| Fibroblast | Adapted, reused 70% subset | 750 | 65 | 0.032 | 12.02 | 8.47 |
| Endothelial | Adapted, reused 70% subset | 350 | 54 | 0.043 | 11.67 | 9.91 |

*Source record: [C0031](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Even the reused 70% subset was involved in ALL100-based model selection. Its R² of 0.032 for fibroblasts and 0.043 for endothelial cells therefore remains selection-affected. In contrast to panel D, MAD in the current model-fit and external panels denotes median absolute prediction error, here reported as MedAE. This difference in metric definitions is retained explicitly rather than assigning one meaning to all MAD labels.

*Source record: [C0012](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Protocol-specific evidence for the transfer gap

An earlier unreduced ElasticNet branch provides a separate comparison of transfer before and after adaptation. It used 20,000 fibroblast and 10,000 endothelial features, selected after the family-specific benchmark, and achieved conditional five-fold validation R² of 0.742 and 0.633 with MAE 11.86 and 14.18 years. Without adaptation, the full-development models yielded external R² −2.943 and −1.720 and MAE 29.32 and 24.89 years in 93 and 77 samples. Negative R² indicates worse squared-error prediction than the mean of the corresponding external age labels.

*Source record: [C0013](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

For that earlier branch, adaptation used 28 fibroblast and 23 endothelial samples, leaving 65 and 54 group-disjoint samples for evaluation. On matched held-out observations, MAE fell from 29.02 to 13.14 years for fibroblasts and from 25.65 to 13.18 years for endothelial cells, while adapted R² remained −0.001 and −0.355. These results document a substantial transfer gap despite reduced error after target-cohort refitting. They are not the validation results of the current 750/350-feature models. Alternate integrations of the same external samples and per-GSE representations are likewise not independent biological replications.

*Source record: [C0014](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Cellular senescence scores answer a distinct question

The supporting CellAge snapshot contained hUSI scores for 1,665,482 fibroblast-lineage and 536,017 endothelial-lineage cells, with SenCID identity recommendations available for 1,372,858 and 352,228 cells. hUSI represents senescence-associated transcription rather than age in years [7]. The available recSID dimensions characterize senescence identities and cannot substitute for genuine SenCID senescence-intensity outputs [8]. Dataset-wise score normalization also limits interpretation of absolute cross-dataset or cross-lineage differences.

*Source record: [C0015](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The historical snapshot used an older, partially evaluated clock branch, so its clock–score associations cannot validate the current performance-first models. We therefore aligned the current selected cell states to the adult ModelFeature_CellAge analysis reported in Results Part 3 and Fig. 6. That analysis distinguishes within-sample relative hUSI from within-state age associations and composition changes. A complete evaluation linking current-model out-of-fold predictions or independently validated age acceleration to senescence scores remains unavailable.

*Source record: [C0016](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).
