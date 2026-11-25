# Results Part 2. Celltype-resolved clocks distinguish internal age information from external performance

## Clock overview and grouped model comparison

TabulaStromal predicts recorded chronological age from expression aggregated within each sample and annotated cell state, producing Celltype_Gene features (Fig. 2A). Explicit cell-type proportions or ratios were excluded. The development sets contained 640 fibroblast samples assigned to 634 conservative GSE/GSM proxy groups and 437 endothelial samples assigned to 431 groups. Grouped five-fold cross-validation compared fourteen model families at requested feature counts of 8,000, 10,000, 15,000 and 20,000. Each lineage had all 280 model–feature-count–fold combinations. Proxy groups were disjoint within saved training/validation splits, although unresolved cross-study donor reuse could not be excluded.

*Source record: [C0012](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Across the full benchmark, fibroblast XGBoost at 20,000 features ranked first by mean validation R², achieving 0.803 ± 0.019 with MAE 9.81 ± 0.61 years. Endothelial XGBoost at 20,000 features ranked first with R² 0.727 ± 0.092 and MAE 11.73 ± 1.69 years (mean ± sample standard deviation across five folds). These are descriptive model rankings; no paired superiority test established that either configuration outperformed every alternative. The three IDEAL-Age-inspired implementations were local variants, and their results do not constitute reproduction of the published immune-aging model [9].

*Source record: [C0013](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The heatmaps in Fig. 2B–D present a different summary: source fold 4, with feature count selected for each model using that fold's validation result. They contain training, validation and FT30→ALL100 external values and are neither five-fold means nor a benchmark of the later compact models. Moreover, MAD in panel D is the mean absolute deviation of observed ages from their mean. It describes age-label dispersion rather than prediction error.

*Source record: [C0014](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Current performance-first ElasticNet fits

The current Fig. 2E,F and the model-fit elements in Fig. 2A use the subsequent performance-first ElasticNet source models with 750 fibroblast and 350 endothelial Celltype_Gene features. Applied to the development data used for fitting, these models achieved R² 0.950 and 0.950, MAE 4.98 and 5.03 years, and median absolute error 3.66 and 3.50 years, respectively (Table 1). Each plotted point represents one development sample; the displayed 640 and 437 observations are in-sample predictions. These panels do not establish out-of-fold performance of the compact feature sets.

*Source record: [C0015](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The adapted performance-first models were evaluated on 93 fibroblast and 77 endothelial external samples, yielding R² 0.377 and 0.391, MAE 8.71 and 9.28 years, and median absolute error 4.69 and 6.42 years. These are the values in the lower external panels of Fig. 2A. The ALL100 designation identifies all eligible external observations rather than a literal sample count of one hundred. External outcomes participated in selecting the feature representation and optimization settings, in addition to supplying labeled adaptation samples. Accordingly, these results are outcome-informed external-cohort performance, not independent external validation.

*Source record: [C0016](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Table 1. Current performance-first models. Values describe distinct evaluation populations and must not be compared as equivalent validation estimates.

| Lineage | Evaluation population | Features | n | R² | MAE, years | MedAE, years |
|---|---|---:|---:|---:|---:|---:|
| Fibroblast | Development, in-sample source-model fit | 750 | 640 | 0.950 | 4.98 | 3.66 |
| Endothelial | Development, in-sample source-model fit | 350 | 437 | 0.950 | 5.03 | 3.50 |
| Fibroblast | Adapted, ALL100 used in model selection | 750 | 93 | 0.377 | 8.71 | 4.69 |
| Endothelial | Adapted, ALL100 used in model selection | 350 | 77 | 0.391 | 9.28 | 6.42 |
| Fibroblast | Adapted, reused 70% subset | 750 | 65 | 0.032 | 12.02 | 8.47 |
| Endothelial | Adapted, reused 70% subset | 350 | 54 | 0.043 | 11.67 | 9.91 |

*Source record: [C0046](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Even the reused 70% subset was involved in ALL100-based model selection. Its R² of 0.032 for fibroblasts and 0.043 for endothelial cells therefore remains selection-affected. In contrast to panel D, MAD in the current model-fit and external panels denotes median absolute prediction error, here reported as MedAE. This difference in metric definitions is retained explicitly rather than assigning one meaning to all MAD labels.

*Source record: [C0017](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Protocol-specific evidence for the transfer gap

An earlier unreduced ElasticNet branch provides a separate comparison of transfer before and after adaptation. It used 20,000 fibroblast and 10,000 endothelial features, selected after the family-specific benchmark, and achieved conditional five-fold validation R² of 0.742 and 0.633 with MAE 11.86 and 14.18 years. Without adaptation, the full-development models yielded external R² −2.943 and −1.720 and MAE 29.32 and 24.89 years in 93 and 77 samples. Negative R² indicates worse squared-error prediction than the mean of the corresponding external age labels.

*Source record: [C0018](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

For that earlier branch, adaptation used 28 fibroblast and 23 endothelial samples, leaving 65 and 54 group-disjoint samples for evaluation. On matched held-out observations, MAE fell from 29.02 to 13.14 years for fibroblasts and from 25.65 to 13.18 years for endothelial cells, while adapted R² remained −0.001 and −0.355. These results document a substantial transfer gap despite reduced error after target-cohort refitting. They are not the validation results of the current 750/350-feature models. Alternate integrations of the same external samples and per-GSE representations are likewise not independent biological replications.

*Source record: [C0019](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Cellular senescence scores answer a distinct question

The supporting CellAge snapshot contained hUSI scores for 1,665,482 fibroblast-lineage and 536,017 endothelial-lineage cells, with SenCID identity recommendations available for 1,372,858 and 352,228 cells. hUSI represents senescence-associated transcription rather than age in years [7]. The available recSID dimensions characterize senescence identities and cannot substitute for genuine SenCID senescence-intensity outputs [8]. Dataset-wise score normalization also limits interpretation of absolute cross-dataset or cross-lineage differences.

*Source record: [C0020](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The historical snapshot used an older, partially evaluated clock branch, so its clock–score associations cannot validate the current performance-first models. We therefore aligned the current selected cell states to the adult ModelFeature_CellAge analysis reported in Results Part 3 and Fig. 6. That analysis distinguishes within-sample relative hUSI from within-state age associations and composition changes. A complete evaluation linking current-model out-of-fold predictions or independently validated age acceleration to senescence scores remains unavailable. The subsequent annotation-transfer update does not retrospectively change these saved prediction cohorts.

*Source record: [C0021](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Exploratory clock–score alignment does not establish age acceleration

A subsequent analysis matched frozen-model development predictions or cell-state coefficient contributions to hUSI, retaining 99 candidate associations, of which 46 were estimable (Extended Data Fig. 6C,D). Three fibroblast associations passed correction within their prespecified scope-by-metric families across both lineages. The Myfibro_CALD1_pan contribution was negatively associated with hUSI (effect per ten contribution-years, −0.0539; 95% CI, −0.0795 to −0.0282; q = 0.00562; 73 sample proxies across 23 studies). Myfibro_CALD1_kidney showed a positive development clock-error association and a negative contribution association. No endothelial association survived correction. Adjustment for actual age and source covariates does not remove dependence introduced by in-sample prediction, selection and shared transcriptional inputs. These results are exploratory clock–score alignments, not validated biological-age acceleration or evidence that a cell state drives aging.

*Source record: [C0022](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Updated process-based comparator remains separate from the compact clock

The fourteen-family benchmark additionally included JungClock-TabulaStromal-adapted, a process-based native-H2O GLM adaptation of the published fibroblast transcriptomic framework [22]. Its best mean validation R² was 0.582 ± 0.100 in fibroblasts at requested K = 20,000, with MAE 12.85 years, and 0.497 ± 0.159 in endothelial cells at K = 8,000, with MAE 14.50 years. Each summary included all five folds. Requested K bounded candidate Celltype_Gene pairs rather than the final number retained after process selection. This local adaptation used training-fold process selection and normalized cell-state expression, omitted bulk TPM/fsva processing and retained a native H2O backend; it was not a numerical reproduction of the original bulk study. The saved Figure 2 B–D heatmaps contain thirteen models and omit this newly added comparator. The full fourteen-model benchmark and that figure snapshot therefore represent distinct result versions.

*Source record: [C0023](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).
