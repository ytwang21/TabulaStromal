# Results Part 3 — Stable expression features nominate lineage-associated molecular hypotheses

## Group-bootstrap stability and subtype support

The interpretable branch used 300 successful group bootstrap fits within each of five outer-training partitions and in the full development set, giving 1,800 checked records per lineage. Filtering, ranking and scaling were repeated within resamples. At selection frequency ≥0.70, full-development stability identified 438 fibroblast Celltype_Gene pairs representing 416 genes in five detailed cell types, and 23 endothelial pairs representing 23 genes in one cell type ([Fig. 3C,D,G–L](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.pdf)). Each held-out stable-model comparison used its corresponding outer-training stable set; full-development stability was used for exploratory explanation.

Fibroblast pair counts were eFibro_COL1A1_pan, 205 pairs; Myfibro_TAGLN_pan, 153 pairs; Myfibro_CALD1_pan, 75 pairs; eFibro_DCN_pan, 4 pairs; eFibro_COL1A2_pan, 1 pair. The endothelial set was confined to CapEC_TM4SF1_pan. Candidate-normalized fractions and source cell/sample coverage accompany the raw counts. Unequal observation and candidate-gene opportunity limit interpretation of these counts as intrinsic subtype importance.

## Stable-feature restriction did not improve age prediction

For fibroblasts, restricting EN to stable features reduced conditional mean validation R² from 0.742 to 0.636 and increased MAE from 11.86 to 14.65 years. For endothelial cells, R² decreased from 0.633 to 0.123 and MAE increased from 14.18 to 25.68 years (Fig. 3A,B). Neither stable version met the promotion criterion. The stable signatures therefore provide an interpretation of selected expression features, whereas the original-feature EN remains the final predictor. Selection stability does not itself imply superior prediction.

## Shared membership is not matched causal direction

Nine genes were shared between the two final stable sets: CD59, CYSTM1, DDAH2, NREP, RPS27L, SOD2, STMN1, WTAP, ZFP36 ([shared-gene supplement](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_FinalClock_shared_specific_genes.png)). Shared membership means recurrence under the recorded selection procedures; it does not establish equal coefficient direction or independent causal replication. Apparent lineage-restricted membership can also reflect differences in sample support and feature observability.

## Functional annotation separates matrix-associated and vesicle-associated signals

Global fibroblast enrichment included Collagen-Containing Extracellular Matrix (GO:0062023; 27/416 query genes; 11,337-gene background; BH FDR=6.06×10⁻⁷) and the Reactome term Extracellular Matrix Organization (24/416 query genes; the same background; FDR=6.59×10⁻⁶). Endothelial enrichment included ER to Golgi Transport Vesicle Membrane (GO:0012507; 3/23 query genes; 10,949-gene background; FDR=0.0133; Fig. 3O,P). These are observed overlaps of stable genes with annotated gene sets, not demonstrations of pathway activation or causal aging mechanisms.

Across the specified global GO, Reactome and KEGG libraries, 56 fibroblast and 12 endothelial terms met FDR<0.05. Adding the separately retained Hallmark library increased the fibroblast count to 75. The exported global tables retain zero-overlap tested terms and the original within-scope/library BH values; the multiplicity family was not recomputed after filtering. Overlapping annotations driven by the same genes do not represent independent mechanistic discoveries.

The eight supplied aging lists remain provided_unverified. Their measured overlaps in Fig. 3Q,R are statistical comparisons with those supplied lists, not validation against authenticated official signatures.

## Expression–age associations define follow-up hypotheses

The saved covariate-adjusted expression–age analyses contained FDR<0.05 results for 10 of 438 fibroblast pairs and 4 of 23 endothelial pairs. These associations remain distinct from model-coefficient signs. Per-study direction resources can share samples and lack a prespecified universal replication criterion, so they do not independently validate every stable feature. The older CellAge signature branch is also distinct from these final stable sets. No perturbation experiment, authenticated GWAS result or transcription-factor network analysis establishes a molecular driver here. The supported interpretation is a set of predictive features and molecular hypotheses for further testing.

The source-to-claim mapping is available in [the evidence registry](../../evidence/evidence_registry.csv) and [the paragraph evidence report](../../audit/claim_to_evidence_report.csv).
