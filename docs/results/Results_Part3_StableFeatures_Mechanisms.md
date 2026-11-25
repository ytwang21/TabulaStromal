# Results Part 3. Selected features nominate context-dependent stromal molecular programs

## Performance-first features span multiple stromal states

Downstream analyses used the current performance-first representation: 750 fibroblast Celltype_Gene pairs encoding 678 unique genes and 350 endothelial pairs encoding 337 genes, spanning 27 detailed states per lineage (Fig. 3A,B). Selected features combined global and celltype-balanced candidate sources with recorded bootstrap denominators of 1,500. Their minimum retained selection frequencies were 0.600 and 0.602. These sets differ from the earlier 300-bootstrap, frequency-0.70 interpretation branch and were not interchanged with it. The current gene sets shared 45 genes, leaving 633 fibroblast-only and 292 endothelial-only memberships; these membership differences do not establish intrinsic lineage specificity.

*Source record: [C0024](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Feature counts were uneven across states. eFibro_COL1A1_pan contributed 70 fibroblast pairs, followed by Myfibro_CALD1_pan with 65 and qFibro_CLU_pan with 62. VenEC_CLU_eye contributed 64 endothelial pairs, followed by EndoMT_VIM_pan with 37 and ArtEC_FLT1_pan with 35. Counts describe representation within the selected model, not cell abundance, independently tested subtype importance or evidence of a lineage transition. Because external outcomes influenced model selection, these downstream annotations are exploratory descriptions of the selected predictors.

*Source record: [C0025](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Expression-pattern displays retained 520 fibroblast and 161 endothelial features after coverage and dynamic-range filtering (Fig. 3C,D). Of the fibroblast features, 183, 43 and 294 had their highest display-window mean in the source-defined young, middle and old windows; the corresponding endothelial counts were 34, 21 and 106. These windows group original A–H labels for visualization and do not establish equal-width chronological-age intervals. Display filtering and row scaling are not differential-expression tests. Likewise, a positive or negative ElasticNet coefficient is a conditional predictive weight, not a demonstrated pro-aging or anti-aging effect.

*Source record: [C0026](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Display-window patterns show uneven reliability and shared functions

The updated pattern analysis quantified uncertainty in the 520 fibroblast and 161 endothelial displayed features rather than treating all window assignments as equally reliable (Extended Data Fig. 5A–D). Bootstrap assignment support was at least 0.8 for 445 fibroblast and 137 endothelial features. However, 149 and 80 features had at least one observed age group represented by only one sample proxy. Requiring bootstrap support, contributing-study deletion support and adequate per-group coverage retained conservative subsets of 303 and sixteen features. Alternative window boundaries preserved 490 and 142 assignments, whereas excluding the earliest A group preserved 444 and 152. Because the source A group includes childhood observations, early-window expression cannot be interpreted specifically as an adult anti-aging program.

*Source record: [C0027](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Functional interpretation depended on the comparison universe. Measured-gene backgrounds identified annotations such as MYC targets in the fibroblast young-high set and metal-response terms in the endothelial old-high set. In contrast, no term in either lineage met the specified significance and overlap criteria when the background was restricted to displayed selected genes, using complete library-by-three-pattern correction families. Thus, the available enrichments do not establish functions preferentially associated with one pattern beyond those shared by the displayed selected representation. Coefficient signs and pattern composition also describe model features rather than age-dependent changes in cell abundance.

*Source record: [C0028](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Twenty-two unique genes were shared between the two lineages' displayed feature sets, compared with 45 shared genes across the complete selected sets (Extended Data Fig. 5E). Seventeen shared displayed genes had matching pattern memberships and five differed. Four genes belonged to the young-high set in both lineages, exceeding the gene-identity-shuffle expectation of 1.0754 (nine-comparison BH q = 0.0300). The thirteen old-high/old-high memberships did not survive that correction (q = 0.0639). Permutation preserved multi-pattern memberships, so sums over the overlap matrix are not counts of unique genes. These results identify limited shared pattern structure without establishing a formal lineage-by-age interaction.

*Source record: [C0029](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Functional annotation distinguishes lineage-associated themes

Functional analysis used measured lineage-specific backgrounds of 11,700 fibroblast and 11,521 endothelial genes and complete library-specific correction families (Fig. 3E–J). In fibroblasts, all selected genes were enriched for the Hallmark inflammatory response (18 overlapping genes; BH q = 0.0185), TNFα signaling via NF-κB (20 genes; q = 0.0339) and the collagen-containing extracellular matrix cellular-component term (26 genes; q = 0.0269). These overlaps connect predictive features to inflammatory and matrix-associated annotations without demonstrating pathway activation or its cause.

*Source record: [C0030](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The endothelial set was enriched for cytoskeleton (30 genes; q = 0.0194), KEGG mineral absorption (seven genes; q = 0.00791) and the Reactome 2022 metallothionein metal-binding term (four genes; q = 0.0183). The latter overlap comprised MT1E, MT1M, MT1X and MT2A and also contributed to related metal-response annotations. Such overlapping terms are not independent replications and do not establish exposure to a particular metal. The current main network uses Reactome 2022; separate older signed-feature outputs using Reactome 2024 were not merged into this correction family.

*Source record: [C0031](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Gene–theme networks summarized significant memberships rather than inferred causal interactions. A separate KEGG context display showed 28 fibroblast pathways covering 132 genes and 25 endothelial pathways covering 52 genes (Extended Data Fig. 3). Inclusion in that display did not require every row to be significant; the fibroblast global all-gene KEGG analysis had no FDR-significant term. Theme membership and selected network edges therefore serve as annotation aids rather than additional statistical discoveries.

*Source record: [C0032](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Regulator inference identifies candidates with limited adult age support

LISA2 analysis integrated public chromatin information with selected gene sets for 14 of 27 fibroblast states and five of 27 endothelial states (Fig. 4A–C). Other states lacked sufficient mapped input genes and were retained as unavailable. After accounting for the number of ChIP datasets per factor and correcting across all factor–state tests within each lineage, examples included TP53 in Myfibro_TUBB2B_pan (q = 7.51 × 10⁻⁷), SOX4 in eFibro_COL1A1_pan (q = 8.44 × 10⁻⁸) and TEAD1 in Myfibro_TAGLN_pan (q = 0.0123). Endothelial candidates included HSF1 in EndoMT_PTN_pan (q = 1.78 × 10⁻¹²), SMARCC1 in ArtEC_FLT1_brain (q = 0.000146) and RUNX1T1 in EndoMT_VIM_pan (q = 0.00186). These results indicate regulatory association with the input gene set; they do not establish factor activation with age or direct binding in the sampled stromal cells.

*Source record: [C0033](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

A separate CollecTRI-prior-constrained RNA analysis was evaluable in 21 fibroblast and 15 endothelial states. It tested 352,675 and 241,986 candidate edges and retained 13,994 and 4,747 supported edges, including 65 and seven edges whose targets matched selected Celltype_Gene features (Fig. 4D,E). These are edge counts rather than unique target-gene counts. Their interpretation depends on prior knowledge, covariate adjustment, sign concordance and resampling support. A stricter adult within-cohort sensitivity analysis tested 1,434 fibroblast and 1,054 endothelial factor–state associations; none survived lineage-wide BH correction. Thus, the network candidates were not established as general adult age-associated regulators.

*Source record: [C0034](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Adult raw-count analysis reveals tissue-dependent matrix-associated contrasts

To examine matrix-associated transcription in adults, raw counts were aggregated within training-sample proxies and selected states after age-concordance and coverage filtering. The analysis retained 300,431 fibroblast cells in 1,748 proxy-by-state pseudobulks and 166,664 endothelial cells in 905 pseudobulks. Counts with fractional entries were not rounded; noninteger genes were excluded separately for each fitted profile. Negative-binomial aged-versus-young models were estimable for 39 cell-state or cell-state-by-tissue profiles (Fig. 5A–D). AGRN was the only shared selected ECM-associated gene, a membership observation rather than a replicated age effect.

*Source record: [C0035](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Using full-gene Wald statistics as inputs to signed gene-set regression, 67 of 312 ECM-signature comparisons met FDR ≤ 0.05 after correction across both lineages: 41 fibroblast and 26 endothelial comparisons. These include overlapping signatures and analysis scopes and are not 67 independent biological replications. For eFibro_COL1A1_pan, the collagen-formation signature had a negative contrast in ovary (ULM score −5.81; q = 1.03 × 10⁻⁶; eight sample proxies), whereas ECM organization had a positive contrast in colon (+4.65; q = 0.000116; eleven proxies). ArtEC_FLT1_pan in adipose showed a positive collagen-formation contrast (+5.33; q = 7.98 × 10⁻⁶; twelve proxies). Each example arose from one study, and these separate contrasts are not a formal age-by-tissue interaction test.

*Source record: [C0036](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Among 142 displayed mechanosensitive-regulator comparisons, only TEAD4 in adipose ArtEC_FLT1_pan survived correction across the full 22,346-test transcription-factor family (ULM score +3.95; q = 0.0174; seven younger and five older proxies from one study; Fig. 5E). The results therefore do not support a general claim of YAP/TAZ–TEAD activation across stromal states. At the exact selected-feature level, PLA2G4A in ovary eFibro_COL1A1_pan was the sole significant ECM-associated example (aged-versus-young log₂ fold change 2.60; cross-profile q = 0.0132; Table 2). No endothelial exact-feature association passed that correction. No same-state candidate simultaneously satisfied the LISA2, TF-age and target-gene-age significance criteria in the integrated candidate table. The broader heatmaps in Fig. 5F,G instead summarize globally selected ECM-associated genes across cell-state-level profiles and do not establish exact-feature support in every displayed state.

*Source record: [C0037](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Table 2. Selected adult contrasts. ULM uses full-gene Wald statistics; log₂FC is an individual-gene negative-binomial age effect. All five examples are supported by one study each. ECM and TF q-values are corrected across both lineages in separate complete families; the gene q-value uses all tested gene-by-profile comparisons within its lineage.

| Context | Readout | Statistic | BH q | n proxies | Studies |
|---|---|---:|---:|---:|---:|
| Fibroblast: eFibro_COL1A1_pan, ovary | Collagen formation | ULM −5.81 | 1.03 × 10⁻⁶ | 8 | 1 |
| Fibroblast: eFibro_COL1A1_pan, colon | ECM organization | ULM +4.65 | 0.000116 | 11 | 1 |
| Endothelial: ArtEC_FLT1_pan, adipose | Collagen formation | ULM +5.33 | 7.98 × 10⁻⁶ | 12 | 1 |
| Endothelial: ArtEC_FLT1_pan, adipose | TEAD4 regulon | ULM +3.95 | 0.0174 | 12 | 1 |
| Fibroblast: eFibro_COL1A1_pan, ovary | PLA2G4A | log₂FC +2.60 | 0.0132 | 8 | 1 |

*Source record: [C0047](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Selected cellular states differ in relative senescence-associated transcription

The current ModelFeature_CellAge analysis aligned the 27 selected states per lineage to adult, age-concordant samples, encompassing 731 fibroblast sample proxies from 87 studies and 529 endothelial proxies from 80 studies before state-specific primary coverage requirements (Fig. 6). The primary within-sample comparison estimated each state's mean hUSI difference from the cell-weighted mean of other scored cells in the same lineage and sample-by-tissue unit. Seventeen fibroblast and fourteen endothelial states had estimable relative-score contrasts; seven and eight, respectively, survived correction across both lineages within this metric.

*Source record: [C0038](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Fibroblast eFibro_COL1A2_pan showed a positive relative hUSI difference of 0.0315 (95% CI, 0.0216–0.0414; q = 1.29 × 10⁻⁶; 185 sample proxies across 45 studies), as did Myfibro_CALD1_pan (0.0287; 0.0186–0.0388; q = 1.46 × 10⁻⁵; 134 proxies across 34 studies). Endothelial VenEC_CLU_pan showed a negative difference (−0.0166; −0.0237 to −0.00955; q = 0.000218; 125 proxies across 34 studies). These three directions remained significant in healthy-only, training-only, outside-training and alternate cell-threshold sensitivity analyses. Those analyses reuse overlapping observations and are sensitivity checks rather than independent replication.

*Source record: [C0039](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

No primary adjusted within-state hUSI–age, younger-versus-older hUSI, SSE-fraction or composition–age comparison survived its metric-specific correction. In particular, continuous hUSI–age models were estimable for sixteen fibroblast and thirteen endothelial states, with no FDR-significant slope. Model feature count and coefficient-derived importance also showed no FDR-significant exploratory correlation with the tested CellAge summaries. Moreover, 553 of 678 fibroblast selected genes and 276 of 337 endothelial genes overlapped the hUSI signature. Relative hUSI therefore describes a dependent transcriptomic state contrast, not independent validation of the clock or evidence that a selected cell state ages faster.

*Source record: [C0040](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Nonlinear age associations prioritize a subset of selected features

Adult analyses tested all 1,100 selected Celltype_Gene features using restricted cubic splines and source-aware covariates, with 649 estimable models. Joint tests of nonlinear components identified 119 features after correction across both lineages: 85 fibroblast and 34 endothelial features. Of these, 28 and ten, respectively, also showed positive held-out-study relative-shape prediction improvement (Fig. 7A,B). For example, eFibro_DCN_heart_PTGER2 had q = 8.77 × 10⁻¹¹ and relative MSE improvement 0.826 in thirty sample proxies from nine studies, but its fitted maximum was at the upper observed age boundary. ArtEC_THY1_pan_TAF15 had an interior fitted maximum near 60.35 years, although its relative MSE improvement was only 0.00163. Neither a boundary maximum nor a small positive improvement confirms a biological transition age. The analysis remains conditional on features selected by the original clock workflow.

*Source record: [C0041](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

A separate formal tissue-interaction analysis tested 215 state-by-module candidates, of which ten had adequate coverage and identifiable joint contrasts. Two eFibro_COL1A1_pan comparisons passed the complete cross-lineage correction: the HAECA ECM/glycoprotein module (q = 0.000176) and Reactome ECM organization (q = 0.00243), each supported by forty sample proxies across eight studies (Extended Data Fig. 6E). These tests address age-by-tissue heterogeneity directly and are distinct from the earlier single-study aged-versus-young contrasts. They used mean log-CPM module expression, rather than Wald-statistic-derived ULM activity, and do not measure ECM stiffness, deposition or a universal tissue-independent aging program.

*Source record: [C0042](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Paired lineage and metal-response analyses identify observational associations

Exact sample-proxy, tissue and age matching produced 203 fibroblast/endothelial pairs across 39 studies. Twenty-four primary cross-lineage associations were estimable; ten passed FDR correction, while fifteen of the separate twenty-four composition-PC sensitivity tests passed their own correction (Fig. 7C). Fibroblast TNFα/NF-κB module expression was positively associated with endothelial Hallmark EMT expression after age/source adjustment (effect 0.425; 95% CI, 0.294–0.557; q = 2.32 × 10⁻⁶). The pair identifier represents a shared sample proxy, not an independently verified donor. Association direction denotes the regression specification and does not establish signaling direction, ligand–receptor binding or endothelial-to-mesenchymal conversion. Residual cell mixture, doublets and shared study context remain alternative explanations.

*Source record: [C0043](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The MT1E/MT1M/MT1X/MT2A expression module was examined in 258 primary state-by-outcome candidates, of which 126 were estimable and twenty-four survived cross-lineage FDR correction (Fig. 7D,E). Positive ROS-module associations included ArtEC_THY1_pan (effect 0.226; 95% CI, 0.200–0.252; q = 1.95 × 10⁻⁵; thirty proxies across eight studies) and Myfibro_CALD1_kidney (0.621; 0.540–0.703; q = 9.73 × 10⁻⁵; 39 proxies across seven studies). Expression outcome modules excluded the MT genes and HMOX1 to reduce direct overlap. Healthy-only analyses constituted a separate, overlapping sensitivity family, with 126 estimable tests and twenty-four significant results; they were not independent replication. Without measured metal exposure or perturbation in these states, the associations cannot distinguish a protective response from damage-associated expression.

*Source record: [C0044](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Public perturbation data provide a bounded test of gene-set response

We additionally reanalyzed the human fibroblast component of GSE325735, associated with an SRSF1 reprogramming study [22], using 136 matched libraries from six donor-derived cell lines (Extended Data Fig. 6A,B). After averaging technical/batch records within donor-condition units, 64 of 192 candidate module-by-condition-by-lineage projections were estimable. None passed BH correction using the primary two-sided donor sign-flip test. Endothelial-selected genes were projected into the same fibroblast experiment and therefore did not constitute an endothelial intervention. Bulk expression lacked the full cell-state-resolved input required for the frozen clock. The reanalysis consequently provides neither significant module-level intervention confirmation under this test nor direct validation of clock rejuvenation; it does not negate findings obtained with different endpoints in the original study.

*Source record: [C0045](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).
