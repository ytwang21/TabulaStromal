# Results Part 1. A human stromal atlas resolves heterogeneous cellular states across studies

## Atlas scale and metadata coverage

The reference contained 5,916,760 fibroblast-lineage cells associated with 4,818 sample labels and 848 study labels, and 2,239,277 endothelial-lineage cells associated with 3,553 sample labels and 748 study labels (Fig. 1A,B). Samples and studies could contribute to both lineage collections. These counts represent recorded labels rather than independently verified donors, because a harmonized cross-study donor identifier was unavailable. The atlas constituted the broad expression reference; clock training used a separately filtered subset.

*Source record: [C0001](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Original tissue metadata contained 39 observed labels per lineage, including categories such as ascites and carcinoma that do not constitute a clean anatomical classification. The tissue-like categories displayed in Fig. 1C–E instead derive from the terminal suffix of each Celltype label: twelve fibroblast and eleven endothelial suffixes, including pan. These suffixes describe annotation names and must not be interpreted as original tissue assignments. The full atlas also included disease-associated and incompletely specified samples, precluding a uniformly healthy-cohort description.

*Source record: [C0002](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Hierarchical annotations and descriptive support

Broad fibroblast annotations comprised Myfibro, eFibro, iFibro and qFibro. The endothelial-lineage collection retained ArtEC, CapEC, VenEC, EndoMT, MCs and LEC. These categories resolved into 32 and 38 detailed Celltype labels, respectively (Extended Data Fig. 1A,B). MCs was retained as a source-defined category within the endothelial-lineage collection; its inclusion does not imply that all cells in that collection are conventional endothelial cells. Original UMAP coordinates were preserved, with 160,000 stratified representative cells displayed per lineage. Full-source counts, expression summaries and composition denominators were not reduced to the plotted subset.

*Source record: [C0003](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Marker detection and expression summaries supported broad state-associated patterns (Fig. 1G,H), and signature means provided a complementary description of the annotated states (Fig. 1I,J). These summaries were derived from the same annotated atlas and therefore represent descriptive support rather than independent identity validation. PI16, LUM and IL33 were unavailable in the fibroblast feature universe and were omitted from the corresponding marker panels. Composition across stored age groups and annotation suffixes was calculated from pooled cells, retaining unknown ages (Fig. 1E,F; Extended Data Fig. 1C–F). These cell-weighted proportions were descriptive and were not included as clock features.

*Source record: [C0004](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Retrospective annotation agreement

A separate marker-combination analysis measured agreement with existing broad annotations (Extended Data Fig. 2). The figure-specific selected combination achieved GSE-equal top-signature concordance of 79.27% for fibroblasts (95% study-bootstrap confidence interval, 75.66–82.56%; 246 retained studies) and 64.27% for endothelial-lineage cells (58.63–68.85%; 208 retained studies). At the annotation-by-Leiden-cluster level, the highest signature agreed with the reference in 51 of 55 fibroblast records and 32 of 44 endothelial records. Discordant clusters were retained rather than relabeled to force agreement.

*Source record: [C0005](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

These agreement estimates were retrospective. Reference annotations were marker-informed, and the retained studies had been inspected in earlier analyses; their designation as a source holdout did not establish a newly blinded accuracy test. Additional filtering-strategy figures preserved evaluations from 24 TISCH and 32 HUSCH datasets, comparing recorded Marker, DEG(only) and DEG(like) procedures against their corresponding reference labels (Supplementary Figs. 1,2). Together, the atlas and supporting analyses establish a heterogeneous stromal reference with explicit annotation and metadata limitations, rather than a donor-balanced population sample.

*Source record: [C0006](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Marker-name curation preserves annotation structure

A subsequent annotation review revised marker tokens while preserving the 32 fibroblast and 38 endothelial annotation units, their broad assignments and tissue suffixes. Marker-name changes affected 20 fibroblast units and 32 endothelial units, covering 3,020,209 and 2,000,407 cells. Candidate selection and review used Leiden-stratified discovery/validation samples of 27,500/13,750 fibroblast cells and 21,481/10,740 endothelial cells. The revised labels had no remaining between-unit marker-name collision or tissue-suffix mismatch. Positivity in this sampled review supports the choice of naming markers; it does not establish independent biological purity of every annotated cell.

*Source record: [C0007](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

The annotation revision was propagated by exact cell-identity matching to training-reference objects containing 849,095 fibroblast and 314,925 endothelial cells. This operation updated annotation fields while retaining the target objects' existing clustering and expression. The saved clock and downstream feature lists remain indexed by their frozen, earlier Celltype names. We therefore retain those identifiers in model-specific results and use the annotation crosswalk to distinguish them from the updated atlas vocabulary. A changed naming marker is not evidence that a new biological state was discovered or that a frozen model was retrained.

*Source record: [C0008](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## Integration diagnostics reveal residual study structure

Inspection of the full reference graphs showed that stored BBKNN construction used tissue, rather than GSE, as the batch key. Conditional diagnostics therefore evaluated residual study association within tissue-by-broad-state strata (Extended Data Fig. 4A,B). For matched queries with at least five compatible neighbors, same-study enrichment was 23.82-fold above the conditional random-mixing reference in fibroblasts (95% study-block bootstrap CI, 19.19–28.41) and 16.26-fold in endothelial cells (12.60–21.39). PCA reference values were 24.81 and 16.84, respectively. The modest reduction after graph integration did not remove substantial study-associated structure. Such structure may reflect technical effects and unmodeled biology; these diagnostics cannot uniquely separate them.

*Source record: [C0009](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

Broad-label neighborhood retention in the stored graph was 0.959 and 0.953 at the strongest thirty neighbors, decreasing to 0.772 and 0.695 when as many as 390 neighbors were included. Integration assessment therefore depended on graph scale. Standard scIB metrics were additionally computed on matched reconstructed subsets of 32,863 fibroblast and 33,984 endothelial cells (Extended Data Fig. 4C), rather than on the original full graph. No scVI model was trained, because a verified full-reference raw-count input was unavailable. These results support reporting the atlas as an annotated integrated reference with residual study structure, rather than as a demonstrably batch-free population.

*Source record: [C0010](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

## External annotation processing is distinct from clock validation

An updated query-label-transfer workflow retained 79,285 fibroblast cells from fourteen studies and 47,847 endothelial cells from eleven studies after numerical-age filtering. GSE-aware query BBKNN graphs were kept separate from a reference/query graph used for reference-neighbor voting. Cluster-level assignments yielded seventeen reference labels per lineage across 48 fibroblast and 93 endothelial Leiden clusters; 25 and 35 clusters, respectively, were flagged for low support and were nevertheless assigned their best available reference label. Vocabulary checks across sixty external objects confirmed label compatibility, not independent identity accuracy. These updated cell-level objects differ from the saved 93/77-sample clock-evaluation populations below; no new clock performance is inferred from the annotation transfer alone.

*Source record: [C0011](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*

References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).
