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

References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).
