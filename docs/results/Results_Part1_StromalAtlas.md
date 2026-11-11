# Results Part 1 — A human stromal atlas resolves cellular heterogeneity across studies

## Atlas scale and metadata coverage

The integrated collections contained 5,916,760 fibroblast-lineage and 2,239,277 endothelial-lineage cells, associated with 4,818 and 3,553 sample labels and 848 and 748 GSE study labels, respectively ([Fig. 1A,B](../../F1_StromalAtlas/figures/Fig1_StromalAtlas_main.pdf)). These are lineage-specific counts; overlapping contributions across lineages were not treated as additional independent participants. The source objects did not provide a harmonized donor identifier. Sample labels therefore describe the available sampling units, rather than a verified number of unique individuals.

The full-atlas tissue metadata contained 39 observed labels per lineage, including nonstandard categories such as ascites and carcinoma. In contrast, the Tissue_from_Celltype variable used in Fig. 1C–E comprised 12 fibroblast and 11 endothelial annotation suffixes, including pan. These suffixes describe annotation naming and are not interchangeable with anatomical source metadata. The full atlas also contains heterogeneous disease annotations and cannot be described as uniformly healthy. Age-clock development was performed in separately filtered cohorts.

## Hierarchical annotation and full-source expression summaries

Fibroblast annotations comprised Myfibro, eFibro, iFibro and qFibro, resolving into 32 detailed Celltype labels. The endothelial-lineage collection retained ArtEC, CapEC, VenEC, EndoMT, MCs and LEC, resolving into 38 detailed labels ([Extended Data Fig. 1](../../F1_StromalAtlas/figures/ExtendFig1_subtype_umap_composition.pdf)). MCs is a source-defined category within that collection; the retained label does not imply that every included cell is a conventional endothelial cell.

The existing UMAP coordinates and broad-state marker patterns described cellular heterogeneity without generating a new embedding. Marker panels summarize mean stored expression and detection fractions; signature panels show full-source state means, with separate quantitative and annotation color scales (Fig. 1G–J). These summaries support descriptive annotation and share the reference data used for labeling; they are not an independent cell-identity validation. Composition summaries retained unknown ages and described pooled-cell distributions. Cell-type proportions were not used as age-clock features.

## Retrospective signature agreement and exploratory resources

A separate marker-combination analysis evaluated agreement with the existing broad reference annotations ([Extended Data Fig. 2](../../F1_StromalAtlas/figures/ExtFig1_SignatureScore_Fibro_Endo.pdf)). Fibroblast GSE-equal top-signature concordance was 79.27% (95% study-bootstrap CI 75.66–82.56%; 246 retained studies). Endothelial GSE-equal top-signature concordance was 64.27% (95% study-bootstrap CI 58.63–68.85%; 208 retained studies). At the annotation-by-Leiden level, the highest signature agreed with the reference label in 51/55 fibroblast and 32/44 endothelial records. Discordant records were retained in the figure.

These measurements are retrospective annotation agreement: reference labels were marker-informed and the retained studies had been inspected during earlier work. They are not a new blinded accuracy test. The stored filtering-strategy resources and metacell trajectories provide additional descriptive context. In particular, assumed trajectory roots and study-dependent transcriptomic structure do not establish a temporal aging sequence or experimentally observed lineage fate.

The source-to-claim mapping is available in [the evidence registry](../../evidence/evidence_registry.csv) and [the paragraph evidence report](../../audit/claim_to_evidence_report.csv).
