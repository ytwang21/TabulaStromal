# Result audit and release scope

The four original result trees were actually scanned: **1,413,618 file entries**. Discovery read file metadata; it did not imply that every model, cache or expression matrix was loaded. Excluding vendor/runtime/test/backup/source-snapshot copies left **847 figure candidates**. The release selects **117 original PNG/PDF files and 81 aggregate tables**, with explicit versions and evaluation scopes.

| Record | Scope |
| --- | --- |
| [source_tree_summary.csv](source_tree_summary.csv) | Counts and bytes for every scanned file, grouped by module/category/extension |
| [scan_summary.json](scan_summary.json) | Actual scan receipt, roots and local full-inventory location |
| [file_inventory.csv](file_inventory.csv) | Core figure/table files and additional sources actually read for evidence verification |
| [figure_inventory.csv](figure_inventory.csv) | All 847 result-figure candidates; published state, version/protocol, source tables and hashes where read |
| [released_assets.csv](released_assets.csv) | Every selected asset, source SHA256, exported SHA256, transformation and protocol |
| [numeric_and_statistical_checks.csv](numeric_and_statistical_checks.csv) | Fresh numerical checks, including direct metric recomputation and split isolation |
| [claim_to_evidence_report.csv](claim_to_evidence_report.csv) | Results paragraphs mapped to exact evidence IDs and source files |
| [docx_validation_report.md](docx_validation_report.md) | Editable Word text/table/package checks |
| [docx_rendering_report.md](docx_rendering_report.md) | Actual Word-to-PDF rendering, text/link/page checks and hash receipts |
| [analysis_status.md](analysis_status.md) | Completed analyses and limits |
| [inconsistency_report.md](inconsistency_report.md) | Source-version and interpretation differences |

The **complete path-level** inventory is kept locally at `outputs/source_audit/full_file_inventory.csv.gz` and is ignored by Git. It can contain sample-like cache paths. Its public summary accounts for every scanned file. The public file inventory deliberately distinguishes content read, metadata read, and exported assets; it is not falsely labeled as every private cache path.

No expression matrix was densified. H5AD verification read observation metadata in bounded chunks. Existing individual-level predictions were read only to recompute metrics and test split membership; their rows were not released. Figures are unaltered copies of the selected original files. The full-resolution Figure 1 composite is available as PNG/PDF, while its smaller original panels and existing extended-figure previews are used for web navigation.

See the [evidence registry](../evidence/evidence_registry.csv), [three Results sections](../docs/results/), and [four-module gallery navigation](../README.md).
