# Current result audit and public release

The completed five-root metadata discovery contains **1,469,578 entries**, including dependencies, historical runs and caches. This is metadata discovery, not a claim to read every file. The full path-level inventory remains local in ignored `outputs/source_audit/full_file_inventory.csv.gz`; the public summary accounts for every entry.

The content-read inventory records **433 files** with source hashes and modification times. The initial export rechecked 432 audited content files; the subsequent [workflow refresh](workflow_refresh.json) reads one additional source image and updates both published copies. The complete metadata scan remains the original snapshot; this targeted image refresh is not a new five-root scan. The public payload contains 229 original image/PDF assets, 201 aggregate/evidence tables and the portable-text version of the 40-slide PPTX. No cell/sample-level observations, predictions, model bundles or matrices are published.

| Record | Meaning |
|---|---|
| [scan_summary.json](scan_summary.json), [source_tree_summary.csv](source_tree_summary.csv) | Five-root complete metadata counts and categories |
| [file_inventory.csv](file_inventory.csv) | Actually read core content, hashes and source-read scope |
| [figure_inventory.csv](figure_inventory.csv) | Result figure candidates after dependency/private-path filtering; published status is separate |
| [released_assets.csv](released_assets.csv) | Reviewed source/export hashes, transformations and protocol scope |
| [workflow_refresh.json](workflow_refresh.json) | Incremental original workflow replacement, dimensions, hashes and unchanged metadata-scan scope |
| [figure_manifest.csv](figure_manifest.csv) | Current manuscript Figure 1–6 / Extended Data / supplementary panel assignments |
| [claim_to_evidence_report.csv](claim_to_evidence_report.csv) | 32 current Results paragraph/table mappings |
| [Evidence registry](../evidence/evidence_registry.csv) | 46,088 source-verified aggregate evidence entries; current and historical branches retained separately |
| [numeric_and_statistical_checks.csv](numeric_and_statistical_checks.csv) | 3,749 distinct saved-output checks, zero failures |
| [manuscript_numeric_consistency.csv](manuscript_numeric_consistency.csv) | 121 current main-text rounded-value checks |
| [references_verified.csv](references_verified.csv) | Twenty verified bibliography records |
| [docx_validation_report.md](docx_validation_report.md), [docx_rendering_report.md](docx_rendering_report.md) | Current Results Word content, actual PDF rendering and hash receipts |

Both excluded F3 audit subtrees and all descendants are absent from published paths. Original scientific images are unchanged. Presentation XML paths are generalized, with all embedded image bytes retained. PPT summaries are not replacements for statistical source tables. Counts do not equate sample proxies with independently verified donors.
