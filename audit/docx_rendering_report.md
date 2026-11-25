# Actual Word rendering check

The released DOCX files were opened by LibreOffice and exported to PDF. These are renders of the Word files, not separately typeset Markdown.

Renderer: **LibreOffice 26.2.7.3**. PDF inspection: PyMuPDF 1.26.7.

| Document | Pages | Text blocks | Links | Page bounds | Visual review |
| --- | ---: | ---: | ---: | --- | --- |
| [Figure_Legends.docx](../docs/results/Figure_Legends.docx) | 6 | 34/34 | 0/0 | passed | all pages inspected; no clipping or table overflow observed |
| [Results_Part1_StromalAtlas.docx](../docs/results/Results_Part1_StromalAtlas.docx) | 3 | 30/30 | 4/4 | passed | all pages inspected; no clipping or table overflow observed |
| [Results_Part2_ClockPerformance_CellAge.docx](../docs/results/Results_Part2_ClockPerformance_CellAge.docx) | 4 | 83/83 | 4/4 | passed | all pages inspected; no clipping or table overflow observed |
| [Results_Part3_StableFeatures_Mechanisms.docx](../docs/results/Results_Part3_StableFeatures_Mechanisms.docx) | 6 | 93/93 | 4/4 | passed | all pages inspected; no clipping or table overflow observed |
| [TabulaStromal_Results_Combined.docx](../docs/results/TabulaStromal_Results_Combined.docx) | 13 | 206/206 | 4/4 | passed | all pages inspected; no clipping or table overflow observed |

Every nonempty source paragraph and table cell was found in extracted body text after Unicode/whitespace normalization, using native and geometric PDF reading order. Every external hyperlink target survived PDF export. Pages were checked for blank bodies, expected dimensions and text outside the page.

Times New Roman was requested by the DOCX. This server used the embedded PDF fonts listed in the JSON receipt; the available Roman substitute was Liberation Serif. Pagination can differ on a machine with other fonts or a different Word processor.

The JSON receipt binds each check to the exact DOCX/PDF SHA256. Page PNGs are kept locally under ignored outputs/document_qa/. A passed geometry/text check alone is not visual approval.
