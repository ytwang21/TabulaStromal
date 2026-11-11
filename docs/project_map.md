# Mapping the original project

This repository is a curated source export, not a recursive copy of the original HPC workspace. Original analysis outputs are preserved outside the repository. The export manifest records original relative paths and SHA256 hashes, without disclosing home directories.

| Original project component | Repository location | Treatment |
| --- | --- | --- |
| `F1_StromalAtlas/code` | `workflows/atlas/code` | Python/R source, configs, output-cleared notebooks |
| `F1_StromalAtlas/extended_figure1_*` | `workflows/atlas/extended_figure1_*` | Signature, subtype and filtering figure source; caches omitted |
| `F2_ModelPerformance` | `workflows/model_performance`, `results/historical_figure2` | Figure source and small aggregate tables |
| `AgePrediction/GeneModel/TabulaStromal_Rebenchmark_2026/tabulastromal` | `workflows/benchmark/tabulastromal` | Main reusable benchmark implementation |
| Rebenchmark final-clock scripts | `workflows/benchmark/final_clock_*`, `run_final_clock.py` | Archive-dependent final EN continuation and standalone inference helpers |
| Fibro accelerated implementation | `workflows/benchmark/Fibro_Accelerated_20261102` | Historical source variant for provenance; no vendor packages or result caches |
| Endothelial acceleration / FT30 variants | `workflows/benchmark` | Historical runtime source; primary supported entry remains the standard benchmark |
| `AgePrediction/CellAge_Exploration` | `workflows/cellage` | Modular pipeline, scripts and tests |
| `wutongdao_TabulaStromal` | `legacy/wutongdao_TabulaStromal` | Annotation and early model notebooks; outputs removed |
| Root `code` | `legacy/code` | Historical integration, annotation inference and exploratory notebooks |
| `data`, `AgePrediction/Valid_data` | Not exported | Expression matrices and individual-level input data |
| `Manuscript` | Not exported | Unpublished drafts, Word documents and large audit files |
| Vendored Python/R packages, model caches, run backups | Not exported | Obtain dependencies separately; keep research archives externally |

Older GeneModel/IDEAL-Age/downstream scripts are superseded for this source release by the explicit rebenchmark and final-clock implementations. Historical notebooks reference some of those external project artifacts; those references do not make them included or portable. Server process watchdogs, scheduler helpers, PID files and automatic restart scripts are omitted because they manage the original running environment.

The primary supported commands are documented in README and `docs/reproducing.md`. Presence of archived code does not imply that its private inputs, optional trajectory tools, old model objects or historical figures were reproduced by the release tests.
