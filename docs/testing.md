# Verification

These checks verify source packaging and computational behavior on synthetic inputs. They do not independently reproduce the entire study or authenticate private data provenance.

## Release audit

```bash
python scripts/validate_repository.py
python scripts/validate_repository.py --verify-checksums
python -m unittest discover -s tests -v
```

The audit checks Python syntax, notebook JSON and cleared outputs, YAML/JSON validity, relative Markdown links, obvious credential patterns, excluded research-file formats, file-size bounds, source-manifest references and consistency of aggregate counts and metrics. A checksum manifest detects post-export changes. An audit result is not a formal security certification.

The synthetic atlas release test also runs the metadata audit from another working directory, resolves input paths against the config location, verifies exact cell/sample counts, checks that donor identity is not invented, confirms unchanged input hashes and reuses completed audit output on resume. It does not execute optional trajectories.

## Benchmark and inference

```bash
python -m pip install -r requirements/benchmark.txt
python -m unittest discover -s workflows/benchmark/tests -p 'test_*.py' -v
```

Checks include: group-separated CV; train-only feature preprocessing and frozen scalers; separate adaptation/test groups; metric formulas and consistent evaluation pairs; corruption detection, including atomic replacement under coarse timestamps; all 13 models producing finite predictions and surviving serialization; all seven stages on small synthetic H5AD; resume without rewriting a completed final model; inference feature ordering, C/Fortran/strided array-layout invariance, missing features, underscore-containing identifiers and unlabeled H5AD inputs.

## CellAge

Run this suite in its module directory so that its historical `scripts` package resolves correctly:

```bash
python -m pip install -r requirements/cellage.txt
(cd workflows/cellage && python -m pytest -q)
python scripts/run_workflow.py cellage --phase smoke-test --output-root outputs/cellage-smoke
```

The CellAge tests exercise ID mapping, score semantics, sample aggregation, gene-feature parsing, atomic writes, checkpoint corruption, interruption and resume. The synthetic pipeline constructs explicit fake samples/scores and checks expected associations; it does not validate real senescence labels.

## Scope not tested in CI

Full private atlas processing, original model-weight equality, actual data access, R/STREAM/Monocle/CellRank trajectories, optional historical accelerators, private archive-dependent final-clock audits, and external hUSI/SenCID training are outside generic CI. Local execution evidence is recorded in `provenance/validation.json`; CI success is established only by the GitHub run itself.

## Results Word files and PDF reading copies

The Markdown, editable Word text and Word table cells are compared directly. Figure and evidence links remain clickable. The included PDF files were rendered from the actual DOCX with LibreOffice; [the page report](../audit/docx_rendering_report.md) records the renderer, fonts, text/link checks and page review. The release validator checks the recorded DOCX/PDF hashes in CI; CI does not launch LibreOffice or repeat visual review.

To rebuild after editing a Results section:

```bash
python -m pip install -r requirements/manuscript.txt
python scripts/build_results_documents.py
python scripts/render_results_documents.py --soffice /path/to/soffice
```

Install LibreOffice separately; the scripts do not download it. Rendered page PNGs stay in ignored `outputs/document_qa/`. Inspect every page, including continued tables and superscript exponents, before recording visual approval in `audit/docx_rendering_report.json` and its Markdown report. Rebuilding either Word or PDF invalidates the previous hash receipt. Fonts can change pagination across machines.

The four gallery pages can be rebuilt using only public release metadata:

```bash
python scripts/write_result_galleries.py
python scripts/validate_results_release.py
```


## Updated fourteen-family results and optional native H2O comparator

Current saved real-data benchmark tables contain fourteen families, including `JungClock`, while the saved Figure 2 heatmaps and self-contained generic CI cover thirteen. The public default retains those thirteen families. To opt in to `JungClock`, append it to `models.names`, install the official H2O Python package 3.44.0.3 and compatible Java, and supply the reviewed level-3 process ontology and author reference resources via `jungclock.ontology`/`reference_dir`. Those external resources are not silently downloaded or relicensed. Configure `python_vendor`, JVM heap and `ice_root` for your installation. The source implementation uses native H2O GLMs and grouped process selection; normalized cell-state inputs omit bulk TPM/fsva, so this is an adaptation rather than exact reproduction of the original bulk model.

Source: https://github.com/saschajung/Clock and https://doi.org/10.1038/s43587-026-01234-8 . Generic synthetic CI does not claim to rerun the native-H2O real-data study. The full fourteen-family saved metrics are independently recomputed in the manuscript evidence audit.
