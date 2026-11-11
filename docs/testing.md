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

Checks include: group-separated CV; train-only feature preprocessing and frozen scalers; separate adaptation/test groups; metric formulas and consistent evaluation pairs; corruption detection; all 13 models producing finite predictions and surviving serialization; all seven stages on small synthetic H5AD; resume without rewriting a completed final model; inference feature ordering, missing features, underscore-containing identifiers and unlabeled H5AD inputs.

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
