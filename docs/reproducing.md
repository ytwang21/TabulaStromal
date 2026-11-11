# Reproducing analyses

Run from the repository root with Python 3.11. Use separate virtual environments when historical module dependencies conflict. `requirements/*-historical-lock.txt` record the original runtime; they are separate provenance files, not a claim that every historic package combination is supported on every platform. The active requirements use bounded compatible ranges and are exercised by synthetic tests.

## Benchmark from user-supplied H5AD

```bash
python -m pip install -r requirements/benchmark.txt
# On the original server, or a local copy with the same input tree:
python scripts/configure_local.py --project-root /path/to/LLM_Aging
python scripts/run_workflow.py benchmark --config configs.local/benchmark.yaml --dry-run
python scripts/run_workflow.py benchmark --config configs.local/benchmark.yaml --lineage Fibro --stage audit
python scripts/run_workflow.py benchmark --config configs.local/benchmark.yaml --lineage Fibro --stage all --resume
```

For a different input layout, edit the generated config. Enrichment additionally needs the actual GMT/HGNC resources described in `resources/README.md`. The synthetic integration test supplies its own clearly labelled toy GMT and does not validate the biological libraries.

Stages are `audit`, `features`, `benchmark`, `external`, `stability`, `enrichment`, and `figures`. The workflow enforces prerequisites and records completed artifacts and signatures. `--resume` reuses matching, complete, integrity-checked units. Changed configs, source files or inputs may require a new output root. Do not manually relabel stale checkpoints as valid.

The configured benchmark has 13 models × 4 K values × 5 outer folds per lineage, plus nested linear-model tuning, external analyses and stability fits. Resource requirements depend on cells, samples and feature matrices; no fixed memory/runtime guarantee is made. The full project inputs are much larger than the synthetic test.

## Atlas

```bash
python -m pip install -r requirements/atlas.txt
python scripts/run_workflow.py atlas --config configs.local/atlas.json --run-dir outputs/atlas/run01 --stage audit
python scripts/run_workflow.py atlas --config configs.local/atlas.json --run-dir outputs/atlas/run01 --stage atlas --resume
python scripts/run_workflow.py atlas --config configs.local/atlas.json --run-dir outputs/atlas/run01 --stage markers --resume
```

The `trajectory` stage additionally uses R `Rscript`, Monocle/DDRTree, STREAM/ElPiGraph components and CellRank dependencies. Inspect `workflows/atlas/code/setup_stream_R.R`, `monocle_runner.R`, `metacell_trajectory.py` and the trajectory requirements before installing. These optional dependencies are not installed by the benchmark requirements or tested in generic CI. Trajectories describe transcriptomic states, not longitudinal aging or lineage tracing.

`08_make_fig1.py` renders the later Figure 1 layout from a completed atlas source run. Its JSON `source_run` and output path must be set to the chosen local run. Supplementary signature workflows have their own config under `extended_figure1_SS`; historical configs with missing source files must be reviewed before use.

## Final ElasticNet: archived continuation

This branch preserves the final study analysis, including conservative audit checks. It is not interchangeable with `models/final_model/`, which belongs to the benchmark-wide winning estimator. The final code loads the validated original benchmark/bootstrap archive; a new benchmark does not automatically reproduce every archived figure dependency or checkpoint layout.

Required archive components include `reports/run_identity.json`, feature matrices and schemas, model and preprocessing bundles with matching completion markers, `stability/selection_decision.json`, all bootstrap records, enrichment resources, and the historical `Stromal_Combined_Figure/outputs/source_manifest_5.json` with its image assets. `final_clock_audit.py` and `final_clock_stability.py` enumerate and validate the exact dependencies. They intentionally fail if the archive is incomplete. Checkpoint signatures and paths are part of the archive and must remain consistent.

```bash
python -m pip install -r requirements/final-clock.txt
export TABULASTROMAL_FIBRO_ROOT=/path/to/verified/Fibro
export TABULASTROMAL_ENDO_ROOT=/path/to/verified/Endo
python scripts/run_workflow.py final-clock --lineage Fibro --stage audit --resume
```

Use a **writable working copy** of the original result archive; this continuation writes `final_clock_*` outputs there. Historical figure assets must be restored under `workflows/benchmark/Stromal_Combined_Figure/outputs/` or the corresponding source references adapted explicitly. These archive-specific dependencies are disclosed rather than replaced with synthetic research artifacts. Worker defaults are reduced to 4 for the release; the bootstrap count (300), stability threshold (0.70), seed and scientific rules are retained.

For inference alone, an authentic compatible bundle is sufficient:

```bash
python scripts/predict_age.py --model /path/to/final_clock.joblib \
  --input /path/to/new_cells.h5ad --expression X \
  --output outputs/inference/predictions.csv
```

Named pseudobulk CSV inputs use the sample index in the first column and raw, unscaled Celltype_Gene columns. The saved schema orders features, zero-fills missing features and applies the frozen training scaler. Extra ratio columns do not enter prediction. A raw prediction residual is not a clinically validated measure of biological aging.

## CellAge

```bash
python -m pip install -r requirements/cellage.txt
python scripts/run_workflow.py cellage --config configs.local/cellage.yaml --dry-run
python scripts/run_workflow.py cellage --config configs.local/cellage.yaml --phase all --lineage both
```

CellAge uses its own run IDs, join audits and resumable stages. Its historical `genemodel` inputs deliberately identify K=8,000 results. Changing them to final-clock predictions defines a new analysis; record the change and do not present the archived aggregate tables as the new result. Synthetic verification is available with `--phase smoke-test --output-root outputs/cellage-smoke`.

## Historical figures and notebooks

`workflows/model_performance/make_fig2.py` needs its archived fold-4 inputs and overview PDF. This release includes aggregate Figure 2 metrics, not sample-level plot data or an asserted final-model figure. The source is retained for inspection; its asset paths require the original archive. `scripts/plot_summary.py` independently regenerates the README figure from released aggregate values.

The `legacy/` notebooks retain annotation/integration logic with outputs cleared and personal paths generalized. They require external inputs and their historical environments. They are not covered by a claim of end-to-end execution in this source release.
