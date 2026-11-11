# Source-export changes

The original project is not modified. `provenance/source_manifest.csv` records the SHA256 and size of each source file before export. `provenance/export_checksums.csv` records the final exported payload after the following changes.

- Absolute project-home paths are generalized to relative paths. Other researchers' home paths become explicit `external_inputs/<basename>` placeholders; those assets are not supplied.
- Notebook outputs, execution counts and cell metadata are removed. Historical notebook code remains historical; it is not asserted to run without external inputs.
- Text files use LF line endings. Trailing whitespace in Python sources is removed only where the parsed syntax tree remains identical.
- Primary benchmark and CellAge output roots point to ignored `outputs/` directories. Benchmark input/statistical parameters are retained; enrichment paths point to local resource slots and optional historical code-reference paths are removed.
- Atlas default H5AD paths resolve relative to the relocated config, and `Rscript` is taken from the local executable path rather than a personal Conda path.
- Final-clock helpers use installed dependencies rather than a vendored celer directory. Lineage checkpoints can be selected using `TABULASTROMAL_FIBRO_ROOT` and `TABULASTROMAL_ENDO_ROOT`; defaults are `outputs/benchmark/Fibro` and `outputs/benchmark/Endo`. Scientific source audits and signature checks are retained.
- Standalone final-clock inference normalizes equal-valued input buffers to contiguous C order before frozen scaling and prediction. This prevents column reordering or pandas memory layout from selecting different float32 accumulation paths. The remote synthetic test exposed a 1.91×10⁻⁶-year numerical discrepancy before this correction; feature values, scaler parameters, trained weights and released study results are unchanged. Exact-equality tests cover reordered columns, Fortran-order arrays and strided arrays.
- The benchmark digest cache also keys by filesystem device/inode. Atomic replacement with the same byte count and timestamp tick therefore invalidates the cached digest. A regression test simulates coarse timestamp metadata and requires checkpoint corruption detection after replacement.
- Final-clock workers and I/O workers default to 4 and reserve memory to 4 GB; the original bootstrap count, seed, thresholds and adaptation strategy are retained. This changes execution resources, not a claim of an exact bitwise reproduction.
- The old CellAge donor-identity note is replaced with an explicit limitation: the exported repository cannot independently verify that sample labels identify unique donors.
- Numeric result values are unchanged. Path fields are generalized and final-model paths describe the expected local output layout.
- Small helper scripts, documentation, CI and release audits are new. The README displays the supplied project overview and the original Figure 2 B–D panel exports, with their evaluation scope identified in the caption. The separate `clock_performance.png` asset is generated from the released final-clock aggregate metrics.

Checkpoint identities include source/config/input fingerprints. A new export cannot silently reuse an old run under a different identity. Preserve the verified original archive when inspecting historical results, and start a new output directory for new analyses. Archived accelerated runners retain source-layout assumptions and are included for method provenance, not as the default portable command.
