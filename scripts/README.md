# Project scripts

The primary analysis entry points remain `run_workflow.py` and `predict_age.py`; see the project homepage and data contract. The original-EN continuation requires a validated archive.

## Current result publication

Run from the repository root with a completed, source-verified manuscript audit in the local project:

```bash
python scripts/update_results_release.py --project-root /path/to/LLM_Aging
python scripts/write_current_result_galleries.py
python scripts/build_results_documents.py
python scripts/render_results_documents.py --soffice /path/to/soffice
python scripts/validate_results_release.py
python scripts/validate_repository.py --write-checksums
python scripts/validate_repository.py --verify-checksums
```

The updater rechecks all audited source hashes, rejects the two excluded F3 trees, exports reviewed aggregate figures/tables, generalizes presentation XML paths while preserving embedded images, and keeps the complete path inventory in ignored outputs. Individual rows and model/matrix files are not exported. The five result documents require actual PDF rendering and manual page review before the release validator accepts their receipts.

The November 9 audit/export/Results-writing helpers are preserved as historical source and their direct CLI entry points are disabled to prevent overwriting the current five-module release with stale protocols. The current updater and gallery writer replace those publication operations. `publish_github.py` is the original empty-repository creator; updates to this existing repository use a reviewed fast-forward Git push.
