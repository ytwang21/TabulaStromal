# Release and analysis helpers

[Analysis entry points](../docs/reproducing.md) · [Result galleries](../docs/results/)

The source-result release is rebuilt from authorized original inputs in this order:

```bash
python scripts/audit_source_results.py --project-root /path/to/LLM_Aging
python scripts/publish_result_assets.py --project-root /path/to/LLM_Aging
python scripts/build_result_evidence.py --project-root /path/to/LLM_Aging
python scripts/write_results_sections.py
python scripts/write_result_galleries.py
python scripts/build_results_documents.py
python scripts/validate_repository.py --write-checksums
python scripts/validate_results_release.py
```

The first script stores a complete metadata inventory locally under ignored outputs/. Figure files are copied without alteration. Only aggregate tables are released. The evidence builder reads observation metadata and existing predictions for checks; it does not export individual prediction rows or densify expression. Markdown and Word share a single source. `publish_result_assets.py` exports local files and does not itself publish to GitHub.
