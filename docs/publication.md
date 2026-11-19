# Publishing reviewed repository updates

An update is complete when the reviewed local commit is on the remote main branch and the fetched remote head matches it. The current public payload contains original figures, aggregate tables, a portable-text downstream PPTX and three Results drafts; it does not publish private matrices, models or the remaining manuscript sections.

Refresh the source assets using the scripts documented in [scripts/README.md](../scripts/README.md), build/render the Results files, review the pages, and run both repository and result-release checks. Then commit the reviewed payload and push main with normal fast-forward Git semantics. Fetch before pushing; do not force-push over other work. Verify the remote commit and its GitHub Actions checks after upload.

The original `scripts/publish_github.py` only creates an empty repository and intentionally refuses an existing nonempty repository. It is retained for provenance and is not the update procedure. Existing GitHub authorization can be reused; never include authentication secrets in repository files or publication receipts.
