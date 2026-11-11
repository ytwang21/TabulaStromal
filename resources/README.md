# External analysis resources

`provided_aging_sets.json` contains the original project-supplied lists. Their original label is **provided_unverified**. Do not treat a list name as authoritative source verification.

For actual enrichment, supply the original, appropriately obtained versioned GMT files under `gene_sets/` and record source URL, version, retrieval date and SHA256 in a `download_manifest.json`. The pipeline validates hashes when the manifest provides `Version` and `SHA256` fields. GO/KEGG/Reactome term names in aggregate outputs do not imply that the full gene-set libraries are redistributed here.

An optional `hgnc_complete_set.tsv` should contain approved symbols, HGNC IDs, aliases and previous symbols. Pin its version and hash in your run archive. If absent, the original code reports that syntax-only symbol checks are unverified; do not silently claim approved mapping.

Resource slots are ignored by Git. No download script fabricates the original library version or substitutes a newer library into an archived result.
