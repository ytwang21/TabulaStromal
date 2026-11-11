# Data contract and availability

Research inputs are not distributed in this Git repository. No public archive accession or model-download URL has been supplied for this release, so none is invented. Original study-level GSE identifiers in the local project do not establish redistribution permission for a combined atlas. Obtain authorized data from the project maintainers or prepare equivalent inputs.

## Clock inputs

The default `workflows/benchmark/config.yaml` retains the original relative project layout under untracked `data/` and `AgePrediction/Valid_data/` paths. Edit a copy under ignored `configs.local/`; `python scripts/configure_local.py --project-root /path/to/LLM_Aging` can point configs to a local copy of the original project. The helper writes new analysis outputs under this repository's ignored `outputs/` tree, not into the original result directories.

| Input | Required content |
| --- | --- |
| Development H5AD | Cells × genes; `obs` fields `Age` (numeric years), `sample`, `GSE`, `GSM`, `Celltype`; gene symbols in `var_names` |
| Combined external H5AD | Compatible gene identifiers, subtype annotations and expression scale; chronological ages required for evaluation |
| Per-GSE external directory | Optional dataset files matching the discovery rules in `tabulastromal/data.py`; BBKNN and alternate Scanorama sources are distinguished |
| Donor crosswalk, if available | Explicit independently verified donor mapping; sample and GSM/GSE proxies do not prove unique donor identity |

The source pipeline resolves conflicting labels, canonicalizes sample keys, rejects incompatible inputs and records schema audits. Review these audit outputs before interpreting a run. Development pseudobulk averages expression within each sample × subtype; the explicit Celltype × Gene schema is retained rather than parsing names at underscores. Missing subtype/gene blocks are zero-filled under the documented model convention. This is a modeling choice, not evidence of zero biological expression.

`X`, `.raw.X`, normalized expression and counts are not interchangeable. The historical benchmark uses configured sources: BBKNN external files use `raw`, Scanorama uses `X`. The final inference CLI requires `--expression raw` where appropriate. Reusing weights requires the same gene identifiers, subtype naming, preprocessing scale and feature schema as training.

## Atlas inputs

`workflows/atlas/code/config.json` resolves H5AD paths relative to the config file. Full atlases use pre-existing `X_umap`, annotations such as `Celltype_general`, `Celltype_sub`, `Celltype`, and available tissue/study/sample metadata. `fields` can override discovered columns. Trajectory analyses use **separate existing metacell inputs** and require the additional Python/R dependencies in `requirements/trajectory.txt` and `docs/reproducing.md`.

The full atlas and age-filtered clock cohorts are different inputs. `Tissue_from_Celltype`-style suffix groups (including `pan`) must not be counted as distinct anatomical organs. Historical tissue labels require harmonization before anatomy-specific claims.

## CellAge inputs

`workflows/cellage/config.yaml` expects full/age/raw atlas files and **precomputed** hUSI and SenCID output tables. It does not train these external scoring methods. Match cell IDs exactly and inspect the audited join rate (`min_match_rate=0.95`).

| Score source | Configured columns | Meaning used by this release |
| --- | --- | --- |
| hUSI | `hUSI`, `husi_SSE`, `husi_GMM` | Senescence-associated outputs; not chronological-age estimates |
| SenCID recommendations | `SenCID_rec_SID1` … `SenCID_rec_SID6` | Recommended SID assignments; not validated senescence intensities |
| Historical clock predictions | model/K/fold/split/sample/age/prediction fields | K=8,000 ElasticNet snapshot with only available held-out folds; distinct from final EN |

Original CellAge files and column names use “donor” for the configured sample aggregation unit. This release does not turn that historical label into evidence of independently verified donor identity. Supply a verified crosswalk before making donor-independence claims.

## External resources and weights

Place appropriately obtained versioned GMT libraries under ignored `resources/gene_sets/`, together with a provenance manifest. Supply HGNC mapping as `resources/hgnc_complete_set.tsv` if needed. The small project-supplied lists in `resources/provided_aging_sets.json` remain explicitly **unverified** and are not redistributed as authoritative CellAge, SASP, or aging databases.

Model bundles are excluded because they contain fitted Python objects and may embed training identifiers. Only load joblib/pickle bundles from a trusted source. The inference command is tested with a synthetic final-EN bundle; no pretrained download is claimed.

## What is public in this release

Public payload: source code, stripped notebooks, configuration templates, original result figures, aggregate performance, bootstrap feature summaries, aggregate CellAge associations, three Results drafts with editable Word versions, and provenance hashes. Excluded: expression matrices, cell/sample mapping files, sample-level predictions and clinical metadata, fitted model bundles, runtime/vendor directories, unreleased manuscript sections, and secrets. These exclusions are enforced by the release audit and `.gitignore`.
