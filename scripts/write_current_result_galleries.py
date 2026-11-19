#!/usr/bin/env python3
"""Write current result navigation and galleries from the reviewed public asset manifest."""
import csv
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def rows(path):
    with (ROOT/path).open() as f:return list(csv.DictReader(f))


def write(path,text):
    p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text.strip()+'\n')


def main():
    release=json.loads((ROOT/'audit/current_release.json').read_text())
    assets=rows('audit/released_assets.csv')
    f3='F3_ds_analysis'
    feature=f3+'/F3_features'
    fn=feature+'/functional_network_20261118'
    lisa=f3+'/TF_analysis/stable_celltype_LISA2_20261118'
    ecm=f3+'/2026_Nat Aging_HAECA/stromal_ECM_20261118'
    cell=feature+'/ModelFeature_cellAge'
    scope='Current snapshot: **19 November 2026**. Counts refer to sample/proxy units unless stated otherwise. Selected expression features are not explicit cell-type proportions; predictive weights and regulatory annotations are not causal evidence.'
    write('README.md',f'''# TabulaStromal

**Human stromal atlases, celltype-resolved chronological-age prediction and downstream molecular associations.**

TabulaStromal integrates fibroblast- and endothelial-lineage references containing **5,916,760** and **2,239,277 cells**. Clock development uses separately filtered sets of **640/437 samples**, grouped into **634/431 GSE/GSM proxies**. These labels are not independently verified donor counts, and the full atlas is not uniformly healthy.

This release updates all five result modules from the source-verified 19 November snapshot, including the downstream presentation, original figures, aggregate tables, evidence records and three English Results sections.

## Model Overview

![TabulaStromal model overview](docs/assets/TabulaStromal_Workflow_V6.png)

*Original workflow schematic. Illustrative numbers, healthy-only/donor labels and GWAS drawings are not additional empirical findings. See the current complete Figure 2 and source-specific definitions in [F2_ModelPerformance](F2_ModelPerformance/).*

## Model Performance

<p align="center">
  <a href="docs/assets/Fig2B_R2_heatmap.png"><img src="docs/assets/Fig2B_R2_heatmap.png" width="32%" alt="Source-fold-4 R squared" /></a>
  <a href="docs/assets/Fig2C_MAE_heatmap.png"><img src="docs/assets/Fig2C_MAE_heatmap.png" width="32%" alt="Source-fold-4 MAE" /></a>
  <a href="docs/assets/Fig2D_MAD_heatmap.png"><img src="docs/assets/Fig2D_MAD_heatmap.png" width="32%" alt="Source-fold-4 age-label MAD" /></a>
</p>

*Original Figure 2 B–D exports. Each model's K is selected using source fold 4. External heatmap values use FT30→ALL100, including adaptation samples. Panel D MAD is mean absolute age-label deviation, not prediction error. These heatmaps are not five-fold means or compact-model validation.*

The current Figure 2 E/F instead shows **development-set in-sample fits** of performance-first ElasticNet models with **750 fibroblast / 350 endothelial Celltype_Gene features**. The current external panels use models selected using **ALL100 outcomes**:

| Current compact-model evaluation | Fibroblast n | Fibroblast R² | MAE, years | Endothelial n | Endothelial R² | MAE, years |
|---|---:|---:|---:|---:|---:|---:|
| Development, in-sample source-model fit | 640 | 0.950 | 4.98 | 437 | 0.950 | 5.03 |
| Adapted ALL100, used in model selection | 93 | 0.377 | 8.71 | 77 | 0.391 | 9.28 |
| Adapted reused 70% subset | 65 | 0.032 | 12.02 | 54 | 0.043 | 11.67 |

**Neither the development fits nor the outcome-selected external estimates establish independent generalization of the current compact models.** The reused 70% outcomes also participated in ALL100 selection. MAD in these current fit/external panels means median absolute prediction error (MedAE), unlike panel D. Exact values and evaluation scopes are linked in [Results Part 2](docs/results/Results_Part2_ClockPerformance_CellAge.md).

The full grouped five-fold benchmark ranked XGBoost first for fibroblasts (mean validation R² 0.803) and ExtraTrees first for endothelial cells (0.716). These are different configurations from the compact ElasticNet fits. The [earlier unreduced EN branch](TabulaStromal_Rebenchmark_2026/) retains its original unadapted and held-out adaptation metrics as historical protocol evidence.

## Explore the five result modules

| Module | Current content |
|---|---|
| [F1_StromalAtlas](F1_StromalAtlas/) | Atlas overview, detailed states, metadata coverage, retrospective annotation agreement and filtering figures |
| [F2_ModelPerformance](F2_ModelPerformance/) | Current complete Figure 2, original overview, B–F panels and exact aggregate metric tables |
| [F3_ds_analysis](F3_ds_analysis/) | Selected features, functional annotation, LISA2, RNA networks, adult raw-count ECM, current CellAge, and the **40-slide PPTX** |
| [TabulaStromal_Rebenchmark_2026](TabulaStromal_Rebenchmark_2026/) | Current 750/350-feature definitions and optimization protocol, with the earlier unreduced EN branch labeled separately |
| [CellAge_Exploration](CellAge_Exploration/) | Historical score provenance and its distinction from the current selected-state CellAge analysis |

Current selected features span **27 states per lineage**, encode **678/337 unique genes**, and share **45 genes**. Functional and regulatory results are exploratory annotations of the selected representation. Relative hUSI state differences are distinct from age slopes; no current primary adjusted hUSI–age or composition–age test survives its applicable correction. Read [Results Part 3](docs/results/Results_Part3_StableFeatures_Mechanisms.md) for complete correction families and coverage limits.

## Read or download the updated Results

| Section | English text | Editable Word | Actual Word-rendered PDF |
|---|---|---|---|
| Stromal atlas | [Part 1](docs/results/Results_Part1_StromalAtlas.md) | [DOCX](docs/results/Results_Part1_StromalAtlas.docx) | [PDF](docs/results/Results_Part1_StromalAtlas.pdf) |
| Clock overview, performance and CellAge context | [Part 2](docs/results/Results_Part2_ClockPerformance_CellAge.md) | [DOCX](docs/results/Results_Part2_ClockPerformance_CellAge.docx) | [PDF](docs/results/Results_Part2_ClockPerformance_CellAge.pdf) |
| Features and molecular associations | [Part 3](docs/results/Results_Part3_StableFeatures_Mechanisms.md) | [DOCX](docs/results/Results_Part3_StableFeatures_Mechanisms.docx) | [PDF](docs/results/Results_Part3_StableFeatures_Mechanisms.pdf) |

[Combined Results Word](docs/results/TabulaStromal_Results_Combined.docx) · [Combined PDF](docs/results/TabulaStromal_Results_Combined.pdf) · [Figure legends](docs/results/Figure_Legends.md) · [Verified references](docs/results/References.md) · [Panel/source manifest](audit/figure_manifest.csv)

## Audit and provenance

The complete five-root metadata discovery contains **{release['source_files']:,} file entries**, including historical and environment files; it does not imply reading all file contents. The public release contains **{release['original_figure_files']} original scientific image/PDF files**, **{release['aggregate_tables']} aggregate/evidence tables** and the downstream PPTX. Core content hashes were rechecked for **{release['source_reads_reverified']} audited source files** before export. See [audit scope and inventories](audit/README.md), [evidence registry](evidence/evidence_registry.csv), [release checks](audit/results_release_validation.json) and [Word rendering receipts](audit/docx_rendering_report.md).

All descendants of the F3 audit directory named `VenEC_CLU_eye_audit` and `F3_features/Age10_20_audit_20261117` are excluded. A legitimate Celltype named `VenEC_CLU_eye` in current aggregate results is not that excluded audit directory. Source images are unaltered. The PPTX retains all 40 slides and embedded image bytes; personal server paths in its XML text are generalized. Full cache-path inventories remain local under ignored `outputs/`.

## Run the source workflows

Use Python **3.11**. No private inputs or pretrained weights are downloaded automatically.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/benchmark.txt
python scripts/run_workflow.py benchmark --help
python scripts/run_workflow.py benchmark --dry-run --lineage both
python scripts/validate_repository.py --verify-checksums
python scripts/validate_results_release.py
python -m unittest discover -s tests -v
```

The source benchmark and inference workflows retain their documented data contracts. The historical final-EN continuation requires its validated checkpoint archive. Current F3 is a curated result release, not a claim that private raw inputs or all execution environments are distributed. See [data requirements](docs/data.md), [reproduction](docs/reproducing.md) and [testing](docs/testing.md). Synthetic checks are software verification, not biological evidence.

## Citation and permissions

[CITATION.cff](CITATION.cff) provides software attribution; no manuscript DOI or published-paper status is asserted. A project-wide license has not yet been selected; see [NOTICE](NOTICE). Third-party resources retain their terms. Contributions should preserve explicit protocol and inference-unit distinctions: [CONTRIBUTING.md](CONTRIBUTING.md).

[中文说明](README.zh-CN.md) · [Publication procedure](docs/publication.md)
''')
    write(f3+'/README.md',f'''# F3 — Selected features and downstream molecular associations

{scope}

[Download the 40-slide downstream presentation](TabulaStromal_F3_Downstream_Results_20261119.pptx) · [Slide index](TabulaStromal_F3_PPT_slide_index.csv) · [Original image provenance](TabulaStromal_F3_PPT_sources.json) · [Results Part 3](../docs/results/Results_Part3_StableFeatures_Mechanisms.md) · [Figure legends](../docs/results/Figure_Legends.md)

The presentation has **40 slides, 54 image placements and 52 unique original images**. Its source images and available PDF partners are reproduced below the corresponding analysis directories. Public PPTX XML paths are generalized; embedded scientific image bytes are unchanged. The PPT is a navigation aid; statistical claims use the underlying tables.

| Analysis | Results and scope |
|---|---|
| [Features and function](F3_features/) | 750/350 selected pairs, 678/337 genes, 27 states each; 45 shared genes; coverage-filtered age displays and complete-library enrichment |
| [LISA2 and RNA networks](TF_analysis/) | Public-chromatin gene-set associations and prior-constrained RNA edges are separate evidence layers; adult TF-age sensitivity has no FDR-significant finding |
| [Adult raw-count ECM](2026_Nat%20Aging_HAECA/) | 39 eligible NB profiles; 67/312 ECM comparisons significant across both lineages; 1/142 displayed mechanosensitive comparisons after full-factor correction |
| [Current ModelFeature CellAge](F3_features/ModelFeature_cellAge/) | Relative hUSI state differences, adjusted age associations, composition and sensitivity analyses using current selected states |

![Fibroblast selected features](F3_features/results/Fibroblast/figures/stable_features_subtype_barplot.png)
![Endothelial selected features](F3_features/results/Endothelial/figures/stable_features_subtype_barplot.png)

These bars count selected model features, not cells or independently tested subtype aging importance. Bootstrap frequency sources use 1,500 resamples, and the representation was selected using ALL100 outcomes. The older 438/23-pair branch is not the current F3 representation.

![Current senescence-associated integration](F3_features/ModelFeature_cellAge/figures/Fig_ModelFeature_CellAge_Integration.png)

Relative hUSI compares a state with other cell-weighted lineage cells in the same sample×tissue unit. It is not age acceleration. No primary adjusted age-dependent hUSI/SSE/composition comparison survives its metric-specific BH family. Shared genes with hUSI limit independence.

The two user-excluded audit trees and all their descendants are absent from this public tree. [Evidence](evidence/evidence_registry.csv), [selected assets](audit/released_assets.csv) and [figure inventory](audit/figure_inventory.csv) identify the actual scope. No expression matrices, fitted models or individual-level score/prediction rows are released.
''')
    write('F2_ModelPerformance/README.md','''# F2 — Clock overview and current Figure 2

![Current complete Figure 2](figures/Fig2_StromalClock_ModelPerformance.png)

[Full-resolution PNG](figures/Fig2_StromalClock_ModelPerformance.png) · [PDF](figures/Fig2_StromalClock_ModelPerformance.pdf) · [Current source manifest](tables/source_manifest.json) · [Results Part 2](../docs/results/Results_Part2_ClockPerformance_CellAge.md)

- **B–D:** source fold 4, K selected separately per model; external FT30→ALL100 includes adaptation samples. Panel D MAD is mean absolute age-label deviation.
- **E/F:** development-set in-sample fits of the current 750/350-feature ElasticNet source models, n=640/437, R²≈0.950, MAE=4.98/5.03 years. MAD means MedAE here.
- **A updated external elements:** n=93/77, R²=0.377/0.391, MAE=8.71/9.28 after adaptation and ALL100-informed selection. These are not independent external-validation estimates.

The full grouped five-fold benchmark and earlier unreduced EN branch are separate protocols. See [the legend](../docs/results/Figure_Legends.md) for remaining conceptual workflow labels and the metric definitions.

![R² heatmap](figures/panels/Fig2B_R2_heatmap.png)
![MAE heatmap](figures/panels/Fig2C_MAE_heatmap.png)
![Age-label MAD heatmap](figures/panels/Fig2D_MAD_heatmap.png)

[Development-fit metrics](tables/Fig2_scatter_metrics.csv) · [Outcome-selected external metrics](tables/Fig2_workflow_external_metrics.csv) · [Source-fold-4 metrics](tables/Fig2_performance_long.csv) · [Evidence](evidence/evidence_registry.csv)
''')
    write('F1_StromalAtlas/README.md','''# F1 — Human stromal atlas

![Stromal atlas](figures/Fig1_StromalAtlas_main.png)

The full reference contains **5,916,760 fibroblast** and **2,239,277 endothelial-lineage cells**, with **4,818/3,553 sample labels** and **848/748 study labels**. These are not verified unique donor counts. The source atlas includes heterogeneous disease annotations.

The UMAPs display 160,000 stratified representative cells per lineage; counts and expression summaries use the full source. Detailed annotation resolves 32/38 states. The `pan`/organ-like Celltype suffixes are naming categories, not original anatomical tissue metadata. Marker/signature agreement is retrospective rather than blinded identity validation.

[Results Part 1](../docs/results/Results_Part1_StromalAtlas.md) · [Full PDF](figures/Fig1_StromalAtlas_main.pdf) · [Figure gallery](figures/) · [Aggregate source tables](tables/) · [Evidence](evidence/evidence_registry.csv) · [File inventory](audit/file_inventory.csv)
''')
    write('TabulaStromal_Rebenchmark_2026/README.md','''# Rebenchmark — Current compact model and earlier EN protocols

The current performance-first representation retains **750/350 Celltype_Gene pairs**, **678/337 unique genes**, and **27 states per lineage**. Feature masks and optimization settings were selected using ALL100 outcomes. Development-set R²≈0.950 is in-sample; adapted ALL100 R²=0.377/0.391 is outcome-selected, not independent validation. The reused 70% subset also participated in selection.

[Current selection summary](tables/current_performance_first_summary.csv) · [Fibroblast features/protocol](tables/Fibro/current_performance_first/) · [Endothelial features/protocol](tables/Endo/current_performance_first/) · [Current downstream interpretation](../F3_ds_analysis/) · [Results Part 2](../docs/results/Results_Part2_ClockPerformance_CellAge.md)

![Current compact EN source composite](figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.png)

*The source composite now depicts the 750/350-feature models. Its nested-GSE display is retrospective and is not a new independent cohort. The reused filename must not be read using the older original-EN panel-source CSV.*

## Earlier unreduced ElasticNet branch

The preserved `final_clock_*` filenames and figures below describe the earlier original-EN branch (K=20,000/10,000), not the compact 750/350-feature model. Its conditional internal R² was 0.742/0.633; unadapted external R² was −2.943/−1.720, with MAE 29.32/24.89 years. Adaptation on 28/23 samples reduced matched held-out MAE to 13.14/13.18 in 65/54 observations, while R² remained −0.001/−0.355. Its 300-bootstrap, frequency-0.70 stable sets contained 438/23 pairs and are not the current F3 feature sets.

[Earlier tables](tables/) · [Earlier figure archive](figures/) · [Evidence](evidence/evidence_registry.csv) · [Registered assets](audit/released_assets.csv)
''')
    write('CellAge_Exploration/README.md','''# CellAge — Score provenance and current-model state analysis

The historical 20261030 snapshot contains hUSI scores for 1,665,482/536,017 fibroblast/endothelial-lineage cells and recSID recommendations for 1,372,858/352,228 cells. Its K8000 partial clock branch is not the current compact clock. hUSI is a senescence-associated transcriptomic score, not age in years; recSID dimensions are identity recommendations, not verified senescence intensities.

The current analysis of the 27 selected states per lineage is published under [F3 ModelFeature_CellAge](../F3_ds_analysis/F3_features/ModelFeature_cellAge/). It covers 731/529 adult sample proxies before state-specific primary coverage filtering. Relative hUSI state contrasts and adjusted age slopes are different tests; no primary adjusted age-dependent hUSI/SSE/composition comparison survives its metric-specific correction.

![Current selected-state CellAge](../F3_ds_analysis/F3_features/ModelFeature_cellAge/figures/Fig_ModelFeature_CellAge_Integration.png)

## Historical snapshot

![Historical CellAge figure](figures/Fig_CellAge_Stromal_Integration.png)

The historical figure is retained as provenance. Its equal-state relative comparator and older clock cannot be substituted for the current cell-weighted selected-state analysis. [Historical aggregate tables](tables/) · [Current Results Part 3](../docs/results/Results_Part3_StableFeatures_Mechanisms.md) · [Evidence](evidence/evidence_registry.csv)
''')
    # Add useful navigation to every released module subdirectory.
    for module in release['modules']:
        paths=[Path(r['repository_file']) for r in assets if r['module']==module]
        directories={ROOT/module/'audit',ROOT/module/'evidence'}
        for p in paths:
            directory=(ROOT/p).parent
            while directory!=ROOT/module and directory.is_relative_to(ROOT/module):directories.add(directory);directory=directory.parent
        for directory in sorted(directories):
            prefix=directory.relative_to(ROOT)
            direct=sorted(p.name for p in directory.iterdir() if p.is_file() and p.name!='README.md')
            child=sorted(p.name for p in directory.iterdir() if p.is_dir())
            content='# '+prefix.name+'\n\n'+scope+'\n\n[Module overview]('+os.path.relpath(ROOT/module/'README.md',directory)+') · [Results and legends]('+os.path.relpath(ROOT/'docs/results',directory)+')\n\n'
            if module=='TabulaStromal_Rebenchmark_2026' and 'current_performance_first' not in str(prefix):
                content+='Files named `final_clock_*` and the 438/23-pair signatures describe the earlier unreduced EN branch; the current compact representation is linked from the module overview.\n\n'
            if module=='CellAge_Exploration':content+='This folder belongs to the historical score snapshot. The current selected-state analysis is linked from the module overview.\n\n'
            images=[p for p in paths if p.suffix=='.png' and (ROOT/p).is_relative_to(directory)]
            for p in images[:8]:content+='!['+p.stem+']('+os.path.relpath(ROOT/p,directory)+')\n\n'
            if child:content+='Subdirectories: '+', '.join('['+n+']('+n.replace(' ','%20')+'/)' for n in child)+'.\n\n'
            if direct:content+='\n'.join('- ['+n+']('+n.replace(' ','%20')+')' for n in direct)+'\n'
            write(prefix/'README.md',content)
    # Override topic landing pages with the relevant statistical interpretation.
    for path,title,note,image in [
        (feature,'Selected features and function','750/350 current pairs with 1,500-resample source frequencies. Age-group displays are descriptive and not differential-expression tests. Main functional networks use Reactome 2022; older signed Reactome 2024 plots are a separate family.',feature+'/results/Fibroblast/figures/top_dynamic_stable_features.png'),
        (f3+'/TF_analysis','Regulatory candidates','LISA2 gene-set association, prior-constrained RNA edges and stricter adult TF-age sensitivity are distinct. No strict adult factor-age association survives lineage-wide FDR; candidates do not establish direct binding or causal regulation.',lisa+'/figures/LISA2_celltype_eligibility.png'),
        (f3+'/2026_Nat Aging_HAECA','Adult count-based ECM analysis','Raw X counts are summed by training proxy×state. Noninteger genes are excluded per profile without rounding. Primary younger/older labels use tissue-specific median ages capped at sixty; ULM uses NB Wald z, not clock coefficients.',ecm+'/Endothelial/figures/ECM_signature_age_contrasts.png'),
        (cell,'Current selected-state CellAge','Within-sample relative hUSI differs from within-state age slopes. BH is metric-specific across both lineages. Heatmaps show row-standardized hUSI and mean sample proportions; gray cells are unavailable. hUSI gene overlap limits independence.',cell+'/figures/Fig_ModelFeature_CellAge_Integration.png')]:
        directory=ROOT/path;children=sorted(p.name for p in directory.iterdir() if p.is_dir())
        write(Path(path)/'README.md','# '+title+'\n\n'+note+'\n\n![Current source figure]('+os.path.relpath(ROOT/image,directory)+')\n\n'+', '.join('['+n+']('+n.replace(' ','%20')+'/)' for n in children)+'\n\n[Full downstream overview]('+os.path.relpath(ROOT/f3/'README.md',directory)+')')
    write('audit/README.md',f'''# Current result audit and public release

The completed five-root metadata discovery contains **{release['source_files']:,} entries**, including dependencies, historical runs and caches. This is metadata discovery, not a claim to read every file. The full path-level inventory remains local in ignored `outputs/source_audit/full_file_inventory.csv.gz`; the public summary accounts for every entry.

Immediately before export, **{release['source_reads_reverified']} audited content files** were rechecked against their stored SHA256 and modification times. The public payload contains {release['original_figure_files']} original image/PDF assets, {release['aggregate_tables']} aggregate/evidence tables and the portable-text version of the 40-slide PPTX. No cell/sample-level observations, predictions, model bundles or matrices are published.

| Record | Meaning |
|---|---|
| [scan_summary.json](scan_summary.json), [source_tree_summary.csv](source_tree_summary.csv) | Five-root complete metadata counts and categories |
| [file_inventory.csv](file_inventory.csv) | Actually read core content, hashes and source-read scope |
| [figure_inventory.csv](figure_inventory.csv) | Result figure candidates after dependency/private-path filtering; published status is separate |
| [released_assets.csv](released_assets.csv) | Reviewed source/export hashes, transformations and protocol scope |
| [figure_manifest.csv](figure_manifest.csv) | Current manuscript Figure 1–6 / Extended Data / supplementary panel assignments |
| [claim_to_evidence_report.csv](claim_to_evidence_report.csv) | {release['claim_mappings']} current Results paragraph/table mappings |
| [Evidence registry](../evidence/evidence_registry.csv) | {release['evidence_records']:,} source-verified aggregate evidence entries; current and historical branches retained separately |
| [numeric_and_statistical_checks.csv](numeric_and_statistical_checks.csv) | 3,749 distinct saved-output checks, zero failures |
| [manuscript_numeric_consistency.csv](manuscript_numeric_consistency.csv) | 121 current main-text rounded-value checks |
| [references_verified.csv](references_verified.csv) | Twenty verified bibliography records |
| [docx_validation_report.md](docx_validation_report.md), [docx_rendering_report.md](docx_rendering_report.md) | Current Results Word content, actual PDF rendering and hash receipts |

Both excluded F3 audit subtrees and all descendants are absent from published paths. Original scientific images are unchanged. Presentation XML paths are generalized, with all embedded image bytes retained. PPT summaries are not replacements for statistical source tables. Counts do not equate sample proxies with independently verified donors.
''')
    write('README.zh-CN.md',f'''# TabulaStromal

本次更新为 **2026-11-19** 经审计结果快照，包含 [F1 图谱](F1_StromalAtlas/)、[F2 模型表现](F2_ModelPerformance/)、[F3 下游机制关联](F3_ds_analysis/)、[Rebenchmark](TabulaStromal_Rebenchmark_2026/) 与 [CellAge](CellAge_Exploration/) 五个模块。

[英文首页及模型图](README.md) · [40 页下游 PPT](F3_ds_analysis/TabulaStromal_F3_Downstream_Results_20261119.pptx) · [三段英文 Results 与 Word/PDF](docs/results/) · [审计清单](audit/README.md)

![当前 Figure 2](F2_ModelPerformance/figures/Fig2_StromalClock_ModelPerformance.png)

当前 compact EN 保留 750/350 个 Celltype_Gene 特征、678/337 个基因、各 27 个状态。开发集 R²≈0.950 是训练内拟合；ALL100 外部 R²=0.377/0.391 的模型选择使用了外部结局，不能视为独立验证；REUSED70 也参与过选择。Figure 2 B–D 是旧 fold-4 benchmark，E/F 是当前开发集拟合。MAD 的两种实际定义已分别说明。

F3 包含稳定入选特征、功能富集、LISA2、RNA prior network、成人 raw-count ECM 和当前选定状态的 CellAge 结果。图和汇总表均已加入网页。共有基因、富集和调控关联不能直接写成共同因果衰老机制；相对 hUSI 状态差异不是随年龄增加或 age acceleration。

明确排除 `VenEC_CLU_eye_audit` 与 `F3_features/Age10_20_audit_20261117` 审计子树及全部后代。合法结果中的 Celltype 名称不等于被排除的审计目录。

完整元数据发现 {release['source_files']:,} 条（含环境/历史文件），实际读取并复核哈希的核心文件 {release['source_reads_reverified']} 个。公开 {release['original_figure_files']} 个原始图像/PDF、{release['aggregate_tables']} 个汇总/证据表和经路径泛化的原始 PPTX。完整 cache 路径清单只留在本地 ignored outputs；不上传表达矩阵、模型或单样本预测/评分行。

运行方法、依赖、可复现性边界和授权信息见[英文首页](README.md)、[数据规范](docs/data.md)、[复现说明](docs/reproducing.md)及 [NOTICE](NOTICE)。本仓库结果草稿不表示论文已经发表。
''')
    write('docs/results.md','''# Current results and interpretation

The 19 November release supersedes the old Results narrative. Read [Part 1](results/Results_Part1_StromalAtlas.md), [Part 2](results/Results_Part2_ClockPerformance_CellAge.md), [Part 3](results/Results_Part3_StableFeatures_Mechanisms.md), [panel legends](results/Figure_Legends.md), and [verified references](results/References.md).

The current compact EN uses 750/350 features; its development fits and ALL100-informed external estimates are not independent validation. The earlier original-EN K20k/K10k branch, 438/23-pair stable sets and old K8000 CellAge join remain historical protocol evidence. They are not substituted for current model metrics or F3 identities. Complete testing families, missing estimates and sample-proxy limitations are retained.

All five modules are linked from the [project homepage](../README.md). The [audit](../audit/README.md) records metadata/content distinctions, excluded trees, source hashes and current paragraph mappings.
''')
    write('results/README.md','''# Aggregate result tables

These existing tables preserve the earlier benchmark, original-EN and CellAge branches for reproducibility. Their `final_clock` names refer to that historical branch. Current compact-model and downstream results are in [F2](../F2_ModelPerformance/), [current Rebenchmark definitions](../TabulaStromal_Rebenchmark_2026/), and [F3](../F3_ds_analysis/). See [the current Results](../docs/results/) for protocol-specific comparisons.
''')
    write('evidence/README.md','''# Source-linked aggregate evidence

The registry preserves verified values and stable evidence IDs from the completed 19 November audit, including explicitly separate current and historical branches. It contains aggregate results, not individual prediction or score rows. Personal server prefixes are generalized. Current paragraph/table mappings are in [the claim report](../audit/claim_to_evidence_report.csv); source/output hashes are in [released assets](../audit/released_assets.csv). The registry can be downloaded when GitHub's CSV preview limit is exceeded.
''')
    write('audit/analysis_status.md','''# Release status — 19 November 2026

Five result roots and the 40-slide F3 presentation were read through the completed manuscript audit; core content SHA256 values were rechecked before public export. Current Results 1–3, source-specific legends, references, figures, aggregate tables and evidence were updated. Excluded F3 trees do not contribute published paths or claims. Word/PDF validity is bound to exact hashes in the rendering report.

Remaining limitations are scientific/reporting boundaries: unresolved donor identity, outcome-informed external selection, unavailable models due to coverage, overlapping signatures and lack of causal/longitudinal outcome validation. These do not imply unperformed experiments were completed.
''')
    write('audit/inconsistency_report.md','''# Source versions and interpretation boundaries

Figure 2 B–D remains source fold 4; E/F is the current compact-model development fit; updated A external elements use ALL100-informed selection. MAD has two actual definitions. Current F3 uses 750/350 pairs and 1,500-resample source frequencies, whereas the earlier 438/23-pair branch used 300 bootstraps. Current main networks use Reactome 2022; older signed Reactome 2024 figures remain distinct. The combined `FINAL_v2` filename was reused for a current compact-model composite; its old panel-source CSV still names the unreduced original EN and is not used to label that image.

Public-chromatin LISA2, RNA-prior edges, adult factor-age sensitivity and raw-count ULM are separate evidence layers. Adult ECM uses tissue-median age cutoffs capped at sixty and unrounded counts with profile-specific noninteger-gene exclusions. Globally selected ECM gene heatmaps show Wald z, not log2FC or exact-feature support in every state. Current CellAge uses cell-weighted other-cell comparisons; hUSI heatmaps are row-standardized, composition is mean sample proportion, and recSID is not senescence intensity.
''')
    write('docs/results/README.md','''# Current English Results and editable documents

[Part 1](Results_Part1_StromalAtlas.md) · [Part 2](Results_Part2_ClockPerformance_CellAge.md) · [Part 3](Results_Part3_StableFeatures_Mechanisms.md) · [Figure legends](Figure_Legends.md) · [Verified references](References.md)

The three Results sections and legends are supplied as separate editable Word files, with a combined Results Word and actual Word-rendered PDFs. Source/evidence links remain clickable. These are research drafts with explicit protocol and statistical limitations, not a claim of paper publication. [Rendering report](../../audit/docx_rendering_report.md).
''')
    print('Updated five module galleries, subdirectory navigation, current scientific scope and public audit documentation')


if __name__=='__main__':main()
