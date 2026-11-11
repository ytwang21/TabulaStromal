#!/usr/bin/env python3
"""Build navigable result-directory READMEs from the released asset manifest."""
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    with (ROOT / 'audit/released_assets.csv').open() as f:
        assets = list(csv.DictReader(f))
    facts = json.loads((ROOT / 'audit/verified_facts.json').read_text())
    source_map = json.loads((ROOT / 'audit/scan_summary.json').read_text())['source_roots']
    result_files = {
        'F1_StromalAtlas': 'Results_Part1_StromalAtlas',
        'F2_ModelPerformance': 'Results_Part2_ClockPerformance_CellAge',
        'TabulaStromal_Rebenchmark_2026': 'Results_Part3_StableFeatures_Mechanisms',
        'CellAge_Exploration': 'Results_Part2_ClockPerformance_CellAge',
    }
    workflows = {'F1_StromalAtlas':'atlas','F2_ModelPerformance':'model_performance',
                 'TabulaStromal_Rebenchmark_2026':'benchmark','CellAge_Exploration':'cellage'}

    def pair(a, b, alt_a, alt_b):
        return f'<p>\n  <a href="{a}"><img src="{a}" width="49%" alt="{alt_a}" /></a>\n  <a href="{b}"><img src="{b}" width="49%" alt="{alt_b}" /></a>\n</p>\n'

    heroes = {
        'F1_StromalAtlas':'F1_StromalAtlas/figures/ExtendFig1_subtype_umap_composition_preview.png',
        'F2_ModelPerformance':'F2_ModelPerformance/figures/TabulaStromal_Workflow_V6.png',
        'TabulaStromal_Rebenchmark_2026':'TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.png',
        'CellAge_Exploration':'CellAge_Exploration/figures/Fig_CellAge_Stromal_Integration.png',
    }
    for module in source_map:
        own = [r for r in assets if r['module'] == module]
        title = {'F1_StromalAtlas':'Stromal Atlas','F2_ModelPerformance':'Clock Overview and Source-Fold Performance',
                 'TabulaStromal_Rebenchmark_2026':'Rebenchmark, Final ElasticNet and Stable Features',
                 'CellAge_Exploration':'CellAge: Senescence-Associated Scores and Clock Integration'}[module]
        text = [f'# {title}', '',
                f'Source project directory: `{source_map[module]}`. This page presents the selected verified result versions, original figure files, and aggregate evidence.', '',
                f'**Results:** [English text](../docs/results/{result_files[module]}.md) · [Editable Word](../docs/results/{result_files[module]}.docx) · [PDF preview](../docs/results/{result_files[module]}.pdf) · [Figure legends](../docs/results/Figure_Legends.md)', '',
                '**Audit:** [File inventory](audit/file_inventory.csv) · [Figure inventory](audit/figure_inventory.csv) · [Evidence registry](evidence/evidence_registry.csv) · [Released asset hashes](audit/released_assets.csv)', '',
                f'**Code:** [Workflow source](../workflows/{workflows[module]}/) · [Reproduction instructions](../docs/reproducing.md)', '']
        if module == 'F1_StromalAtlas':
            text += ['## Atlas and annotation coverage', '',
                     '| Lineage | Cells | Sample labels | GSE labels | Detailed Celltype labels |',
                     '| --- | ---: | ---: | ---: | ---: |']
            for lin,row in facts['atlas'].items():
                text.append(f'| {lin} | {row["n_cells"]:,} | {row["n_sample_labels"]:,} | {row["n_GSE_labels"]} | {row["n_Celltypes"]} |')
            text += ['', 'Sample labels are not a verified count of unique donors. The full atlas includes heterogeneous disease metadata. Figure 1C–E use annotation suffixes including pan; these are different from original tissue labels.', '',
                     '## Figure 1: full atlas and expression support', '',
                     '[Original complete PNG](figures/Fig1_StromalAtlas_main.png) · [Complete vector PDF](figures/Fig1_StromalAtlas_main.pdf)', '',
                     'The full-resolution composite is available above. Original marker, metadata and signature panels are displayed separately for navigation.', '',
                     pair('figures/Fig1G_Fibroblast_markers.png','figures/Fig1H_Endothelial_markers.png','Figure 1G fibroblast markers','Figure 1H endothelial markers'),
                     pair('figures/Fig1E_tissue_summary.png','figures/Fig1F_age_summary.png','Figure 1E annotation suffix counts','Figure 1F age-group counts'),
                     pair('figures/Fig1I_Fibroblast_signatures.png','figures/Fig1J_Endothelial_signatures.png','Figure 1I fibroblast signatures','Figure 1J endothelial signatures'),
                     '## Extended Data Figure 1: detailed subtypes and composition', '',
                     '![Detailed annotation and composition](figures/ExtendFig1_subtype_umap_composition_preview.png)', '',
                     '[Original PDF](figures/ExtendFig1_subtype_umap_composition.pdf). Pooled-cell proportions are descriptive and are not clock inputs.', '',
                     '## Extended Data Figure 2: marker-signature agreement', '',
                     '![Signature agreement preview](figures/ExtFig1_SignatureScore_preview.png)', '',
                     '[Original PDF](figures/ExtFig1_SignatureScore_Fibro_Endo.pdf). This is retrospective reference-label agreement, not a new blinded annotation test.', '',
                     '## Supplementary filtering figures', '',
                     '[Supplementary Figure 1 — TISCH](figures/Extended_Figure1_TISCH.pdf) · [Supplementary Figure 2 — HUSCH](figures/Extended_Figure1_HUSCH.pdf)', '']
        elif module == 'F2_ModelPerformance':
            text += ['## Model Overview', '', '![Supplied model overview](figures/TabulaStromal_Workflow_V6.png)', '',
                     'The supplied overview is schematic. Its illustrative values, healthy-only label and GWAS/TF sketches are not substituted for audited study evidence.', '',
                     '## Model Performance: panels B–D', '', '<p>',
                     '  <img src="figures/panels/Fig2B_R2_heatmap.png" width="32%" alt="Figure 2B R2" />',
                     '  <img src="figures/panels/Fig2C_MAE_heatmap.png" width="32%" alt="Figure 2C MAE" />',
                     '  <img src="figures/panels/Fig2D_MAD_heatmap.png" width="32%" alt="Figure 2D age-label MAD" />', '</p>', '',
                     'These are source-fold-4 values with internally selected K, not fivefold averages. External uses FT30 ALL100, including adaptation samples. MAD is mean absolute deviation of age labels, not median prediction error.', '',
                     '## ElasticNet example: panels E–F', '',
                     pair('figures/panels/Fig2E_Fibro_squareplot.png','figures/panels/Fig2F_Endo_squareplot.png','Source-fold-4 fibroblast EN','Source-fold-4 endothelial EN'),
                     'Both displayed EN examples use K=20,000. They differ from the later full-development final clocks. The historical biological-age axis wording refers to the recorded chronological-age target.', '',
                     '[Full Figure 2 PNG](figures/Fig2_StromalClock_ModelPerformance.png) · [Full Figure 2 PDF](figures/Fig2_StromalClock_ModelPerformance.pdf) · [Final-clock gallery](../TabulaStromal_Rebenchmark_2026/)', '']
        elif module == 'TabulaStromal_Rebenchmark_2026':
            text += ['## Figure 3: final EN protocols and interpretation', '',
                     '![Final EN clocks and stable features](figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.png)', '',
                     '[Original PDF](figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.pdf). Fibroblast/endothelial final EN uses K=20,000/10,000. Stable features remain explanatory because reduced EN models predict less well.', '',
                     '## Independent external evaluation and adaptation', '',
                     pair('figures/Fibro/final_clock_external_combined_noadapt.png','figures/Endo/final_clock_external_combined_noadapt.png','Fibroblast no-adaptation external','Endothelial no-adaptation external'),
                     pair('figures/Fibro/final_clock_external_ft30_comparison.png','figures/Endo/final_clock_external_ft30_comparison.png','Fibroblast distinct external scopes','Endothelial distinct external scopes'),
                     'No-adaptation merged external n=93/77. FT30 uses 28/23 fitted and 65/54 held-out samples. ALL100 contains fitted observations and must not be described as independent validation.', '',
                     '## Stable cell types and functional overlaps', '',
                     pair('figures/Fibro/final_clock_stable_celltypes.png','figures/Endo/final_clock_stable_celltypes.png','Fibroblast candidate-normalized stable celltypes','Endothelial candidate-normalized stable celltypes'),
                     pair('figures/Fibro/final_clock_GO_Reactome_KEGG_enrichment.png','figures/Endo/final_clock_GO_Reactome_KEGG_enrichment.png','Fibroblast enrichment','Endothelial enrichment'),
                     '![Shared stable-gene membership](figures/combined/Stromal_FinalClock_shared_specific_genes.png)', '',
                     'Bootstrap frequency, model coefficients and enrichment are exploratory associations. No causal or perturbation-validated mechanism is claimed. Global enrichment tables retain the original tested terms and FDR family.', '']
        else:
            text += ['## Figure 4: historical CellAge integration', '',
                     '![CellAge integrated result](figures/Fig_CellAge_Stromal_Integration.png)', '',
                     '[Original PDF](figures/Fig_CellAge_Stromal_Integration.pdf). Source run: 20261030_220917_407544; large_text is a display revision of that run, not a new analysis.', '',
                     '## Subtype scores and adjusted age effects', '',
                     pair('figures/panels/C_fibro_ranking.png','figures/panels/D_endo_ranking.png','Fibroblast hUSI subtype means','Endothelial hUSI subtype means'),
                     pair('figures/panels/E_fibro_age.png','figures/panels/F_endo_age.png','Fibroblast adjusted hUSI age effects','Endothelial adjusted hUSI age effects'),
                     'The primary adjusted models contain 15 identifiable fibroblast and 6 endothelial subtype fits; none has BH q<0.05. Source donor terminology refers to configured sample aggregation and is not independent identity verification.', '',
                     '## SID recommendations and historical clock context', '',
                     pair('figures/panels/G_fibro_sid.png','figures/panels/H_endo_sid.png','Fibroblast SID recommendation landscape','Endothelial SID recommendation landscape'),
                     'recSID outputs are identity recommendations, not genuine SenCID senescence intensities. hUSI is not age in years. Panels I–L use older K=8,000 EN/stable sets and do not validate the current final clocks. Significant-only historical pathway tables preserve the source q values and their original correction family.', '']
        text += ['## Released file catalogue', '', '<details>', '<summary>Original figure and aggregate-table files</summary>', '',
                 '| File | Type | Size (MB) |', '| --- | --- | ---: |']
        for row in own:
            rel = Path(row['repository_file']).relative_to(module)
            text.append(f'| [{rel}]({rel.as_posix()}) | {row["kind"]} | {int(row["bytes"])/1e6:.2f} |')
        text += ['', '</details>', '',
                 'All released images/PDFs are original copies. The figure inventory also lists unselected historical candidates as metadata-only; their presence does not establish a verified or published result.', '']
        (ROOT/module/'README.md').write_text('\n'.join(text))

    # Give figure, table, audit and evidence subdirectories their own useful indices.
    for module in source_map:
        directories = {p.parent for p in (ROOT/module).rglob('*') if p.is_file()}
        for folder in sorted(directories):
            if folder == ROOT/module or (folder/'README.md').exists():
                continue
            parent_link = os.path.relpath(ROOT/module/'README.md',folder)
            hero = os.path.relpath(ROOT/heroes[module],folder)
            lines = [f'# {module} — {folder.name}', '', f'[Analysis overview and protocol]({parent_link})', '',
                     f'![Related result figure]({hero})', '', '## Files', '']
            for p in sorted(folder.iterdir()):
                if p.is_file():lines.append(f'- [{p.name}]({p.name})')
            lines += ['', 'Figure and data provenance are recorded in the analysis directory audit and evidence registries.', '']
            (folder/'README.md').write_text('\n'.join(lines))
    print('Wrote four analysis galleries and directory indices.')


if __name__ == '__main__':
    main()
