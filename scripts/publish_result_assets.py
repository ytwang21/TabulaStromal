#!/usr/bin/env python3
"""Export verified result figures and aggregate tables from the original project.

No model bundles, cell/sample-level records, or expression matrices are exported.
Figures are copied byte-for-byte; derived tables retain original numeric strings.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil

from audit_source_results import SOURCE_ROOTS, write_csv

ROOT = Path(__file__).resolve().parents[1]
BENCH = Path(SOURCE_ROOTS['TabulaStromal_Rebenchmark_2026'])
CELL = Path(SOURCE_ROOTS['CellAge_Exploration']) / 'runs/20261030_220917_407544'
ATLAS = Path('F1_StromalAtlas/results/runs/meta_full_20261031_232542_2002510')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def clean(text, project):
    text = text.replace(str(project) + '/', '')
    text = re.sub(r'(?<![A-Za-z0-9_])StromalAtlas/', 'F1_StromalAtlas/', text)
    return re.sub(r'/(?:fs|data)/home/[^\s,\x22\x27;<>]+',
                  lambda m: 'external_inputs/' + m.group(0).rstrip('/').split('/')[-1], text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', required=True, type=Path)
    a = parser.parse_args()
    project = a.project_root.resolve()
    assets = []
    sources = []
    figmap = {}
    reads = {}

    def export(source, dest, module, figure_id='', panels='', scope='', table_sources='', predicate=None):
        source, dest = Path(source), Path(dest)
        src, dst = project / source, ROOT / dest
        dst.parent.mkdir(parents=True, exist_ok=True)
        transformation = 'unchanged original bytes'
        if src.suffix.lower() in {'.png', '.pdf', '.jpg', '.svg'}:
            shutil.copyfile(src, dst)
            assert digest(src) == digest(dst)
        elif src.suffix == '.csv':
            with src.open() as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                rows = [row for row in reader if predicate is None or predicate(row)]
            text_fields = [{k: clean(v, project) for k, v in row.items()} for row in rows]
            write_csv(dst, text_fields, fields)
            transformation = ('filtered aggregate rows; numeric strings retained; paths generalized'
                              if predicate else 'aggregate table; numeric strings retained; paths generalized')
        else:
            dst.write_text(clean(src.read_text(), project))
            transformation = 'source text with project paths generalized'
        record = dict(module=module, repository_file=str(dest), source_file=str(source),
                      bytes=dst.stat().st_size, sha256=digest(dst), source_sha256=digest(src),
                      kind='figure' if src.suffix.lower() in {'.png', '.pdf', '.jpg', '.svg'} else 'aggregate_table',
                      figure_id=figure_id, panels=panels, evaluation_scope=scope,
                      source_tables=table_sources, transformation=transformation)
        assets.append(record)
        sources.append(dict(source_path=str(source), repository_path=str(dest),
                            source_bytes=src.stat().st_size, source_sha256=digest(src),
                            transformation=transformation))
        reads[str(source)] = dict(module=module, source_file=str(source), type=src.suffix,
                                 bytes=src.stat().st_size,
                                 mtime_utc=datetime.fromtimestamp(src.stat().st_mtime, timezone.utc).isoformat(),
                                 expected_use=record['kind'], read_status='content_read_and_export_verified',
                                 repository_file=str(dest), source_sha256=digest(src))
        if record['kind'] == 'figure':
            figmap[str(source)] = record
        return record

    # Figure 1 and selected extended/supplementary figure versions.
    module = 'F1_StromalAtlas'
    folder = Path(module) / 'results/fig1_20261103_v2'
    for p in sorted((project / folder).glob('*.png')):
        export(p.relative_to(project), Path(module) / 'figures' / p.name, module, 'Figure 1',
               'A-J' if '_main' in p.name else p.name.split('_')[0].replace('Fig1', ''),
               'Full-atlas descriptive summaries; 160000 display cells per lineage', str(folder))
    export(folder / 'Fig1_StromalAtlas_main.pdf', Path(module) / 'figures/Fig1_StromalAtlas_main.pdf',
           module, 'Figure 1', 'A-J', 'Full-atlas descriptive summaries', str(folder))
    for source, figure_id in [
        ('extended_figure1_subtype/ExtendFig1_subtype_umap_composition_preview.png', 'Extended Data Figure 1'),
        ('extended_figure1_subtype/ExtendFig1_subtype_umap_composition.pdf', 'Extended Data Figure 1'),
        ('extended_figure1_SS/figures/integrated/ExtFig1_SignatureScore_preview.png', 'Extended Data Figure 2'),
        ('extended_figure1_SS/figures/integrated/ExtFig1_SignatureScore_Fibro_Endo.pdf', 'Extended Data Figure 2'),
        ('extended_figure1_filtering_strategy/results/Extended_Figure1_TISCH.pdf', 'Supplementary Figure 1'),
        ('extended_figure1_filtering_strategy/results/Extended_Figure1_HUSCH.pdf', 'Supplementary Figure 2'),
    ]:
        export(Path(module) / source, Path(module) / 'figures' / Path(source).name,
               module, figure_id, 'source-labelled panels',
               'Descriptive annotation support; retrospective agreement; not independent donor validation', source.rsplit('/', 1)[0])
    for p in sorted((project / folder).glob('*.csv')):
        export(p.relative_to(project), Path(module) / 'tables/figure1' / p.name, module)
    for lin in ['fibro', 'endo']:
        for name in ['counts_general.csv', 'counts_subtype.csv', 'counts_celltype.csv', 'counts_tissue.csv',
                     'counts_age_group.csv', 'marker_summary.csv', 'signature_mean.csv']:
            export(ATLAS / lin / name, Path(module) / 'tables' / lin / name, module)
        for name in [f'marker_combination_train_holdout_{lin}.csv', f'plot_reference_umap_signatures_{lin}.csv']:
            export(Path(module) / 'extended_figure1_SS/data' / name, Path(module) / 'tables/signatures' / name, module)
    export(Path(module) / 'extended_figure1_SS/data/plot_reference_cluster_check.csv',
           Path(module) / 'tables/signatures/plot_reference_cluster_check.csv', module)
    for name in ['dataset_audit.csv', 'full_reference_cell_counts.csv', 'full_saved_filter_counts.csv', 'stored_evaluation_all_lineages.csv']:
        export(Path(module) / 'extended_figure1_filtering_strategy/results' / name,
               Path(module) / 'tables/filtering' / name, module)

    # Exact Figure 2 sources; the original overview remains schematic.
    module = 'F2_ModelPerformance'
    for source in ['overview/TabulaStromal_Workflow_V6.png', 'Fig2_StromalClock_ModelPerformance.png',
                   'Fig2_StromalClock_ModelPerformance.pdf']:
        export(Path(module) / source, Path(module) / 'figures' / Path(source).name, module,
               'Figure 2', 'A-F' if 'Fig2_' in source else 'A',
               'Source fold 4; internally selected K; FT30 ALL100 includes adaptation samples',
               'F2_ModelPerformance/data/Fig2_performance_long.csv')
    for p in sorted((project / module / 'panels').glob('*')):
        if p.suffix in ['.png', '.pdf'] and not p.name.startswith('Fig2A'):
            export(p.relative_to(project), Path(module) / 'figures/panels' / p.name, module,
                   'Figure 2', p.name[4], 'Source fold 4; FT30 ALL100', 'F2_ModelPerformance/data')
    for name in ['Fig2_performance_long.csv', 'Fig2_scatter_metrics.csv', 'Fig2_metric_validation.csv']:
        export(Path(module) / 'data' / name, Path(module) / 'tables' / name, module)

    # Final EN figures are distinct from historical FT30/benchmark branches.
    module = 'TabulaStromal_Rebenchmark_2026'
    for lin, rel in [('Fibro', 'Fibro_Accelerated_20261102/results/Fibro'), ('Endo', 'Endo')]:
        root = BENCH / rel
        figures = project / root / 'figures'
        for p in sorted(figures.glob('final_clock_*.png')):
            for suffix in ['.png', '.pdf']:
                src = p.with_suffix(suffix)
                export(src.relative_to(project), Path(module) / 'figures' / lin / src.name,
                       module, 'Figure 3 supporting figure', p.stem,
                       'Final EN; scope specified by filename; full-development stability exploratory',
                       str(root / ('enrichment' if any(x in p.name for x in ['enrichment', 'pathway', 'aging_sets', 'associations']) else 'stability')))
        for name in [f'Fig3_TabulaStromal_{lin}_FINAL.png', f'Fig3_TabulaStromal_{lin}_FINAL.pdf']:
            export(root / 'figures' / name, Path(module) / 'figures' / lin / name,
                   module, 'Figure 3 lineage composite', 'source-labelled panels', 'Final EN protocols', str(root))
        for category, names in {
            'benchmark': ['cv_summary.csv', 'cv_fold_metrics.csv', 'model_comparison.csv'],
            'external': ['final_clock_external_combined_metrics.csv', 'final_clock_external_protocol_comparison.csv', 'final_clock_per_gse_metrics.csv'],
            'stability': ['stable_Celltype_Gene.csv', 'stable_celltypes_ranked.csv', 'stable_genes_unique.csv',
                          'final_clock_internal_cv_summary.csv', 'final_clock_internal_cv_metrics.csv',
                          'final_clock_paired_cv_comparison.csv', 'stable_feature_counts_summary.csv'],
            'enrichment': ['final_clock_enrichment_summary.csv', 'final_clock_gene_set_provenance.csv',
                           'final_clock_stable_gene_age_association.csv'],
        }.items():
            for name in names:
                export(root / category / name, Path(module) / 'tables' / lin / category / name, module)
        export(root / 'enrichment/final_clock_pathway_enrichment.csv',
               Path(module) / 'tables' / lin / 'enrichment/global_pathway_enrichment.csv', module,
               predicate=lambda r: r['scope'] == 'global')
    for name in ['Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2', 'Stromal_FinalClock_shared_specific_genes']:
        for suffix in ['.png', '.pdf']:
            export(BENCH / 'Stromal_Combined_Figure/outputs' / (name + suffix),
                   Path(module) / 'figures/combined' / (name + suffix), module,
                   'Figure 3' if 'FINAL_v2' in name else 'Figure 3 shared-gene supplement',
                   'A-R' if 'FINAL_v2' in name else 'shared and lineage-restricted memberships',
                   'Final EN protocol; paired lineages', str(BENCH / 'final_clock_results_summary.csv'))

    # Latest CellAge result run, using its larger-label display revision.
    module = 'CellAge_Exploration'
    source_figures = CELL / 'figure_revisions/large_text/figures'
    for p in sorted((project / source_figures).rglob('*')):
        if p.is_file() and p.suffix in ['.png', '.pdf']:
            export(p.relative_to(project), Path(module) / 'figures' / p.relative_to(project / source_figures),
                   module, 'Figure 4', 'A-L' if 'Integration' in p.name else p.name.split('_')[0],
                   'Run 20261030_220917_407544; large_text display revision; historical K8000 EN', str(CELL / 'tables'))
    for p in sorted((project / CELL / 'tables').glob('*.csv')):
        if p.name == 'cohort_overview.csv' or any(p.name.endswith(s) for s in [
            '_age_associations.csv', '_clock_associations.csv', '_subtype_contrasts.csv',
            '_subtype_ranking.csv', '_score_feature_gene_overlap.csv', '_score_qc.csv']):
            export(p.relative_to(project), Path(module) / 'tables' / p.name, module)
        elif p.name.endswith('_pathway_enrichment.csv'):
            export(p.relative_to(project), Path(module) / 'tables' / p.name.replace('.csv', '_significant.csv'),
                   module, predicate=lambda r: r.get('q', '') not in ['', 'nan'] and float(r['q']) < .05)

    write_csv(ROOT / 'audit/released_assets.csv', assets, list(assets[0]))
    write_csv(ROOT / 'audit/file_inventory.csv', list(reads.values()), list(next(iter(reads.values()))))
    candidates = []
    with (ROOT / 'outputs/source_audit/figure_candidates.csv').open() as handle:
        for i, row in enumerate(csv.DictReader(handle), 1):
            match = figmap.get(row['source_file'], {})
            candidates.append(dict(module=row['module'], figure_candidate_id=f'FC{i:04d}',
                                   figure_id=match.get('figure_id', 'unselected historical or supporting candidate'),
                                   source_file=row['source_file'], bytes=row['bytes'], mtime_utc=row['mtime_utc'],
                                   panels=match.get('panels', 'not inferred from filename'),
                                   source_tables=match.get('source_tables', ''),
                                   evaluation_scope=match.get('evaluation_scope', 'not promoted without source/protocol review'),
                                   repository_file=match.get('repository_file', ''),
                                   published=bool(match), sha256=match.get('source_sha256', ''),
                                   read_status='original bytes read and copied' if match else 'metadata only'))
    write_csv(ROOT / 'audit/figure_inventory.csv', candidates, list(candidates[0]))
    for module in SOURCE_ROOTS:
        for name, rows in [('file_inventory.csv', [r for r in reads.values() if r['module'] == module]),
                           ('figure_inventory.csv', [r for r in candidates if r['module'] == module]),
                           ('released_assets.csv', [r for r in assets if r['module'] == module])]:
            write_csv(ROOT / module / 'audit' / name, rows, list(rows[0]))
    path = ROOT / 'provenance/source_manifest.csv'
    with path.open() as handle:
        prior = list(csv.DictReader(handle))
    merged = {r['repository_path']: r for r in prior}
    merged.update({r['repository_path']: r for r in sources})
    write_csv(path, list(merged.values()), list(prior[0]))
    print(json.dumps(dict(exported_assets=len(assets), figure_files=sum(r['kind'] == 'figure' for r in assets),
                          total_bytes=sum(r['bytes'] for r in assets), figure_candidates=len(candidates)), indent=2))


if __name__ == '__main__':
    main()
