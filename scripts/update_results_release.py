#!/usr/bin/env python3
"""Refresh the public result release from a completed, source-verified manuscript audit.

Only reviewed figures, presentations, aggregate tables and evidence are exported.
The complete filesystem inventory remains in ignored outputs/. Scientific images
are copied without alteration; presentation text paths are made portable.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
MODULES = ['F1_StromalAtlas', 'F2_ModelPerformance', 'F3_ds_analysis',
           'TabulaStromal_Rebenchmark_2026', 'CellAge_Exploration']
EXCLUSIONS = ['F3_ds_analysis/**/VenEC_CLU_eye_audit/**',
              'F3_ds_analysis/F3_features/Age10_20_audit_20261117/**']
NONRESULT = {'vendor', 'runtime', 'site-packages', 'node_modules', 'source_snapshots',
             '__pycache__', '.venv', '.ipynb_checkpoints', 'tests', 'Rlib', 'backups'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None: fields = list(rows[0]) if rows else []
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def excluded(path):
    parts = Path(path).parts
    return ('F3_ds_analysis' in parts and
            ('VenEC_CLU_eye_audit' in parts or 'Age10_20_audit_20261117' in parts))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, required=True)
    args = ap.parse_args(); project = args.project_root.resolve()
    aliases = {str(project), str(args.project_root.absolute())}
    manuscript = project/'Manuscript'
    audit = json.loads((manuscript/'docx_validation.json').read_text())
    assert audit['failed_numeric_checks'] == 0 and audit['failed_manuscript_checks'] == 0
    assert audit['references_ok'] and audit['figure_sources_ok']
    roots = {'F1_StromalAtlas':'F1_StromalAtlas', 'F2_ModelPerformance':'F2_ModelPerformance',
             'F3_ds_analysis':'F3_ds_analysis',
             'TabulaStromal_Rebenchmark_2026':'AgePrediction/GeneModel/TabulaStromal_Rebenchmark_2026',
             'CellAge_Exploration':'AgePrediction/CellAge_Exploration'}

    def clean(value):
        value = str(value)
        for prefix in sorted(aliases, key=len, reverse=True):
            value = value.replace(prefix+'/', '')
        return re.sub(r'/(?:fs|data)/home/[^\s,\x22\x27;<>]+',
                      lambda m: 'external_inputs/'+m.group().rstrip('/').split('/')[-1], value)

    def rel(path): return clean(str(path))
    def module_of(path):
        p = rel(path)
        for key, prefix in roots.items():
            if p == prefix or p.startswith(prefix+'/'): return key
        return 'audit'

    # Every source snapshot used by the preceding audit must still be current.
    source_reads = read_csv(manuscript/'content_read_log.csv')
    for row in source_reads:
        p = Path(row['source_file'])
        assert not excluded(p) and p.exists(), p.name
        assert p.stat().st_mtime_ns == int(row['mtime_ns']), 'Source changed since audit: '+p.name
        assert sha(p) == row['sha256'], 'Source hash changed: '+p.name

    assets = {r['repository_file']:r for r in read_csv(ROOT/'audit/released_assets.csv')}
    provenance = {r['repository_path']:r for r in read_csv(ROOT/'provenance/source_manifest.csv')}
    source_to_repo = {r['source_file']:r['repository_file'] for r in assets.values()}
    def export(source, destination, module, scope='', figure_id='', panels=''):
        src = Path(source)
        if not src.is_absolute(): src = project/src
        assert src.is_file() and not excluded(src), src
        dst = ROOT/destination; dst.parent.mkdir(parents=True, exist_ok=True)
        suffix = src.suffix.lower(); transformation = 'original image bytes unchanged'
        if suffix in {'.png','.pdf','.jpg','.jpeg','.svg'}:
            shutil.copyfile(src, dst); kind = 'figure'
        elif suffix == '.pptx':
            with zipfile.ZipFile(src) as old, zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED) as new:
                for info in old.infolist():
                    raw = old.read(info.filename)
                    if info.filename.endswith(('.xml','.rels')): raw = clean(raw.decode()).encode()
                    new.writestr(info, raw)
            with zipfile.ZipFile(src) as old, zipfile.ZipFile(dst) as new:
                assert new.testzip() is None
                assert all(old.read(n)==new.read(n) for n in old.namelist() if n.startswith('ppt/media/'))
            kind = 'presentation'; transformation = 'PPTX XML paths generalized; slide content and embedded image bytes retained'
        else:
            if suffix=='.csv' and (dst.name=='global_pathway_enrichment.csv' or dst.name.endswith('_pathway_enrichment_significant.csv')):
                rows=read_csv(src)
                fields=list(rows[0])
                rows=[r for r in rows if r.get('scope')=='global'] if dst.name=='global_pathway_enrichment.csv' else [r for r in rows if r.get('q','') not in ('','nan') and float(r['q'])<.05]
                write_csv(dst,[{k:clean(v) for k,v in r.items()} for r in rows],fields)
            else:dst.write_text(clean(src.read_text()))
            kind = 'aggregate_table' if suffix in {'.csv','.tsv'} else 'source_metadata'
            transformation = 'aggregate/evidence text with personal server paths generalized; numeric strings retained'
        record = dict(module=module, repository_file=str(destination), source_file=rel(src),
                      bytes=dst.stat().st_size, sha256=sha(dst), source_sha256=sha(src), kind=kind,
                      figure_id=figure_id, panels=panels, evaluation_scope=scope,
                      source_tables='', transformation=transformation)
        assets[str(destination)] = record; source_to_repo[rel(src)] = str(destination)
        provenance[str(destination)] = dict(source_path=rel(src), repository_path=str(destination),
            source_bytes=src.stat().st_size, source_sha256=sha(src), transformation=transformation)

    # Refresh all previously reviewed source assets; older branches remain labeled historical.
    for row in list(assets.values()):
        scope = row['evaluation_scope']
        if row['module'] in ['TabulaStromal_Rebenchmark_2026','CellAge_Exploration']:
            scope = 'Historical unreduced EN / older CellAge branch; '+scope
        if 'Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2' in row['repository_file']:
            scope='Current compact EN source composite; ALL100 exploratory; nested-GSE rows retrospective, not a new independent cohort. Reused filename; old panel_sources.csv is superseded.'
        if row['module']=='F2_ModelPerformance':
            scope = 'Current Figure 2: B-D source fold 4; E/F development in-sample compact EN; A ALL100 selection-affected'
        export(row['source_file'], row['repository_file'],row['module'],scope,row['figure_id'],row['panels'])

    f2=project/'F2_ModelPerformance'
    for p in sorted((f2/'data').glob('*.csv')):
        if 'prediction' not in p.stem.lower():
            export(p,Path('F2_ModelPerformance/tables')/p.name,'F2_ModelPerformance','Current Fig2 exact source metrics; distinct protocols')
    for name in ['source_manifest.json','validation_report.json']:
        export(f2/'data'/name,Path('F2_ModelPerformance/tables')/name,'F2_ModelPerformance')
    for name in ['TabulaStromal_Workflow_V6.png']:
        export(f2/'overview'/name,Path('docs/assets')/name,'F2_ModelPerformance','Original overview schematic')
    for letter,metric in [('B','R2'),('C','MAE'),('D','MAD')]:
        name=f'Fig2{letter}_{metric}_heatmap.png'
        export(f2/'panels'/name,Path('docs/assets')/name,'F2_ModelPerformance','Source-fold-4 heatmap; age-label MAD')

    f3=project/'F3_ds_analysis'; ppt_sources=json.loads((f3/'TabulaStromal_F3_PPT_sources.json').read_text())
    slides = ppt_sources if isinstance(ppt_sources,list) else ppt_sources.get('slides',ppt_sources.get('slide_sources'))
    assert slides and len(slides)==56
    placed = []
    for slide in slides:
        for im in slide.get('images',[]):
            p=f3/im['path']; assert sha(p)==im['sha256']; placed.append(str(p))
            for q in [p,p.with_suffix('.pdf')]:
                if q.exists(): export(q,q.relative_to(project),'F3_ds_analysis','Current F3; see complete correction family and eligibility in Results/legends')
    assert len(placed)==74 and len(set(placed))==72
    for name in ['TabulaStromal_F3_Downstream_Results_20261119.pptx','TabulaStromal_F3_PPT_slide_index.csv',
                 'TabulaStromal_F3_PPT_sources.json','TabulaStromal_F3_PPT_validation.json']:
        export(f3/name,Path('F3_ds_analysis')/name,'F3_ds_analysis','56-slide downstream summary; navigation rather than statistical primary evidence')

    feature=f3/'F3_features'; fn=feature/'functional_network_20261118'
    lisa=f3/'TF_analysis/stable_celltype_LISA2_20261118'; ecm=f3/'2026_Nat Aging_HAECA/stromal_ECM_20261118'; cell=feature/'ModelFeature_cellAge'
    selected=[]
    for lin in ['Fibroblast','Endothelial']:
        for name in ['stable_features_annotated.csv','celltype_feature_counts.csv','age_heatmap_display_summary.json',
                     'age_group_definition.json','KEGG_stromal_pathway_annotations.csv','KEGG_pathway_gene_matrix_displayed.csv']:
            selected.append(feature/f'results/{lin}/tables'/name)
        selected.extend((fn/lin/'enrichment/global').glob('*.csv'))
        for name in ['gene_model_directions.csv','gene_theme_edges_top3_per_gene.csv','global_significant_pathways.csv','gene_universe.csv']:
            selected.append(fn/lin/'tables'/name)
        for name in ['raw_count_audit.json','age_contrast_definition.json','stable_ECM_feature_age_associations.csv','stable_ECM_gene_age_associations.csv','stable_ECM_gene_age_heatmap_matrix.csv']:
            selected.append(ecm/lin/'tables'/name)
        selected.append(cell/lin/'tables/stable_celltype_model_summary.csv')
    for base,names in [(lisa,['all_celltype_TF_rankings.csv','celltype_run_status.csv','lisa_parameters.json','software_and_parameters.json']),
        (lisa/'integration',['TF_ECM_candidate_evidence_by_feature.csv','candidate_TF_target_ECM_intersections.csv','shared_specific_stable_ECM_genes.csv']),
        (f3/'TF_analysis/results',['analysis_summary.csv','run_config.json']),
        (f3/'TF_analysis/interpretation',['adult_TF_age_sensitivity.csv','adult_TF_age_candidates.csv']),
        (ecm,['all_ECM_ULM_age_contrasts.csv','all_TF_ULM_age_contrasts.csv','all_age_model_coverage.csv','software_and_parameters.json']),
        (cell/'tables',['all_stable_celltype_associations.csv','all_sensitivity_results.csv','model_importance_cellage_correlations.csv','stable_feature_hUSI_signature_overlap.csv'])]:
        selected.extend(base/name for name in names)
    for p in selected:
        export(p,p.relative_to(project),'F3_ds_analysis','Current aggregate F3 evidence; no cell/sample-level observations')

    # New reviewed aggregate evidence. Individual observations and private caches stay local.
    for dirname in ['CelltypeAnno','batchRemove']:
        root=project/'F1_StromalAtlas'/dirname
        aggregate=[]
        if dirname=='CelltypeAnno':
            aggregate=[root/'combined_summary.tsv',root/'acceptance_checks.json',root/'training_annotation_sync_20261122.json',root/'fibro/annotation_mapping.tsv',root/'endo/annotation_mapping.tsv',root/'exter_valid/bbknn_transfer_20261123/combined_summary.tsv']
        else:aggregate=sorted((root/'tables').glob('*.csv'))
        for source in aggregate:
            export(source,Path('F1_StromalAtlas/tables')/(dirname+'_'+('_'.join(source.relative_to(root).parts))),'F1_StromalAtlas','Annotation naming update / residual-study diagnostics; see Results Part 1')
        if dirname=='batchRemove':
            for source in sorted((root/'figures').glob('*')):
                if source.suffix in {'.png','.pdf'}:export(source,Path('F1_StromalAtlas/figures')/source.name,'F1_StromalAtlas','Full-reference graph or reconstructed matched-subset diagnosis; no scVI')
    pattern=feature/'results/SF_pattern_analysis'
    for source in sorted(pattern.rglob('*')):
        if not source.is_file():continue
        if source.suffix in {'.png','.pdf'} or (source.suffix=='.csv' and source.parent.name=='tables' and not any(term in source.name for term in ['GSE_by','tissue_by','leave_one_GSE','source_displayed_row'])):
            export(source,source.relative_to(project),'F3_ds_analysis','Pattern assignments and complete-background enrichment; postselection reliability analysis')
    exploration=f3/'exploration'
    for source in sorted(exploration.glob('task*/all_tests.csv')):
        export(source,source.relative_to(project),'F3_ds_analysis','Full task candidate family; main and sensitivity scopes separate; no individual observations')
    for source in sorted(exploration.glob('task*/**/*')):
        if source.is_file() and source.suffix in {'.png','.pdf'} and source.parent.name=='figures':
            export(source,source.relative_to(project),'F3_ds_analysis','Exploratory associations; source-specific correction families and unavailable tests retained')
    for name in ['analysis_config.json','validation_report.json','software_versions.json']:
        export(exploration/name,Path('F3_ds_analysis/exploration')/name,'F3_ds_analysis','Six-task configuration/validation metadata; no raw records')

    bench=project/roots['TabulaStromal_Rebenchmark_2026']
    for short,folder in [('Fibro','Fibro_Accelerated_20261102/results/Fibro'),('Endo','Endo')]:
        for name in ['performance_first_features.csv','performance_first_summary.csv','performance_first_protocol.json','performance_first_celltypes.csv']:
            p=bench/folder/'stable_clock_optimization'/name
            export(p,Path('TabulaStromal_Rebenchmark_2026/tables')/short/'current_performance_first'/name,
                   'TabulaStromal_Rebenchmark_2026','Current compact EN selected using ALL100 outcomes; not independent validation')
    p=bench/'StableClock_Optimization/balanced_ALL100_selected_summary.csv'
    export(p,Path('TabulaStromal_Rebenchmark_2026/tables/current_performance_first_summary.csv'),'TabulaStromal_Rebenchmark_2026','ALL100 outcome-informed selection')

    # Evidence is aggregate only. IDs and exact numeric strings remain stable.
    evidence=read_csv(manuscript/'evidence_registry.csv')
    assert not any(re.search(r'(?<![A-Za-z0-9])GSM\d+|sample_id=|cellid=', ' '.join(r.values())) for r in evidence)
    export(manuscript/'evidence_registry.csv',Path('evidence/evidence_registry.csv'),'audit','Aggregate evidence; current/historical sources distinguished')
    for name in ['numeric_and_statistical_checks.csv','manuscript_numeric_consistency.csv','references_verified.csv']:
        export(manuscript/name,Path('audit')/name,'audit','Verified source snapshot 2026-11-25')
    for module in MODULES:
        subset=[{k:clean(v) for k,v in r.items()} for r in evidence if module_of(r['source_file'])==module]
        destination=Path(module)/'evidence/evidence_registry.csv'
        write_csv(ROOT/destination,subset,list(evidence[0]))
        assets[str(destination)]=dict(module=module,repository_file=str(destination),source_file='Manuscript/evidence_registry.csv',
            bytes=(ROOT/destination).stat().st_size,sha256=sha(ROOT/destination),source_sha256=sha(manuscript/'evidence_registry.csv'),
            kind='aggregate_table',figure_id='',panels='',evaluation_scope='Module-filtered aggregate evidence; current/historical source paths explicit',source_tables='',
            transformation='Filtered by original source module; numeric strings retained; server paths generalized')
        provenance[str(destination)]=dict(source_path='Manuscript/evidence_registry.csv',repository_path=str(destination),
            source_bytes=(manuscript/'evidence_registry.csv').stat().st_size,source_sha256=sha(manuscript/'evidence_registry.csv'),
            transformation=assets[str(destination)]['transformation'])
    read_rows=[dict(module=module_of(r['source_file']),source_file=rel(r['source_file']),type=Path(r['source_file']).suffix,
                   bytes=r['bytes'],mtime_ns=r['mtime_ns'],source_sha256=r['sha256'],read_status='content read; source SHA256 rechecked before public export') for r in source_reads]
    write_csv(ROOT/'audit/file_inventory.csv',read_rows)

    current_claims=read_csv(manuscript/'claim_to_evidence_report.csv'); public_claims=[]
    names={'03_Results_Part1_StromalAtlas':'Results_Part1_StromalAtlas',
           '04_Results_Part2_ClockPerformance_CellAge':'Results_Part2_ClockPerformance_CellAge',
           '05_Results_Part3_StableFeatures_Mechanisms':'Results_Part3_StableFeatures_Mechanisms'}
    for old,new in names.items():
        text=(manuscript/(old+'.md')).read_text(); blocks=[]; last_heading=''
        for block in text.strip().split('\n\n'):
            blocks.append(block)
            if block.startswith('#'): last_heading=block.lstrip('# ')
            matches=[r for r in current_claims if r['claim']==block or (r['section'].startswith('Table') and r['claim']==block)]
            for r in matches:
                public_claims.append(dict(section_file='docs/results/'+new+'.md',heading=last_heading,paragraph=len(blocks),
                    claim_text=block,claim_ids=r['evidence_registry_ids'],source_files=clean(r['source_file']),status='source-verified; protocol and exploratory limits retained',claim_id=r['claim_id']))
                blocks.append('*Source record: ['+r['claim_id']+'](../../audit/claim_to_evidence_report.csv); figure files are indexed in the [panel manifest](../../audit/figure_manifest.csv).*')
        blocks.append('References: [verified bibliography](References.md). Figures: [source-matched legends](Figure_Legends.md).')
        (ROOT/'docs/results'/(new+'.md')).write_text('\n\n'.join(blocks)+'\n')
    assert len(public_claims)==len(current_claims),(len(public_claims),'paragraph/table mappings')
    write_csv(ROOT/'audit/claim_to_evidence_report.csv',public_claims)
    for old,new in [('08_Figure_Legends.md','Figure_Legends.md'),('09_References.md','References.md')]:
        (ROOT/'docs/results'/new).write_text(clean((manuscript/old).read_text()))
    fm=[]
    for row in read_csv(manuscript/'figure_manifest.csv'):
        row={k:clean(v) for k,v in row.items()}
        for field in ['figure_file','preview_file']:
            if row[field] not in source_to_repo:
                module=module_of(row[field]);dst=Path(module)/'figures'/Path(row[field]).name
                export(row[field],dst,module,'Current figure source/preview; see assigned panel legend',row['figure'],row['figure_panels'])
        for field in ['figure_file','preview_file']: row[field]=source_to_repo[row[field]]
        row['table_repository_file']=source_to_repo.get(row['source_file'],'')
        fm.append(row)
    write_csv(ROOT/'audit/figure_manifest.csv',fm)

    # Publish complete counts, keep all path-level cache metadata local.
    local=ROOT/'outputs/source_audit';local.mkdir(parents=True,exist_ok=True)
    totals=defaultdict(lambda:[0,0]);count=0
    fields=['module','source_file','type','bytes','mtime_utc','expected_use','read_status']
    with gzip.open(local/'full_file_inventory.csv.gz','wt',compresslevel=1,newline='') as out:
        writer=csv.DictWriter(out,fieldnames=fields,lineterminator='\n');writer.writeheader()
        with (manuscript/'file_inventory.csv').open() as source:
            for r in csv.DictReader(source):
                assert not excluded(r['path']);module=module_of(r['path']);source_file=rel(r['path']);category=r['expected_use']
                if any(p in NONRESULT for p in Path(source_file).parts): category='runtime_vendor_backup_or_test'
                row=dict(module=module,source_file=source_file,type=r['type'],bytes=r['bytes'],mtime_utc=r['mtime_utc'],expected_use=category,read_status=r['read_status'])
                writer.writerow(row);totals[(module,category,r['type'])][0]+=1;totals[(module,category,r['type'])][1]+=int(r['bytes']);count+=1
    summary=[dict(module=k[0],category=k[1],extension=k[2],n_files=v[0],bytes=v[1]) for k,v in sorted(totals.items())]
    write_csv(ROOT/'audit/source_tree_summary.csv',summary)
    scan=json.loads((manuscript/'inventory_summary.json').read_text());assert count==scan['total']
    scan.update(scanned_files=count,source_roots=roots,complete_inventory='outputs/source_audit/full_file_inventory.csv.gz',complete_inventory_public=False)
    (ROOT/'audit/scan_summary.json').write_text(json.dumps(scan,indent=2)+'\n')
    candidates=[]
    for r in read_csv(manuscript/'figure_inventory.csv'):
        path=rel(r['path'])
        if excluded(path) or any(p in NONRESULT for p in Path(path).parts) or re.search(r'GSM\d+',path): continue
        dst=source_to_repo.get(path,'')
        candidates.append(dict(module=module_of(path),source_file=path,format=r['format'],bytes=r['bytes'],mtime_utc=r['mtime_utc'],
            figure_id=r['main_figure_suitability'],panels=r['panel_semantics'],published=bool(dst),repository_file=dst,
            source_sha256=assets[dst]['source_sha256'] if dst else '',scope='Metadata candidate; source evidence is required for a scientific claim'))
    write_csv(ROOT/'audit/figure_inventory.csv',candidates)
    write_csv(ROOT/'audit/released_assets.csv',list(assets.values()))
    write_csv(ROOT/'provenance/source_manifest.csv',list(provenance.values()))
    for module in MODULES:
        write_csv(ROOT/module/'audit/file_inventory.csv',[r for r in read_rows if r['module']==module],list(read_rows[0]))
        write_csv(ROOT/module/'audit/figure_inventory.csv',[r for r in candidates if r['module']==module],list(candidates[0]))
        write_csv(ROOT/module/'audit/released_assets.csv',[r for r in assets.values() if r['module']==module],list(next(iter(assets.values()))))
    release=dict(date='2026-11-25',modules=MODULES,exclusions=EXCLUSIONS,source_files=count,
        source_reads_reverified=len(source_reads),evidence_records=len(evidence),claim_mappings=len(public_claims),
        figure_candidates=len(candidates),original_figure_files=sum(r['kind']=='figure' for r in assets.values()),
        aggregate_tables=sum(r['kind']=='aggregate_table' for r in assets.values()),presentations=sum(r['kind']=='presentation' for r in assets.values()),
        downstream_slides=len(slides),downstream_image_placements=len(placed),downstream_unique_original_images=len(set(placed)),
        raw_data_or_individual_predictions_published=False,complete_inventory_public=False)
    (ROOT/'audit/current_release.json').write_text(json.dumps(release,indent=2)+'\n')
    print(json.dumps(release,indent=2))


if __name__=='__main__': main()
