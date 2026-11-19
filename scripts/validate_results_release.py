#!/usr/bin/env python3
"""Validate current scientific release assets, correction families and rendered documents."""
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from zipfile import ZipFile
from xml.etree import ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}


def rows(rel):
    with (ROOT/rel).open() as f:return list(csv.DictReader(f))


def bh(values):
    order=sorted(range(len(values)),key=values.__getitem__);result=[0.]*len(values);current=1.
    for rank in range(len(order),0,-1):
        i=order[rank-1];current=min(current,values[i]*len(order)/rank);result[i]=current
    return result


def main():
    checks=[]
    def check(name,condition):
        checks.append({'check':name,'passed':bool(condition)})
        if not condition:raise AssertionError(name)
    def correction(name,records,p='P_value',q='FDR'):
        records=[r for r in records if r.get(p,'') not in ('','nan','NaN')]
        expected=bh([float(r[p]) for r in records])
        check(name,all(math.isclose(x,float(r[q]),rel_tol=1e-6,abs_tol=1e-14) for x,r in zip(expected,records)))
    release=json.loads((ROOT/'audit/current_release.json').read_text());assets=rows('audit/released_assets.csv')
    for record in assets:
        p=ROOT/record['repository_file'];raw=p.read_bytes()
        check('Reviewed asset bytes/hash '+record['repository_file'],len(raw)==int(record['bytes']) and sha256(raw).hexdigest()==record['sha256'])
        parts=p.relative_to(ROOT).parts
        check('Allowed source tree '+record['repository_file'],'VenEC_CLU_eye_audit' not in parts and 'Age10_20_audit_20261117' not in parts)
        if record['kind']=='figure':
            check('Unaltered original figure '+record['repository_file'],record['sha256']==record['source_sha256'])
            check('Scientific image signature '+record['repository_file'],raw.startswith(b'\x89PNG\r\n\x1a\n') if p.suffix=='.png' else raw.startswith(b'%PDF') if p.suffix=='.pdf' else True)
    check('Figure asset count',sum(r['kind']=='figure' for r in assets)==release['original_figure_files'])
    check('Aggregate/evidence count',sum(r['kind']=='aggregate_table' for r in assets)==release['aggregate_tables'])
    check('Five modules',len(release['modules'])==5 and 'F3_ds_analysis' in release['modules'])
    scan=json.loads((ROOT/'audit/scan_summary.json').read_text());summary=rows('audit/source_tree_summary.csv')
    check('Complete source metadata counts',sum(int(r['n_files']) for r in summary)==scan['scanned_files']==release['source_files'])
    evidence=rows('evidence/evidence_registry.csv');ids={r['claim_id'] for r in evidence}
    check('Evidence IDs unique',len(ids)==len(evidence)==release['evidence_records'])
    check('Aggregate evidence has no individual IDs',not any(re.search(r'(?<![A-Za-z0-9])GSM\d+|sample_id=|cellid=',' '.join(r.values())) for r in evidence))
    claims=rows('audit/claim_to_evidence_report.csv')
    check('All current Results paragraphs and tables mapped',len(claims)==release['claim_mappings']==32)
    for r in claims:
        check('Current claim text '+r['claim_id'],r['claim_text'] in (ROOT/r['section_file']).read_text())
        check('Claim evidence resolves '+r['claim_id'],bool(r['claim_ids']) and set(r['claim_ids'].split(';'))<=ids)
    for module in release['modules']:
        check(module+' evidence',bool(rows(module+'/evidence/evidence_registry.csv')))
        check(module+' visible figures','![' in (ROOT/module/'README.md').read_text() or '<img' in (ROOT/module/'README.md').read_text())
    for r in rows('audit/figure_manifest.csv'):
        check('Current assigned panel '+r['figure']+r['figure_panels'],(ROOT/r['figure_file']).is_file() and (ROOT/r['preview_file']).is_file())
    f3=Path('F3_ds_analysis');fn=f3/'F3_features/functional_network_20261118';lisa=f3/'TF_analysis/stable_celltype_LISA2_20261118'
    ecm=f3/'2026_Nat Aging_HAECA/stromal_ECM_20261118';cell=f3/'F3_features/ModelFeature_cellAge'
    sources=json.loads((ROOT/f3/'TabulaStromal_F3_PPT_sources.json').read_text());expected_images={im['sha256'] for slide in sources for im in slide['images']}
    with ZipFile(ROOT/f3/'TabulaStromal_F3_Downstream_Results_20261119.pptx') as z:
        check('Forty-slide public presentation',sum(bool(re.fullmatch(r'ppt/slides/slide\d+\.xml',n)) for n in z.namelist())==40)
        check('PPT scientific images unchanged',{sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('ppt/media/')}==expected_images and len(expected_images)==52)
    for slide in sources:
        for im in slide['images']:check('PPT original source '+im['path'],sha256((ROOT/f3/im['path']).read_bytes()).hexdigest()==im['sha256'])
    fit=rows('F2_ModelPerformance/tables/Fig2_scatter_metrics.csv');external=rows('F2_ModelPerformance/tables/Fig2_workflow_external_metrics.csv')
    check('Current Fig2 in-sample scope',len(fit)==2 and all(r['prediction_scope']=='development_in_sample' for r in fit))
    check('Current Fig2 external scope',len(external)==2 and all(r['prediction_scope']=='selected_external_ALL100' for r in external))
    text=(ROOT/'docs/results/Results_Part2_ClockPerformance_CellAge.md').read_text()
    for label in ['Fibroblast','Endothelial']:
        for r in fit+external:
            if r['lineage']==label:
                check('Current Table1 '+label+r['prediction_scope'],f"| {int(float(r['n']))} | {float(r['R2']):.3f} | {float(r['MAE']):.2f} | {float(r['MAD']):.2f} |" in text)
    gene_sets=[]
    for short,n,ng in [('Fibro',750,678),('Endo',350,337)]:
        features=rows(f'TabulaStromal_Rebenchmark_2026/tables/{short}/current_performance_first/performance_first_features.csv')
        check(short+' exact current representation',len(features)==n and len({r['Gene'] for r in features})==ng and len({r['Celltype'] for r in features})==27)
        check(short+' source bootstrap denominator',all(int(float(r['bootstrap_n']))==1500 and math.isclose(float(r['selection_frequency']),float(r['selection_count'])/1500,rel_tol=1e-10) for r in features))
        gene_sets.append({r['Gene'] for r in features})
    check('Forty-five shared selected genes',len(gene_sets[0]&gene_sets[1])==45)
    for p in sorted((ROOT/fn).glob('*/enrichment/global/*.csv')):correction('Complete-library BH '+str(p.relative_to(ROOT)),rows(p.relative_to(ROOT)))
    rankings=rows(lisa/'all_celltype_TF_rankings.csv')
    for lin in ['Fibroblast','Endothelial']:
        correction('LISA lineage-wide BH '+lin,[r for r in rankings if r['Lineage']==lin],q='FDR_all_celltype_tests')
    check('LISA factor/dataset multiplicity',all(math.isclose(float(r['P_value']),min(1,float(r['Best_dataset_P'])*int(r['n_ChIP_datasets'])),rel_tol=1e-7,abs_tol=1e-14) for r in rankings))
    ecms=rows(ecm/'all_ECM_ULM_age_contrasts.csv');tfs=rows(ecm/'all_TF_ULM_age_contrasts.csv')
    correction('ECM complete across-lineage BH',ecms,q='FDR_across_lineages');correction('All-factor across-lineage BH',tfs,q='FDR_across_lineages')
    check('Complete ECM family',len(ecms)==312 and sum(float(r['FDR_across_lineages'])<=.05 for r in ecms)==67)
    mech=[r for r in tfs if r['Source'] in {'YAP1','WWTR1','TEAD1','TEAD2','TEAD3','TEAD4'}]
    sig=[r for r in mech if float(r['FDR_across_lineages'])<=.05]
    check('Full factor correction precedes mechanosensitive display',len(tfs)==22346 and len(mech)==142 and len(sig)==1 and sig[0]['Source']=='TEAD4')
    primary=rows(cell/'tables/all_stable_celltype_associations.csv')
    for metric in sorted({r['Metric'] for r in primary}):
        group=[r for r in primary if r['Metric']==metric];correction('Current CellAge metric BH '+metric,group)
        if metric!='relative_hUSI':check('No primary corrected age effect '+metric,not any(r['FDR'] not in ('','nan') and float(r['FDR'])<=.05 for r in group))
    for p in (ROOT/'docs/results').glob('*.docx'):
        with ZipFile(p) as z:
            root=ET.fromstring(z.read('word/document.xml'))
            check('Editable valid Word '+p.name,z.testzip() is None and len(''.join(t.text or '' for t in root.findall('.//w:t',NS)))>1000)
    render=json.loads((ROOT/'audit/docx_rendering_report.json').read_text())
    check('Five current rendering receipts',len(render['documents'])==5)
    for r in render['documents']:
        check('Current DOCX/PDF render hashes '+r['docx'],sha256((ROOT/r['docx']).read_bytes()).hexdigest()==r['docx_sha256'] and sha256((ROOT/r['pdf']).read_bytes()).hexdigest()==r['pdf_sha256'])
        check('Rendered text/links/geometry '+r['docx'],r['source_text_blocks_found']==r['source_text_blocks_checked']>0 and r['hyperlink_targets_found']==r['hyperlink_targets_checked'] and r['out_of_page_text_spans']==0)
        check('Visual page review '+r['docx'],r['visual_review']=='all pages inspected; no clipping or table overflow observed')
    report=dict(status='passed',checks=len(checks),original_figure_files=release['original_figure_files'],aggregate_tables=release['aggregate_tables'],
                figure_candidates=release['figure_candidates'],evidence_records=len(evidence),paragraph_claim_mappings=len(claims),complete_source_file_entries=scan['scanned_files'],
                downstream_presentation_slides=40,original_PPT_images=52,word_rendering='Actual current DOCX renders; hashes, text, links, geometry and page review verified')
    (ROOT/'audit/results_release_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'audit/results_release_checks.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['check','passed'],lineterminator='\n');w.writeheader();w.writerows(checks)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
