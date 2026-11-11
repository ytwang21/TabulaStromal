#!/usr/bin/env python3
"""Verify figure assets, evidence links, Results tables and editable Word content."""
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
from zipfile import ZipFile
from xml.etree import ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}


def rows(rel):
    with (ROOT/rel).open() as handle:return list(csv.DictReader(handle))


def main():
    checks=[]
    def check(name,condition):
        checks.append({'check':name,'passed':bool(condition)})
        if not condition:raise AssertionError(name)
    assets=rows('audit/released_assets.csv')
    figures=rows('audit/figure_inventory.csv')
    evidence=rows('evidence/evidence_registry.csv')
    identifiers={r['claim_id'] for r in evidence}
    for record in assets:
        p=ROOT/record['repository_file'];raw=p.read_bytes()
        check('Asset bytes/hash: '+record['repository_file'],len(raw)==int(record['bytes']) and sha256(raw).hexdigest()==record['sha256'])
        if record['kind']=='figure':
            check('Original unaltered figure: '+record['repository_file'],record['sha256']==record['source_sha256'])
            check('Figure file signature: '+record['repository_file'],
                  raw.startswith(b'\x89PNG\r\n\x1a\n') if p.suffix=='.png' else raw.startswith(b'%PDF'))
    check('117 original figure files',sum(r['kind']=='figure' for r in assets)==117)
    check('81 released aggregate tables',sum(r['kind']=='aggregate_table' for r in assets)==81)
    check('847 figure candidates',len(figures)==847)
    check('Published figure mapping',sum(r['published']=='True' for r in figures)==117)
    scan=json.loads((ROOT/'audit/scan_summary.json').read_text())
    summary=rows('audit/source_tree_summary.csv')
    check('Complete source summary count',sum(int(r['n_files']) for r in summary)==scan['scanned_files'])
    for claim in rows('audit/claim_to_evidence_report.csv'):
        text=(ROOT/claim['section_file']).read_text()
        check('Claim exists in section: '+claim['heading'],claim['claim_text'] in text)
        check('Evidence IDs resolve: '+claim['heading'],bool(claim['claim_ids']) and set(claim['claim_ids'].split(';'))<=identifiers)
    for module in ['F1_StromalAtlas','F2_ModelPerformance','TabulaStromal_Rebenchmark_2026','CellAge_Exploration']:
        check(module+' module evidence exists',len(rows(module+'/evidence/evidence_registry.csv'))>0)
        check(module+' figures displayed','![' in (ROOT/module/'README.md').read_text() or '<img' in (ROOT/module/'README.md').read_text())
    text=(ROOT/'docs/results/Results_Part2_ClockPerformance_CellAge.md').read_text()
    table=[line for line in text.splitlines() if line.startswith('| Fibro |') or line.startswith('| Endo |')]
    expected=[]
    for lin in ['Fibro','Endo']:
        for r in rows(f'TabulaStromal_Rebenchmark_2026/tables/{lin}/external/final_clock_external_protocol_comparison.csv'):
            expected.append((lin,int(float(r['n'])),f"{float(r['R2']):.3f}",f"{float(r['MAE']):.2f}"))
    actual=[]
    for line in table:
        values=[v.strip() for v in line.strip('|').split('|')]
        actual.append((values[0],int(values[2]),values[3],values[4]))
    check('All eight Results external-table rows match source',actual==expected)
    for p in (ROOT/'docs/results').glob('*.docx'):
        with ZipFile(p) as archive:
            xml=ET.fromstring(archive.read('word/document.xml'))
            texts=[t.text or '' for t in xml.findall('.//w:t',NS)]
            check('Editable Word text: '+p.name,len(''.join(texts))>1000 and archive.testzip() is None)
    report={'status':'passed','checks':len(checks),'original_figure_files':117,'aggregate_tables':81,
            'figure_candidates':len(figures),'evidence_records':len(evidence),
            'paragraph_claim_mappings':len(rows('audit/claim_to_evidence_report.csv')),
            'complete_source_file_entries':scan['scanned_files'],
            'word_rendering':'Not performed; text/table/XML validation is recorded separately.'}
    (ROOT/'audit/results_release_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'audit/results_release_checks.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=['check','passed'],lineterminator='\n')
        writer.writeheader();writer.writerows(checks)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
