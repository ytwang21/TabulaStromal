#!/usr/bin/env python3
"""Build editable Results Word files from the same Markdown displayed on GitHub."""
import json
from pathlib import Path
import re
from zipfile import ZipFile
from xml.etree import ElementTree as ET

from docx import Document
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / 'docs/results'


def plain(text):
    text = re.sub(r'!?\[([^\]]+)\]\([^)]*\)', r'\1', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'\1', text)
    return text.replace('`', '')


def configured():
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = sec.bottom_margin = Inches(.85)
    sec.left_margin = sec.right_margin = Inches(.9)
    normal = doc.styles['Normal']
    normal.font.name = 'Times New Roman'
    normal.font.size = Pt(11)
    normal.paragraph_format.line_spacing = 1.3
    normal.paragraph_format.space_after = Pt(7)
    return doc


def append_markdown(doc, text):
    expected, tables = [], []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1; continue
        if line.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                cells = [plain(c.strip()) for c in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[:\- ]+', c) for c in cells):
                    rows.append(cells)
                i += 1
            t = doc.add_table(rows=1, cols=len(rows[0]))
            t.style = 'Table Grid'
            for index, row in enumerate(rows):
                cells = t.rows[0].cells if index == 0 else t.add_row().cells
                for cell, value in zip(cells, row):cell.text=value
            tables.append(rows)
            continue
        if line.startswith('#'):
            level = len(line)-len(line.lstrip('#'))
            value = plain(line[level:].strip())
            doc.add_heading(value, level=min(level,3))
            expected.append(value)
            i += 1; continue
        paragraph = [line]
        i += 1
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|')):
            paragraph.append(lines[i].strip());i += 1
        value=plain(' '.join(paragraph))
        doc.add_paragraph(value)
        expected.append(value)
    return expected,tables


def verify(path, expected, expected_tables):
    reopened=Document(path)
    actual=[p.text for p in reopened.paragraphs if p.text]
    assert actual==expected, path.name+' text mismatch'
    tables=[[[c.text for c in row.cells] for row in t.rows] for t in reopened.tables]
    assert tables==expected_tables, path.name+' table mismatch'
    with ZipFile(path) as archive:
        assert archive.testzip() is None
        for name in archive.namelist():
            if name.endswith('.xml'):ET.fromstring(archive.read(name))
    return dict(file=str(path.relative_to(ROOT)),text_paragraphs=len(actual),tables=len(tables),
                zip_and_xml='passed',source_text_and_tables='identical',editable_text=True,
                pagination_rendering='not performed: no Word/LibreOffice renderer available')


def main():
    files=sorted(DIR.glob('Results_Part*.md'))
    assert len(files)==3
    reports=[]
    combined=configured();all_text=[];all_tables=[]
    for index,path in enumerate(files):
        text=path.read_text()
        doc=configured();expected,tables=append_markdown(doc,text)
        target=path.with_suffix('.docx');doc.save(target)
        reports.append(verify(target,expected,tables))
        if index:combined.add_page_break()
        expected,tables=append_markdown(combined,text)
        all_text+=expected;all_tables+=tables
    target=DIR/'TabulaStromal_Results_Combined.docx';combined.save(target)
    reports.append(verify(target,all_text,all_tables))
    legend=DIR/'Figure_Legends.md'
    if legend.exists():
        doc=configured();expected,tables=append_markdown(doc,legend.read_text())
        target=legend.with_suffix('.docx');doc.save(target)
        reports.append(verify(target,expected,tables))
    report=dict(status='passed_text_table_and_package_checks',documents=reports,
                limitation='Page layout, font substitution and clipping require a real Word/PDF rendering pass.')
    (ROOT/'audit/docx_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'audit/docx_validation_report.md').write_text(
        '# Word validation\n\n'+
        '\n'.join(f'- {r["file"]}: editable text, source paragraphs/tables and ZIP/XML verified.' for r in reports)+
        '\n\nNo Word, LibreOffice or PDF page renderer is installed. Pagination, font substitution and clipping have not been visually validated. This is an editable Results draft, not a claim of publication-ready typesetting.\n')
    print(json.dumps({'validated_documents':len(reports),'editable_text':'verified','pagination_rendering':'not available'},indent=2))


if __name__=='__main__':main()
