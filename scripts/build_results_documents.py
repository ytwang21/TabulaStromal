#!/usr/bin/env python3
"""Build editable Results Word files from the same Markdown displayed on GitHub."""
import json
from pathlib import Path
import re
import posixpath
from urllib.parse import quote, urlsplit
from zipfile import ZipFile
from xml.etree import ElementTree as ET

from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / 'docs/results'
REPOSITORY = 'https://github.com/ytwang21/TabulaStromal/blob/main/'
INLINE = re.compile(r'\[([^\]]+)\]\(([^)]+)\)|\*\*([^*]+)\*\*|`([^`]+)`')


def rich_text(paragraph, text, source):
    """Retain readable labels and working evidence/figure hyperlinks in Word."""
    offset = 0
    for match in INLINE.finditer(text):
        paragraph.add_run(text[offset:match.start()])
        label, target, bold, code = match.groups()
        if target is not None:
            parsed = urlsplit(target)
            if not parsed.scheme:
                relative = posixpath.normpath(str(source.parent.relative_to(ROOT)) + '/' + parsed.path)
                if relative.startswith('../') or relative.startswith('/'):
                    raise ValueError('Link escapes repository: ' + target)
                if not (ROOT / relative).is_file():
                    raise FileNotFoundError(relative)
                target = REPOSITORY + quote(relative, safe='/')
                if parsed.fragment:target += '#' + parsed.fragment
            elif parsed.scheme not in {'http', 'https'}:
                raise ValueError('Unsupported hyperlink scheme: ' + target)
            hyperlink = OxmlElement('w:hyperlink')
            hyperlink.set(qn('r:id'), paragraph.part.relate_to(target, RT.HYPERLINK, is_external=True))
            run = OxmlElement('w:r');properties = OxmlElement('w:rPr')
            color = OxmlElement('w:color');color.set(qn('w:val'), '155E8A')
            properties.append(color)
            underline = OxmlElement('w:u');underline.set(qn('w:val'), 'single')
            properties.append(underline);run.append(properties)
            value = OxmlElement('w:t');value.text = label
            run.append(value);hyperlink.append(run);paragraph._p.append(hyperlink)
        else:
            run = paragraph.add_run(bold if bold is not None else code)
            if bold is not None:run.bold = True
        offset = match.end()
    paragraph.add_run(text[offset:])


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
    normal.paragraph_format.widow_control = True
    for level, size in [(1,17), (2,13), (3,11)]:
        style=doc.styles[f'Heading {level}']
        style.font.name='Times New Roman';style.font.size=Pt(size)
        style.font.color.rgb=RGBColor.from_string('143D59')
        style.paragraph_format.keep_with_next=True
    for name in ['Normal','Heading 1','Heading 2','Heading 3']:
        fonts=doc.styles[name].element.get_or_add_rPr().get_or_add_rFonts()
        for attribute in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme']:
            fonts.attrib.pop(qn('w:'+attribute),None)
        for attribute in ['ascii','hAnsi','eastAsia','cs']:
            fonts.set(qn('w:'+attribute),'Times New Roman')
    doc.core_properties.title='TabulaStromal — evidence-linked Results'
    doc.core_properties.author='TabulaStromal project'
    footer=sec.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run('TabulaStromal  |  ')
    field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE')
    footer._p.append(field)
    return doc


def append_markdown(doc, text, source):
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
            t.autofit = False
            available=doc.sections[0].page_width-doc.sections[0].left_margin-doc.sections[0].right_margin
            weights=[.12,.46,.08,.14,.20] if len(rows[0])==5 else [1/len(rows[0])]*len(rows[0])
            for column,weight in zip(t.columns,weights):column.width=int(available*weight)
            for index, row in enumerate(rows):
                cells = t.rows[0].cells if index == 0 else t.add_row().cells
                properties=t.rows[index]._tr.get_or_add_trPr()
                properties.append(OxmlElement('w:cantSplit'))
                if index==0:properties.append(OxmlElement('w:tblHeader'))
                for cell,value,weight in zip(cells,row,weights):
                    cell.width=int(available*weight);cell.text=value
                    for paragraph in cell.paragraphs:
                        paragraph.paragraph_format.line_spacing=1.1
                        paragraph.paragraph_format.space_after=Pt(4)
                        for run in paragraph.runs:run.font.size=Pt(9);run.bold=(index==0)
            tables.append(rows)
            continue
        if line.startswith('#'):
            level = len(line)-len(line.lstrip('#'))
            value = plain(line[level:].strip())
            heading=doc.add_heading(level=min(level,3))
            rich_text(heading,line[level:].strip(),source)
            expected.append(value)
            i += 1; continue
        paragraph = [line]
        i += 1
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|')):
            paragraph.append(lines[i].strip());i += 1
        value=plain(' '.join(paragraph))
        rich_text(doc.add_paragraph(),' '.join(paragraph),source)
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
                hyperlinks=sum(r.reltype==RT.HYPERLINK for r in reopened.part.rels.values()),
                pagination_rendering='See separate docx_rendering_report.json; re-render after rebuilding.')


def main():
    files=sorted(DIR.glob('Results_Part*.md'))
    assert len(files)==3
    reports=[]
    combined=configured();all_text=[];all_tables=[]
    for index,path in enumerate(files):
        text=path.read_text()
        doc=configured();expected,tables=append_markdown(doc,text,path)
        target=path.with_suffix('.docx');doc.save(target)
        reports.append(verify(target,expected,tables))
        if index:combined.add_page_break()
        expected,tables=append_markdown(combined,text,path)
        all_text+=expected;all_tables+=tables
    target=DIR/'TabulaStromal_Results_Combined.docx';combined.save(target)
    reports.append(verify(target,all_text,all_tables))
    legend=DIR/'Figure_Legends.md'
    if legend.exists():
        doc=configured();expected,tables=append_markdown(doc,legend.read_text(),legend)
        target=legend.with_suffix('.docx');doc.save(target)
        reports.append(verify(target,expected,tables))
    report=dict(status='passed_text_table_and_package_checks',documents=reports,
                layout_report='audit/docx_rendering_report.json',
                limitation='Text/package checks alone do not verify pagination. A rendering receipt must match the current DOCX hashes.')
    (ROOT/'audit/docx_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'audit/docx_validation_report.md').write_text(
        '# Word validation\n\n'+
        '\n'.join(f'- {r["file"]}: editable text, source paragraphs/tables and ZIP/XML verified.' for r in reports)+
        '\n\nFigure and evidence links are clickable. Tables repeat their header on continuation pages and retain editable cells. '+
        'Page layout is checked separately in [the rendering report](docx_rendering_report.md); that report must match the current DOCX hashes.\n')
    print(json.dumps({'validated_documents':len(reports),'editable_text':'verified','pagination_rendering':'separate report required'},indent=2))


if __name__=='__main__':main()
