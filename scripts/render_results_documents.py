#!/usr/bin/env python3
"""Render the actual DOCX files with LibreOffice and audit the resulting PDFs.

Requires a local LibreOffice executable and requirements/manuscript.txt.
Preview PNGs stay in ignored outputs/; PDF reading copies and hash receipts can
be committed. This does not regenerate, alter or validate scientific figures.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unicodedata

from docx import Document
from docx.opc.constants import RELATIONSHIP_TYPE as RT
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/results'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def canonical(text):
    text = unicodedata.normalize('NFKC', text)
    return ''.join(c for c in text if not c.isspace() and c not in '\u00ad\u200b\ufeff')


def audit_document(docx_path, pdf_path, preview_dir):
    source = Document(docx_path)
    section = source.sections[0]
    page_width = section.page_width.pt
    page_height = section.page_height.pt
    expected = [p.text for p in source.paragraphs if p.text.strip()]
    for table in source.tables:
        expected.extend(cell.text for row in table.rows for cell in row.cells if cell.text.strip())
    targets = {rel.target_ref for rel in source.part.rels.values() if rel.reltype == RT.HYPERLINK}
    body = []
    spatial_body = []
    pages = []
    fonts = set()
    links = set()
    bad_bounds = []
    embedded_fonts = True
    document_preview = preview_dir / docx_path.stem
    document_preview.mkdir(parents=True, exist_ok=True)
    with pymupdf.open(pdf_path) as pdf:
        for index, page in enumerate(pdf):
            if abs(page.rect.width-page_width) > 1 or abs(page.rect.height-page_height) > 1:
                raise AssertionError(f'Unexpected page dimensions: {docx_path.name}, page {index+1}')
            page_body = []
            span_count = 0
            for block in page.get_text('dict')['blocks']:
                for line in block.get('lines', []):
                    texts = []
                    for span in line['spans']:
                        x0,y0,x1,y1 = span['bbox']
                        if (x0 < -.5 or y0 < -.5 or x1 > page.rect.width+.5 or y1 > page.rect.height+.5):
                            bad_bounds.append({'page':index+1,'text':span['text'],'bbox':span['bbox']})
                        # Exclude page-number footer and retain all body text,
                        # including paragraphs split across a page boundary.
                        if y0 >= section.top_margin.pt-3 and y1 <= page_height-section.bottom_margin.pt+3:
                            texts.append(span['text']);span_count += 1
                    if texts:page_body.append(''.join(texts))
            if not canonical(''.join(page_body)):
                raise AssertionError(f'Blank body page: {docx_path.name}, page {index+1}')
            body.extend(page_body)
            body_rect=pymupdf.Rect(0,section.top_margin.pt-3,page_width,page_height-section.bottom_margin.pt+3)
            spatial_body.append(page.get_text('text',sort=True,clip=body_rect))
            links.update(link['uri'] for link in page.get_links() if link.get('uri'))
            fonts.update(font[3] for font in page.get_fonts())
            embedded_fonts = embedded_fonts and all(bool(pdf.extract_font(font[0])[3]) for font in page.get_fonts())
            image_name=f'page_{index+1:02d}.png'
            page.get_pixmap(matrix=pymupdf.Matrix(1.25,1.25),alpha=False).save(document_preview/image_name)
            pages.append({'page':index+1,'body_spans':span_count,'preview':image_name,
                          'width_pt':page.rect.width,'height_pt':page.rect.height})
    rendered = [canonical(' '.join(body)),canonical(' '.join(spatial_body))]
    # Font fallback glyphs (e.g. superscript minus) may be emitted after a
    # surrounding run in PDF object order; geometric reading order restores
    # them. Both streams are kept because tables need their native cell order.
    missing = [text for text in expected if not any(canonical(text) in stream for stream in rendered)]
    if missing:
        raise AssertionError(f'{docx_path.name}: {len(missing)} source text blocks missing: {missing[:2]}')
    if bad_bounds:raise AssertionError(f'Out-of-page text: {bad_bounds[:2]}')
    if not embedded_fonts:raise AssertionError('PDF has a non-embedded font: '+pdf_path.name)
    if targets-links:raise AssertionError(f'Missing PDF hyperlinks: {sorted(targets-links)}')
    return {
        'docx':str(docx_path.relative_to(ROOT)),'docx_sha256':digest(docx_path),
        'pdf':str(pdf_path.relative_to(ROOT)),'pdf_sha256':digest(pdf_path),
        'page_count':len(pages),'source_text_blocks_checked':len(expected),
        'source_text_blocks_found':len(expected)-len(missing),
        'hyperlink_targets_checked':len(targets),'hyperlink_targets_found':len(targets & links),
        'out_of_page_text_spans':len(bad_bounds),'embedded_pdf_fonts':sorted(fonts),
        'all_pdf_fonts_embedded':embedded_fonts,
        'pages':pages,'visual_review':'pending'
    }


def write_report(report):
    (ROOT/'audit/docx_rendering_report.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Actual Word rendering check','',
           'The released DOCX files were opened by LibreOffice and exported to PDF. '+
           'These are renders of the Word files, not separately typeset Markdown.', '',
           f'Renderer: **{report["renderer"]}**. PDF inspection: PyMuPDF {pymupdf.VersionBind}.', '',
           '| Document | Pages | Text blocks | Links | Page bounds | Visual review |',
           '| --- | ---: | ---: | ---: | --- | --- |']
    for record in report['documents']:
        name=Path(record['docx']).name
        lines.append(f'| [{name}](../{record["docx"]}) | {record["page_count"]} | '+
                     f'{record["source_text_blocks_found"]}/{record["source_text_blocks_checked"]} | '+
                     f'{record["hyperlink_targets_found"]}/{record["hyperlink_targets_checked"]} | passed | '+
                     f'{record["visual_review"]} |')
    lines += ['', 'Every nonempty source paragraph and table cell was found in extracted body text after '+
              'Unicode/whitespace normalization, using native and geometric PDF reading order. '+
              'Every external hyperlink target survived PDF export. '+
              'Pages were checked for blank bodies, expected dimensions and text outside the page.', '',
              'Times New Roman was requested by the DOCX. This server used the embedded PDF fonts listed '+
              'in the JSON receipt; the available Roman substitute was Liberation Serif. '+
              'Pagination can differ on a machine with other fonts or a different Word processor.', '',
              'The JSON receipt binds each check to the exact DOCX/PDF SHA256. Page PNGs are kept locally '+
              'under ignored outputs/document_qa/. A passed geometry/text check alone is not visual approval.', '']
    (ROOT/'audit/docx_rendering_report.md').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--soffice', default=shutil.which('soffice') or shutil.which('libreoffice'))
    args=parser.parse_args()
    if not args.soffice:parser.error('Install LibreOffice or pass --soffice /path/to/soffice')
    executable=str(Path(args.soffice).resolve())
    env=os.environ.copy();env['SAL_USE_VCLPLUGIN']='svp'
    version=subprocess.run([executable,'--version'],env=env,check=True,capture_output=True,text=True).stdout.strip()
    documents=sorted(DOCS.glob('*.docx'))
    if len(documents)!=5:raise AssertionError('Expected three Results files, combined Results and legends')
    records=[]
    with tempfile.TemporaryDirectory(prefix='tabulastromal-word-render-') as temp:
        temp=Path(temp);(temp/'pdf').mkdir()
        command=[executable,'-env:UserInstallation='+(temp/'profile').as_uri(),
                 '--headless','--convert-to','pdf:writer_pdf_Export','--outdir',str(temp/'pdf')]
        completed=subprocess.run(command+[str(path) for path in documents],env=env,
                                 capture_output=True,text=True,timeout=180)
        if completed.returncode:
            raise RuntimeError(f'LibreOffice exit {completed.returncode}: {completed.stdout} {completed.stderr}')
        for docx in documents:
            generated=temp/'pdf'/docx.with_suffix('.pdf').name
            if not generated.is_file():raise RuntimeError('No PDF generated for '+docx.name+': '+completed.stderr)
            pdf=docx.with_suffix('.pdf');shutil.copyfile(generated,pdf)
            records.append(audit_document(docx,pdf,ROOT/'outputs/document_qa'))
    report={'status':'passed_automated_render_checks_visual_review_pending',
            'renderer':version,'pdf_inspector':pymupdf.VersionBind,
            'original_scientific_figures_modified':False,'documents':records}
    write_report(report)
    print(json.dumps({'documents':len(records),'pages':sum(r['page_count'] for r in records),
                      'status':report['status']},indent=2))


if __name__=='__main__':main()
