#!/usr/bin/env python3
"""Bundle code, executed notebook, final figures, group statistics and QC reports.

The several-GB cache and full-cell Parquet scores remain in the server directory.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json
import time


def main():
    r=Path(__file__).resolve().parents[1]
    selected=[r/name for name in ['config.yaml','README.md','requirements.txt','run_nohup.sh',
                                  'StromalAtlas_ExtendFig1_SignatureScore.md']]
    for folder in ['scripts','notebooks','marker_sets','audit','reports']:
        selected += [p for p in (r/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    selected += [p for p in (r/'data').glob('*') if p.is_file() and p.suffix in {'.csv','.json'}]
    selected += [p for p in (r/'figures/integrated').glob('ExtFig1_SignatureScore_Fibro_Endo.*') if p.suffix in {'.png','.pdf'}]
    for lineage in ['Fibroblast','Endothelial']:
        for basename in [f'ExtFig1_SignatureScore_{lineage}',f'{lineage}_Leiden_concordance',
                         f'{lineage}_candidate_versions',f'{lineage}_reference_Leiden_UMAP']:
            selected += [r/f'figures/panels/{basename}{ext}' for ext in ['.png','.pdf']]
    selected += [p for p in (r/'logs').glob('*.json') if p.is_file()]
    selected += [r/'logs/pipeline.log']
    selected=sorted(set(selected));assert all(p.exists() for p in selected)
    manifest=[{'file':str(p.relative_to(r)),'bytes':p.stat().st_size,
               'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in selected]
    mp=r/'DELIVERY_MANIFEST.json'
    mp.write_text(json.dumps({'created_at':time.strftime('%Y-%m-%d %H:%M:%S'),
        'excludes':['data/cache/','data/cell_signature_scores_*.parquet'],
        'reason':'Full cell scores and cached gene/score arrays are preserved on server; compact package includes all group statistics, figures and reproducible code.',
        'files':manifest},ensure_ascii=False,indent=2))
    archive=r/'ExtendedFig1_SignatureScore_Deliverables.zip'
    with ZipFile(archive,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for p in selected+[mp]:z.write(p,arcname='extended_figure1_SS/'+str(p.relative_to(r)))
    with ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps({'archive':str(archive),'size_MB':round(archive.stat().st_size/1e6,2),'n_files':len(selected)+1}))


if __name__=='__main__':main()
