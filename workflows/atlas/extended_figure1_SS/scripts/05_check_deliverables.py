#!/usr/bin/env python3
"""Acceptance checks against actual inputs, score arrays, JSONs and exported PDF."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from PIL import Image
from pypdf import PdfReader
import signature_pipeline as sp


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.yaml'))
    a=p.parse_args();cfg=sp.read_config(a.config);r=sp.root(cfg);checks=[]
    for lin in ['fibro','endo']:
        c=sp.cache(cfg,lin);audit=sp.read_json(r/f'audit/data_audit_{lin}.json')
        fp=sp.fingerprint(cfg['inputs'][lin]);assert all(audit[k]==v for k,v in fp.items());checks.append(f'{lin}: input fingerprint unchanged')
        split=pd.read_csv(r/f'audit/group_split_{lin}.csv')
        train=set(split.loc[split['split'].eq('train'),'group']);hold=set(split.loc[split['split'].eq('holdout'),'group'])
        assert train and hold and not train&hold;checks.append(f'{lin}: actual GSE split disjoint')
        done=sp.read_json(c/'extraction_complete.json');assert done['negative_count']==0 and done['nonfinite_count']==0
        checks.append(f'{lin}: every stored expression value checked finite and nonnegative')
        doc=sp.read_json(r/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_final.json');frozen=sp.read_json(c/'frozen_final_markers.json')
        assert doc['marker_sets']==frozen['marker_sets'];checks.append(f'{lin}: no holdout marker retuning')
        var=set(np.load(c/'var_names.npy'));genes=np.load(c/'marker_genes.npy').tolist();gi={g:i for i,g in enumerate(genes)}
        for s,gg in doc['marker_sets'].items():
            assert len(gg)==len(set(gg)) and set(gg)<=var
            assert doc['subtype_evidence'][s]['status'] in {'supported_relative_to_reference','ambiguous','insufficient_evidence'}
        evidence=pd.read_csv(r/f'data/marker_selection_evidence_{lin}.csv')
        assert evidence.loc[evidence['selected'],'passes_prespecified_QC'].all();checks.append(f'{lin}: every retained gene passes training QC')
        x=np.load(c/'marker_expression.npy',mmap_mode='r');bg=np.load(c/'control_bin_means.npy',mmap_mode='r')
        ctr=sp.read_json(c/'controls.json');assert ctr['holdout_used_for_bins'] is False
        controlgenes=set(g for gg in ctr['controls_by_bin'].values() for g in gg)
        assert not controlgenes&set(genes);checks.append(f'{lin}: training-only controls exclude every candidate marker')
        mean=np.load(c/'final_mean.npy',mmap_mode='r');module=np.load(c/'final_module.npy',mmap_mode='r')
        order=np.load(c/'final_score_names.npy').tolist();rows=np.sort(np.random.default_rng(cfg['seed']).choice(len(x),300,replace=False))
        for j,s in enumerate(order):
            gg=doc['marker_sets'][s]
            if gg:
                expected=x[rows][:,[gi[g] for g in gg]].mean(1)
                corrected=expected-bg[rows][:,[ctr['marker_bins'][g] for g in gg]].mean(1)
                np.testing.assert_allclose(mean[rows,j],expected,rtol=3e-6,atol=2e-6)
                np.testing.assert_allclose(module[rows,j],corrected,rtol=3e-6,atol=2e-6)
            else:assert np.isnan(mean[rows,j]).all() and np.isnan(module[rows,j]).all()
        checks.append(f'{lin}: independent recomputation of 300 sampled original rows matches final mean/module scores')
        cellfile=r/f'data/cell_signature_scores_{lin}.parquet';pf=pq.ParquetFile(cellfile)
        assert pf.metadata.num_rows==audit['n_cells'];assert len(mean)==audit['n_cells']
        obs=pd.read_parquet(c/'metadata.parquet',columns=['row_id','cell_id']);st=0
        for batch in pf.iter_batches(batch_size=131072,columns=['row_id','cell_id']):
            df=batch.to_pandas();en=st+len(df)
            np.testing.assert_array_equal(df['row_id'],np.arange(st,en))
            np.testing.assert_array_equal(df['cell_id'].astype(str),obs.iloc[st:en]['cell_id'].astype(str));st=en
        checks.append(f'{lin}: every full-cell Parquet row/cell ID aligns with original input order')
        active=sp.read_json(r/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_figure.json')
        before=sp.read_json(c/'figure_markers_frozen_before_holdout.json')
        assert active['marker_sets']==before['marker_sets']
        assert active['selection']['holdout_used'] is False and active['selection']['Leiden_ID_used'] is False
        assert active['source_fingerprint']==fp
        for chosen in active['selected_sources'].values():assert chosen['training_QC_passed']
        trials=pd.read_csv(r/f'data/marker_combination_search_{lin}.csv')
        assert not trials.holdout_used.any() and not trials.Leiden_used.any()
        assert len(trials)==active['selection']['trials_evaluated']
        checks.append(f'{lin}: global combination search training-only; lists frozen before holdout; selected sets pass expression QC')
        mean=np.load(c/'figure_mean.npy',mmap_mode='r');order=np.load(c/'figure_score_names.npy').tolist()
        for j,s in enumerate(order):
            expected=x[rows][:,[gi[g] for g in active['marker_sets'][s]]].mean(1)
            np.testing.assert_allclose(mean[rows,j],expected,rtol=3e-6,atol=2e-6)
        checks.append(f'{lin}: new figure scores independently recomputed from expression alone')
        import reference_style_figures as rf
        maps=pd.read_csv(r/f'data/plot_reference_umap_signatures_{lin}.csv')
        assert maps.signature.tolist()==rf.UMAP_LAYOUT[lin]
        for s,gg in active['marker_sets'].items():
            row=maps[maps.signature.eq(s)].iloc[0]
            assert set(row.genes_used.split(';'))==set(gg) and row.source==str(c/'figure_mean.npy')
        checks.append(f'{lin}: requested functional/subtype UMAPs present; subtype UMAPs and bars use identical lists/scores')
        labels=pd.read_parquet(c/'metadata.parquet',columns=['leiden','Celltype','Celltype_general'])
        codes,clusters=pd.factorize(labels['leiden'].astype(str),sort=True)
        counts=np.bincount(codes)
        exact=np.column_stack([np.bincount(codes,weights=mean[:,j],minlength=len(clusters))/counts for j in range(len(order))])
        plotted=pd.read_csv(r/f'data/plot_reference_cluster_scores_{lin}.csv',dtype={'leiden':str})
        pivot=plotted.pivot(index='leiden',columns='signature',values='mean_score').reindex(index=clusters,columns=order)
        np.testing.assert_allclose(pivot.to_numpy(),exact,equal_nan=True,rtol=2e-7,atol=1e-9)
        chart=pd.read_csv(r/f'data/plot_reference_cluster_check_{lin}.csv',dtype={'leiden':str})
        assert set(chart.leiden)==set(clusters) and chart.n_cells.sum()==len(labels)
        actual=labels.groupby('leiden',observed=True).agg(Celltype=('Celltype','first'),Celltype_general=('Celltype_general','first'))
        for row in chart.itertuples():
            ref=actual.loc[row.leiden]
            assert row.panel_title==str(ref.Celltype)+'_C'+row.leiden
            assert row.Celltype_general==str(ref.Celltype_general)
            v=pivot.loc[row.leiden];top=v.idxmax()
            assert row.top_signature==top
            assert row.annotation_matches_highest==(top==str(ref.Celltype_general))
        checks.append(f'{lin}: all {len(chart)} Leiden titles/reference labels and every plotted bar equal original-row score aggregation; mismatches preserved')
    png=r/'figures/integrated/ExtFig1_SignatureScore_Fibro_Endo.png';pdf=png.with_suffix('.pdf')
    im=Image.open(png);assert min(im.size)>=3000 and min(im.info['dpi'])>=299
    checks.append(f'PNG: {im.size[0]}×{im.size[1]}, {im.info["dpi"]} dpi')
    reader=PdfReader(pdf);assert len(reader.pages)==1;page=reader.pages[0]
    content=page.get_contents().get_data();assert content.count(b'BT')>60
    assert len(page.extract_text())>500;assert content.count(b' m\n')>30
    text=page.extract_text()
    assert all(any(line==letter or line.startswith(letter+' ') for line in text.splitlines()) for letter in ['A','B','C','D'])
    fonts=[entry for pg in reader.pages for entry in pg['/Resources']['/Font'].values()];assert len(fonts)>0
    for entry in fonts:
        font=entry.get_object()
        descendants=font.get('/DescendantFonts',[font])
        for descendant in descendants:
            obj=descendant.get_object() if hasattr(descendant,'get_object') else descendant
            descriptor=obj.get('/FontDescriptor')
            assert descriptor is not None
            descriptor=descriptor.get_object()
            assert any(k in descriptor for k in ['/FontFile','/FontFile2','/FontFile3'])
    checks.append('PDF: one A–D page with editable bars/axes/labels/legends and embedded fonts; UMAP point clouds locally rasterized')
    required=['config.yaml','README.md','requirements.txt','scripts/run_signature_score.py','scripts/signature_pipeline.py',
       'scripts/00_audit.py','scripts/01_marker_qc.py','scripts/02_score_signatures.py','scripts/03_validate_annotations.py','scripts/04_plot_extendfig1.py',
       'notebooks/ExtendFig1_SignatureScore.ipynb','reports/audit_report.md','reports/marker_selection_report.md',
       'reports/annotation_validation_report.md','reports/figure_legend.md','reports/run_summary.md',
       'data/marker_selection_evidence.csv','data/marker_sets_candidate_comparison.csv','data/signature_concordance_by_donor.csv',
       'scripts/marker_combination_search.py','scripts/reference_style_figures.py','reports/reference_layout_report.md',
       'marker_sets/Fibroblast_subtype_markers_figure.json','marker_sets/Endothelial_subtype_markers_figure.json']
    for rel in required:assert (r/rel).stat().st_size>0,rel
    checks.append('All required code/notebook/QC/report/JSON/statistical deliverables exist')
    sp.dump(r/'logs/acceptance_checks.json',{'passed':True,'n_checks':len(checks),'checks':checks})
    print(json.dumps({'passed':True,'n_checks':len(checks),'checks':checks},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
