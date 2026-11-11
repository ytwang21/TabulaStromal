#!/usr/bin/env python
"""Validate completed biological outputs without rereading atlas expression matrices."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def validate(root):
    root=Path(root);manifest=json.loads((root/'run_manifest.json').read_text())
    checks=[]
    def check(condition,label):
        checks.append({'check':label,'passed':bool(condition)})
        if not condition:raise AssertionError(label)
    try:
        for stage in ['audit','atlas','markers','trajectory']:
            for lin in ['fibro','endo']:
                check(manifest['stages'].get(f'{stage}:{lin}',{}).get('status')=='complete',f'{stage}:{lin} complete')
        check(manifest['stages'].get('figure:shared',{}).get('status')=='complete','integrated figure stage complete')
        for lin in ['fibro','endo']:
            p=root/lin;info=json.loads((root/'shared'/f'data_audit_{lin}.json').read_text())
            for name in ['summary_tissue.csv','summary_tissue_top.csv','summary_age_group.csv']:
                t=pd.read_csv(p/name)
                check(t.n_cells.sum()==info['n_cells'],f'{lin} {name} cell totals')
                check(t.n_samples.between(0,info['n_samples']).all(),f'{lin} {name} sample bounds')
            t=pd.read_csv(p/'marker_summary.csv')
            check(t.fraction_expressing.between(0,1).all(),f'{lin} marker fraction bounds')
            check(np.isfinite(t.mean_expression).all() and (t.mean_expression>=0).all(),f'{lin} marker expression finite/nonnegative')
            for state,count in info['categories']['general'].items():
                check((t.loc[t.subtype==state,'n_cells']==count).all(),f'{lin} marker denominator {state}')
            scores=np.load(p/'signature_scores.npy',mmap_mode='r')
            check(scores.shape==(info['n_cells'],len(info['categories']['general'])),f'{lin} signature shape')
            check(np.isfinite(scores).all(),f'{lin} all signature scores finite')
            if manifest['config'].get('trajectory_engine')=='metacell_three_method':
                from metacell_trajectory import model_paths,load_model
                obs=pd.read_pickle(p/'trajectory/prepared/metadata.pkl')
                audit=json.loads((p/'trajectory/metacell_audit.json').read_text())
                check(audit['n_meta_cells']>=len(obs)>0,f'{lin} explicit meta-cell source/cohort counts')
                check(obs.index.is_unique,f'{lin} unique fitted meta-cell identity')
                paths=model_paths(p/'trajectory')
                families={path.parent.name if path.name!='cellrank' else 'cellrank' for path in paths}
                check(families=={'monocle','stream','cellrank'},f'{lin} all three requested packages completed')
                for path in paths:
                    raw_ids=pd.read_csv(path/'cells.csv',usecols=['cell_id']).cell_id.astype(str)
                    check(raw_ids.is_unique and set(raw_ids)<=set(obs.index),f'{lin} native fitted IDs belong to common cohort: {path.name}')
                    if path.name!='cellrank':
                        check(set(raw_ids)==set(obs.index),f'{lin} complete native fitted cohort: {path.name}')
                    f=load_model(path,obs)
                    check(np.array_equal(f.index,obs.index),f'{lin} aligned meta-cell IDs: {path.name}')
                    finite=np.isfinite(f.pseudotime)
                    check(finite.any() and f.loc[finite,'pseudotime'].between(0,1.00001).all(),f'{lin} valid pseudotime: {path.name}')
                from scipy import sparse
                t=sparse.load_npz(p/'trajectory/cellrank/transition_matrix.npz')
                check((t.data>=0).all() and np.allclose(np.asarray(t.sum(axis=1)).ravel(),1,atol=1e-6),f'{lin} real stochastic CellRank transition matrix')
                fate=pd.read_csv(p/'trajectory/cellrank/fate_probabilities.csv.gz',index_col=0).to_numpy()
                check(np.isfinite(fate).all() and (fate>=-1e-5).all() and (fate<=1.00001).all() and np.allclose(fate.sum(axis=1),1,atol=1e-4),f'{lin} solved CellRank fate probability simplex')
                decision=json.loads((p/'trajectory/trajectory_selection.json').read_text())
                check(decision['method'] in ['monocle','stream','cellrank'],f'{lin} main panel uses requested method, no PAGA')
                continue
            rows=np.load(p/'trajectory_rows.npy');pt=pd.read_csv(p/'trajectory'/'pseudotime.csv.gz')
            check(np.array_equal(pt.atlas_row,rows),f'{lin} trajectory row alignment')
            check(pt.atlas_row.is_unique and pt.cell_id.is_unique,f'{lin} unique sampled cell identity')
            finite=np.isfinite(pt.pseudotime)
            check(finite.any() and pt.loc[finite,'pseudotime'].between(0,1.00001).all(),f'{lin} DPT range')
            check(pt.loc[finite,'component'].nunique()==1,f'{lin} DPT restricted to root component')
            for f in (p/'trajectory').glob('paga_*.csv'):
                if f.name=='paga_edges.csv':continue
                a=pd.read_csv(f,index_col=0).to_numpy()
                check(np.allclose(a,a.T) and np.isfinite(a).all() and (a>=0).all() and (a<=1.00001).all(),f'{lin} PAGA matrix {f.name}')
            expected=[f'{lin}_umap_celltype',f'{lin}_umap_subtype_dashed',f'{lin}_umap_general_dashed',
                      f'{lin}_dataset_summary',f'{lin}_marker_dotplot',f'{lin}_signature_heatmap',f'{lin}_signature_sample_boxplots']
            for stem in expected:
                for suffix in ['png','pdf']:
                    f=p/f'{stem}.{suffix}';check(f.exists() and f.stat().st_size>1000,f'{lin} {stem}.{suffix}')
            decision=json.loads((p/'trajectory'/'trajectory_selection.json').read_text())
            check(decision['main_panel'] in ['paga','dpt'],f'{lin} explicit main-panel decision')
        for suffix in ['png','pdf']:
            f=root/'integrated_figure'/f'Fig_StromalAtlas_Overview.{suffix}'
            check(f.exists() and f.stat().st_size>1000,f'integrated figure {suffix}')
        for lin,info in manifest['inputs'].items():
            stat=Path(info['path']).stat()
            check(stat.st_size==info['size'] and stat.st_mtime_ns==info['mtime_ns'],f'{lin} original input size/mtime unchanged')
        status='passed'
    except Exception:
        (root/'validation_report.json').write_text(json.dumps(dict(status='failed',checks=checks),indent=2))
        raise
    result=dict(status=status,n_checks=len(checks),checks=checks,
                limits='Computational consistency checks do not validate biological lineage direction, annotation independence, normalization provenance or publication suitability.')
    (root/'validation_report.json').write_text(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run_dir',type=Path)
    result=validate(parser.parse_args().run_dir)
    print(f'{result["status"]}: {result["n_checks"]} checks')
