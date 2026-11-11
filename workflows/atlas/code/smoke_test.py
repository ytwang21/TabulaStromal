#!/usr/bin/env python
"""End-to-end synthetic verification. Synthetic outputs never enter biological results."""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import stromal_atlas_pipeline as pipeline
from utils_data import MARKERS, ORDER, csr_rows


def main():
    from anndata import AnnData
    import h5py
    folder=Path(tempfile.mkdtemp(prefix='stromal_smoke_',dir=pipeline.HERE.parent/'tmp'))
    rng=np.random.default_rng(17);cfg=json.loads((pipeline.HERE/'config.json').read_text())
    cfg.update(chunk_rows=67,plot_cells=300,trajectory_cells=240,trajectory_neighbors=[8,12],
               pilot_blocks=4,pilot_block_rows=32,controls_per_bin=5,signature_bins=8,top_tissues=2,dpi=65)
    cfg['inputs']={}
    cfg['metacell_inputs']={}
    cfg['metacell_trajectory'].update(max_fit_metacells=240,cellrank_landmarks=40,monocle_hvg=[40],monocle_deg=[20])
    for lin in ['fibro','endo']:
        groups=ORDER[lin];genes=sorted({g for gs in MARKERS[lin].values() for g in gs})+[f'CONTROL{i:03d}' for i in range(150)]
        n=600;state=np.array(groups)[np.arange(n)%len(groups)]
        counts=rng.poisson(.25,(n,len(genes)))
        for k in groups:
            rows=np.flatnonzero(state==k);cols=[genes.index(g) for g in MARKERS[lin][k]]
            counts[np.ix_(rows,cols)]+=rng.poisson(4,(len(rows),len(cols)))
        x=np.log1p(counts).astype(np.float32)
        obs=pd.DataFrame(dict(Celltype_general=state,Celltype=[k+'_marker_pan' for k in state],
                              Celltype_sub=[k+'_marker' for k in state],
                              tissue=np.array(['Tissue1','Tissue2','Tissue3'])[np.arange(n)//len(groups)%3],
                              sample=[f'S{i//12:03d}' for i in range(n)],GSE=[f'GSE{i//60}' for i in range(n)],
                              Age_group=np.array(['A','B','C','Unknown'])[np.arange(n)//20%4],
                              age=['Unknown']*n),index=[f'{lin}_cell{i}' for i in range(n)])
        obs['n_cells_in_meta']=30
        for c in obs:
            if c!='n_cells_in_meta':obs[c]=obs[c].astype('category')
        ad=AnnData(x,obs=obs,var=pd.DataFrame({'highly_variable':True},index=genes),dtype=np.float32)
        # Deliberate disconnected state islands test NA outside root component.
        pca=rng.normal(0,.2,(n,12)).astype(np.float32)
        # Continuous PCA keeps the small Markov test graph connected.
        ad.obsm['X_pca']=pca
        ad.obsm['X_umap']=np.column_stack([np.array([groups.index(k)*5 for k in state])+rng.normal(0,.4,n),rng.normal(0,.7,n)]).astype(np.float32)
        ad.uns['log1p']={}
        import scanpy as sc
        sc.pp.neighbors(ad,n_neighbors=30,use_rep='X_pca',random_state=17)
        path=folder/f'{lin}.h5ad';ad.write_h5ad(path);cfg['inputs'][lin]=str(path);cfg['metacell_inputs'][lin]=str(path)
        with h5py.File(path,'r') as f:
            np.testing.assert_allclose(csr_rows(f['X'],11,24).toarray(),x[11:24])
    # Only the explicitly synthetic smoke run overrides production's fixed source guards.
    cfg['overview_plot'] = {'expected_metacell_inputs': cfg['metacell_inputs'], 'dpi': cfg['dpi']}
    config=folder/'config.json';config.write_text(json.dumps(cfg))
    run=folder/'run'
    subprocess.run([sys.executable,str(pipeline.HERE/'stromal_atlas_pipeline.py'),'--config',str(config),'--run-dir',str(run)],check=True)
    for lin in ['fibro','endo']:
        p=run/lin
        assert pd.read_csv(p/'summary_tissue.csv').n_cells.sum()==600
        assert pd.read_csv(p/'summary_tissue_top.csv').n_cells.sum()==600
        t=pd.read_csv(p/'marker_summary.csv')
        assert t.fraction_expressing.between(0,1).all()
        assert np.isfinite(t.mean_expression).all()
        assert np.load(p/'signature_scores.npy').shape==(600,len(ORDER[lin]))
        # Compare a real streamed result to direct matrix computation, including zeros.
        from utils_data import load_metadata
        obs,var,_=load_metadata(cfg['inputs'][lin],lin,cfg)
        with h5py.File(cfg['inputs'][lin],'r') as f:x=csr_rows(f['X'],0,600).toarray()
        controls=json.loads((p/'signature_controls.json').read_text())
        expected=[]
        for k in ORDER[lin]:
            ids=var.index.get_indexer(controls['markers'][k]);ctrl=var.index.get_indexer(controls['controls'][k])
            expected.append(x[:,ids].mean(axis=1)-x[:,ctrl].mean(axis=1))
        np.testing.assert_allclose(np.load(p/'signature_scores.npy'),np.array(expected).T,atol=1e-6)
        first=t.iloc[0];mask=obs.general==first.subtype;g=var.index.get_loc(first.gene)
        np.testing.assert_allclose(first.mean_expression,x[mask,g].mean(),rtol=1e-5)
        np.testing.assert_allclose(first.fraction_expressing,(x[mask,g]>0).mean())
        from metacell_trajectory import model_paths
        paths=model_paths(p/'trajectory')
        assert len(paths)==5  # 2 Monocle + 2 STREAM + 1 CellRank
        cohort=pd.read_pickle(p/'trajectory/prepared/metadata.pkl')
        for model in paths:
            ids=pd.read_csv(model/'cells.csv').cell_id.astype(str)
            assert set(ids)<=set(cohort.index)
        fate=pd.read_csv(p/'trajectory/cellrank/fate_probabilities.csv.gz',index_col=0).to_numpy()
        np.testing.assert_allclose(fate.sum(axis=1),1,atol=1e-4)
    assert (run/'integrated_figure'/'Fig_StromalAtlas_Overview.pdf').stat().st_size>1000
    subprocess.run([sys.executable,str(pipeline.HERE/'stromal_atlas_pipeline.py'),'--config',str(config),'--run-dir',str(run),'--resume'],check=True)
    result={'status':'passed','synthetic_only':True,'directory':str(folder),
            'checks':['bounded dense H5AD row access','full-metadata summary totals','all-stage rendering',
                      'streamed marker means/fractions versus direct matrix','fixed-control signature scores versus direct computation',
                      'all three real tools on identical meta-cell cohorts','genuine solved CellRank fate probability simplex','config-matched resume']}
    (pipeline.HERE.parent/'tmp'/'smoke_test_latest.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
