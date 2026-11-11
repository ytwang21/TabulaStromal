"""Compare the plot-only reader with standard pickle.load on real AnnData."""
import os
from pathlib import Path
import json
import pickle
import tempfile

HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'tmp'/'matplotlib'))
import anndata as ad
import numpy as np
import pandas as pd
from pickle_metadata import MetadataUnpickler,normalize_records


def main():
    rng=np.random.default_rng(42)
    obs=pd.DataFrame({'Celltype_curated':pd.Categorical(['Fibroblast','Other']*128),
                      'Fib_marker':pd.Categorical(['Fibroblast','None']*128)},
                     index=[f'cell_{i}' for i in range(256)])
    a=ad.AnnData(rng.normal(size=(256,512)),obs=obs,dtype=np.float64)
    a.obsm['X_umap']=rng.normal(size=(256,2))
    a.obsm['X_pca']=rng.normal(size=(256,50))
    a.raw=a.copy()
    evaluation=pd.DataFrame({'Marker':['Marker','DEG(only)'],'Recall':[.5,.8]})
    source=[{'adata':a,'recall_df':evaluation}]
    (HERE/'tmp').mkdir(exist_ok=True)
    with tempfile.TemporaryFile(dir=HERE/'tmp') as f:
        pickle.dump(source,f,protocol=4)
        f.seek(0);full=pickle.load(f)
        f.seek(0);loader=MetadataUnpickler(f,threshold=4096);light=loader.load()
        pd.testing.assert_frame_equal(full[0]['recall_df'],light[0]['recall_df'])
        pd.testing.assert_frame_equal(full[0]['adata'].obs,light[0]['adata'].obs)
        np.testing.assert_array_equal(full[0]['adata'].obsm['X_umap'],loader.array(light[0]['adata'].obsm['X_umap']))
        assert loader.skipped_bytes>0 and len(loader.deferred_arrays)>0
        assert len(normalize_records(pd.DataFrame([{'recall_df':evaluation,'adata':a}])))==1
        assert len(normalize_records({'dataset_1':{'recall_df':evaluation,'adata':a}}))==1
        report={'status':'passed','checks':['Exact UMAP equality with standard pickle.load',
            'Exact obs and categorical-label equality','Exact evaluation DataFrame equality',
            'Large expression/raw/PCA buffers deferred','DataFrame and dataset-dictionary adapters'],
            'fixture_skipped_bytes':loader.skipped_bytes,'fixture_deferred_arrays':len(loader.deferred_arrays)}
    out=HERE/'results';out.mkdir(exist_ok=True)
    (out/'loader_validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
