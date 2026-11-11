"""Shared, execution-only helpers for the final ElasticNet clock.

Historical benchmark code and objects are read-only. New artifacts live in the
existing lineage directories, using final_clock names rather than new folders.
"""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import copy,json,logging,os,sys,warnings
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from tabulastromal.storage import file_hash,digest,write_json,write_csv,write_text,dump
from audit_external_metric_scope import independent_metrics

ROOT=Path(__file__).resolve().parent
# Use installed dependencies; checkpoint locations are explicit and overridable.
RESULT_ROOT=Path(os.environ.get('TABULASTROMAL_BENCHMARK_ROOT',str(ROOT.parents[1]/'outputs/benchmark'))).resolve()
ROOTS={lin:Path(os.environ.get('TABULASTROMAL_'+lin.upper()+'_ROOT',str(RESULT_ROOT/lin))).resolve()
       for lin in ['Fibro','Endo']}
METRICS=['R2','MAE','MAD','Age_MAD','MedAE','RMSE','PCC']
STAGES=['audit','bootstrap','final_model','external','enrichment','figures']
LOG=logging.getLogger('final_clock')

def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()

class Context:
    def __init__(self,lineage,config=None,resume=True):
        self.lineage=lineage;self.root=ROOTS[lineage];self.resume=resume
        self.options=yaml.safe_load(Path(config or ROOT/f'final_clock_config_{lineage}.yaml').read_text())
        self.identity=json.loads((self.root/'reports/run_identity.json').read_text())
        self.cfg=self.identity['config'];self.source_signature=self.identity['signature']
        self.K=int(json.loads((self.root/'stability/selection_decision.json').read_text())['EN_K'])
        scientific={k:v for k,v in self.options.items() if k not in ['workers','io_workers','dpi','reserve_memory_gb']}
        self.signature=digest(dict(protocol='FinalEN-Phase3-4-v1',lineage=lineage,source=self.source_signature,options=scientific,K=self.K))
    def path(self,rel):return self.root/rel
    def csv(self,rel):return pd.read_csv(self.path(rel),keep_default_na=False,low_memory=False)
    def json(self,rel):return json.loads(self.path(rel).read_text())

def dataset(ctx,name='train'):
    p=ctx.path('features/'+name)
    return dict(path=p,X=sp.load_npz(p/'pseudobulk.npz'),meta=pd.read_csv(p/'samples.csv',keep_default_na=False),
                schema=json.loads((p/'schema.json').read_text()))

def score(y,p):
    y=np.asarray(y,float);p=np.asarray(p,float)
    if not len(y):return dict(n=0,**{k:np.nan for k in METRICS},metric_note='no_eligible_samples')
    s=independent_metrics(y,p)
    s['metric_note']='insufficient_n_or_constant_age' if len(y)<2 or np.ptp(y)==0 else 'constant_prediction' if np.ptp(p)==0 else ''
    return s

def source_valid(ctx,marker):
    try:
        r=json.loads(Path(marker).read_text())
        return r['status']=='complete' and r['signature']==ctx.source_signature and all(
            Path(x['path']).is_file() and Path(x['path']).stat().st_size==x['size'] and file_hash(x['path'])==x['sha256'] for x in r['outputs'])
    except (OSError,KeyError,ValueError):return False

def own_valid(ctx,marker,code=None):
    try:
        r=json.loads(Path(marker).read_text())
        return ctx.resume and r['signature']==ctx.signature and r['status']=='complete' and (
            code is None or r['code_sha256']==file_hash(code)) and all(
            Path(x['path']).is_file() and Path(x['path']).stat().st_size==x['size'] and file_hash(x['path'])==x['sha256'] for x in r['outputs'])
    except (OSError,KeyError,ValueError):return False

def mark(ctx,marker,outputs,code):
    paths=list(dict.fromkeys(str(Path(p).resolve()) for p in outputs))
    write_json(marker,dict(signature=ctx.signature,status='complete',completed_at=now(),code_sha256=file_hash(code),
        common_sha256=file_hash(__file__),outputs=[dict(path=p,size=Path(p).stat().st_size,sha256=file_hash(p)) for p in paths]))

def fit_checked(model,X,y):
    model=copy.deepcopy(model)
    params=model.get_params();changes={'warm_start':False} if 'warm_start' in params else {}
    model.set_params(**changes)
    X=np.asarray(X,dtype=np.float64,order='F')
    y=np.asarray(y,float);initial=model.max_iter;notes=[]
    for factor in [1,5,25]:
        model.set_params(max_iter=int(initial*factor),**({'warm_start':factor>1} if 'warm_start' in params else {}))
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always');model.fit(X,y)
        notes += [str(w.message) for w in ws]
        bad=any(issubclass(w.category,ConvergenceWarning) for w in ws)
        gap=float(np.max(np.asarray(model.dual_gap_)))
        if not bad and np.isfinite(gap) and gap<=model.tol*max(float(np.var(y)),1e-15)*1.05+1e-8:
            model.final_clock_fit_warnings_=notes;return model
    raise RuntimeError(f'ElasticNet convergence failed: gap={gap}; warnings={notes[-2:]}')

def versions(ctx):return ctx.json('models/model_versions.json')

def load_final_clock(model_path):
    """Load a self-contained full-development bundle, never a best-fold proxy."""
    p=Path(model_path)
    if p.is_dir():p=p/'models/final_clock.joblib' if (p/'models').is_dir() else p/'final_clock.joblib'
    b=joblib.load(p)
    if b.get('format')!='TabulaStromal-FinalEN-v1':raise ValueError('Not a final ElasticNet clock bundle')
    if list(b['feature_schema']['selected_features'])!=list(b['selected_features']):raise ValueError('Feature-schema mismatch')
    return b

def bundle_predict(b,X):
    # Column selection can give equal-valued C/F-order pandas buffers. Use
    # one memory order so float32 BLAS accumulation does not depend on that
    # representation. Frozen scaling, feature values and model stay intact.
    z=b['preprocessor'].transform(np.ascontiguousarray(X,dtype=np.float32))
    return np.asarray(b['model'].predict(np.ascontiguousarray(z)),float)

def selected_h5ad(path,b,expression='X'):
    """Chunked sample×selected-Celltype_Gene means, with NO ratio features.

    Ages are optional for inference. Unknown subtypes are excluded as in the
    source cache. Absent subtype/genes are zeros under the recorded old policy.
    """
    import anndata as ad
    from tabulastromal.data import canonical_sample,matrix_source
    a=ad.read_h5ad(path,backed='r')
    try:
        f=b['fields'];obs=a.obs
        for key in ['sample','gse','celltype']:
            if f[key] not in obs:raise ValueError('Missing inference field '+f[key])
        ids=obs[f['gse']].astype(str)+'|'+obs[f['sample']].astype(str).map(canonical_sample)
        ct=obs[f['celltype']].astype(str);known=set(b['feature_schema']['celltypes'])
        valid=ct.isin(known)&~obs[f['sample']].astype(str).str.lower().isin(['','nan','none','unknown'])
        counts=ids[valid].value_counts();samples=sorted(counts[counts>=b['min_cells_per_sample']].index)
        si=pd.Index(samples).get_indexer(ids);valid=valid.to_numpy()&(si>=0)
        source,gn=matrix_source(a,expression);gn=pd.Index(gn.astype(str))
        if not gn.is_unique:raise ValueError('Duplicate gene symbols in inference input')
        mappings={}
        for cell,rows in pd.DataFrame(b['feature_schema']['selected_mapping']).groupby('Celltype',sort=False):
            gi=gn.get_indexer(rows.Gene);good=gi>=0
            mappings[cell]=(gi[good],rows.column_index.to_numpy(int)[good])
        sums=np.zeros((len(samples),len(b['selected_features'])),np.float32)
        nct={c:np.zeros(len(samples),np.int64) for c in mappings}
        for start in range(0,len(obs),b['chunk_cells']):
            stop=min(start+b['chunk_cells'],len(obs));local=np.flatnonzero(valid[start:stop])
            if not len(local):continue
            chunk=sp.csr_matrix(source[start:stop],dtype=np.float32)[local]
            if not np.isfinite(chunk.data).all():raise ValueError('Nonfinite inference expression')
            cell=ct.iloc[start+local].to_numpy();sample=si[start+local]
            for c in set(cell)&set(mappings):
                mask=np.flatnonzero(cell==c);g,col=mappings[c]
                nct[c]+=np.bincount(sample[mask],minlength=len(samples))
                if not len(g):continue
                membership=sp.csr_matrix((np.ones(len(mask),np.float32),(sample[mask],np.arange(len(mask)))),shape=(len(samples),len(mask)))
                sums[:,col]+=(membership@chunk[mask][:,g]).toarray()
        for c,(_,col) in mappings.items():sums[:,col]/=np.maximum(nct[c],1)[:,None]
        meta=pd.DataFrame(dict(sample_id=samples,n_cells=counts.reindex(samples).to_numpy(),GSE=[s.split('|',1)[0] for s in samples]))
        age_key=next((c for c in [f['age']]+f.get('age_fallbacks',[]) if c in obs),None)
        if age_key:
            ages=pd.DataFrame(dict(sample_id=ids,Age=pd.to_numeric(obs[age_key],errors='coerce'))).groupby('sample_id').Age.first()
            meta['age_true']=ages.reindex(samples).to_numpy()
        return sums,meta
    finally:a.file.close()

def predict_age(h5ad_or_pseudobulk,clock,expression='X'):
    """Predict from .h5ad, cached feature directory, or named raw pseudobulk DataFrame."""
    b=load_final_clock(clock) if not isinstance(clock,dict) else clock
    if isinstance(h5ad_or_pseudobulk,pd.DataFrame):
        frame=h5ad_or_pseudobulk
        if frame.columns.duplicated().any():raise ValueError('Duplicate input feature names')
        X=frame.reindex(columns=b['selected_features'],fill_value=0).to_numpy(np.float32)
        meta=pd.DataFrame(dict(sample_id=frame.index.astype(str)))
    else:
        p=Path(h5ad_or_pseudobulk)
        if p.suffix=='.h5ad':X,meta=selected_h5ad(p,b,expression)
        else:
            schema=json.loads((p/'schema.json').read_text());ix=pd.Index(schema['features']).get_indexer(b['selected_features'])
            raw=sp.load_npz(p/'pseudobulk.npz');X=np.zeros((raw.shape[0],len(ix)),np.float32);good=ix>=0
            X[:,good]=raw[:,ix[good]].toarray();meta=pd.read_csv(p/'samples.csv',keep_default_na=False)
    if not np.isfinite(X).all():raise ValueError('Nonfinite pseudobulk input; no silent imputation')
    out=meta.copy();out['age_pred']=bundle_predict(b,X);out['final_model_id']=b['model_id']
    if 'Age' in out:out['age_true']=pd.to_numeric(out.Age,errors='coerce')
    if 'age_true' in out:out['age_acceleration_raw']=out.age_pred-out.age_true
    return out

def make_bundle(ctx,model,scaler,indices,schema,model_id):
    ng=len(schema['genes']);mapping=[dict(feature=schema['features'][j],Celltype=schema['celltypes'][j//ng],
        Gene=schema['genes'][j%ng],column_index=i,source_column_index=int(j)) for i,j in enumerate(indices)]
    selected=[r['feature'] for r in mapping]
    fs=dict(selected_features=selected,selected_mapping=mapping,celltypes=schema['celltypes'],genes=schema['genes'],
        missing_policy='zero for unobserved subtype/gene BEFORE frozen source scaler; no proportions',age_unit='years')
    return dict(format='TabulaStromal-FinalEN-v1',model=model,preprocessor=scaler,selected_features=selected,
        feature_schema=fs,model_id=model_id,lineage=ctx.lineage,K=ctx.K,source_signature=ctx.source_signature,
        training_scope='full_development',fields=ctx.cfg['fields'],min_cells_per_sample=ctx.cfg['features']['min_cells_per_sample'],
        chunk_cells=ctx.cfg['features']['chunk_cells'],created_at=now())
