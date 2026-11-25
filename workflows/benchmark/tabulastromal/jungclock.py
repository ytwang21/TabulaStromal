"""JungClock-TabulaStromal-adapted: process selection + native H2O GLM.

Uses the audited 2026 author scripts, not the historical ProcessRidge model.
Already-normalized Celltype_Gene inputs retain subtype identity. TPM and fsva
are not applied: raw bulk counts and Gene/Dox/LIBPrep covariates are unavailable.
The Java H2O backend is accessed through its official Python API. All process
selection and grouped inner CV see outer-training samples only.
"""
from pathlib import Path
import importlib.util, os, sys, socket
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from .storage import write_csv, write_json, dump


def audit_method(ctx):
    if 'JungClock' not in ctx.cfg['models']['names']:return []
    import subprocess,shutil,json,importlib.metadata
    from .storage import write_text,file_hash
    root=ctx.root/'audit';source=Path(ctx.cfg['jungclock'].get('reference_dir',str(Path(__file__).resolve().parents[1]/'resources/plesa_clock')))
    provenance=json.loads((source/'method_provenance.json').read_text())
    sources={p.name:file_hash(p) for p in source.glob('*.R')}
    rscript=shutil.which('Rscript') or 'external_inputs/Rscript'
    check=subprocess.run([rscript,'-e','cat(R.version.string,"\\n");for(p in c("h2o","sva","DESeq2"))cat(p,requireNamespace(p,quietly=TRUE),"\\n")'],capture_output=True,text=True)
    write_text(root/'JungClock_dependency_check.md',f'''# JungClock dependency inspection

Rscript: {rscript}
{check.stdout}
{check.stderr}
Java: {shutil.which('java')}
H2O official Python package: {ctx.cfg['jungclock']['python_vendor']}/h2o (3.44.0.3).
Backend: native Java H2O GLM via official Python API, not an sklearn Ridge approximation.
Per worker JVM heap: {ctx.cfg['jungclock']['max_mem_size']}; one JVM per worker, separate local ports.
Official archive attempted: https://zenodo.org/records/23570472/files/Datasets_PlesaEtAl2026_NatAging.zip ; HTTP 403 on 2026-11-23.
Ontology fallback: existing public MBotC annotation CSV with provenance in resources/plesa_clock/ontology_source.json.
Only the ontology resource originates in the author's other repository; no MultiTIMER model is used.
Actual toy training/native save/reload/age-inverse evidence: tests/jungclock_smoke_evidence.json.
''')
    write_text(root/'JungClock_method_audit.md',f'''# JungClock-TabulaStromal-adapted

Official reference: https://github.com/saschajung/Clock ; reviewed commit {provenance['current_code_commit']}.
Publication: https://doi.org/10.1038/s43587-026-01234-8 . Requested bioRxiv v1 access failed.
Official R scripts inspected and hashed: {json.dumps(sources,indent=2)}.
No repository software LICENSE was present in the inspected source inventory; no license is invented.

| Official executed step | Current implementation and reason |
|---|---|
| ProcessFibroblastData.R: bulk count datasets | Supplied sample/subtype means of stored expression; current Age labels, no imported paper training ages |
| TPM with measured gene lengths then log1p | Omitted for normalized supplied matrices; no invented gene lengths or second TPM |
| DESeq2/SVA/fsva and Gene/Dox/LIBPrep design | Unsupported original covariates/assay; replaced by outer-train-only existing scaler, explicitly not numerical reproduction of fsva |
| SelectProcesses.R: level-3 MBotC, per-process PCA + H2O ridge GLM, process Lasso | Candidate genes mapped from explicit schema; each Celltype_Gene remains an individual predictor; train-only PCA, H2O GLM and native process Lasso. Require >=10 distinct genes up front. Filter processes using grouped inner CV R2>=0 rather than optimistic train R2. Outer Valid/external never enter |
| SelectProcesses.R custom age scaling at70 | Process selection uses the final model Horvath transform consistently; a declared difference from the original process screen |
| TrainClock.R seven manually fixed priorProcs | Learn processes from each outer train using the reviewed selector; do not relabel fixed historical processes as data-learned |
| Horvath target transform at20 and inverse | Implemented exactly; predictions in years, including negative values |
| Computed age weights | Reviewed TrainClock.R actually passes a Weight column of ones; retain unit weights, not its unused computed vector |
| H2O GLM alpha=0, lambda_search, standardize, IRLSM | Native H2O backend retained; grouped inner fold_column replaces random sample folds. Final GLM with >4000 pair predictors uses native COORDINATE_DESCENT, same ridge objective/lambda search, to avoid quadratic IRLSM Gram memory |
| LibType predictor | Omitted: not a Celltype_Gene predictor and unavailable consistently |
| ApplyNewData.R native load and inverse | Native H2O save/load, explicit selected columns and inverse age transform |

K is the train-only allowed pair universe, not the final number of model predictors. Actual process/feature counts are exported per fold.
This is a defined native-H2O adaptation; it is neither a full original bulk-data reproduction nor the historical Plesa2023_ProcessRidge baseline.
External adaptation for this model is explicitly residual calibration, not GLM continuation.
''')
    return [root/'JungClock_dependency_check.md',root/'JungClock_method_audit.md']


def backend(cfg):
    vendor=cfg.get('jungclock',{}).get('python_vendor')
    if vendor and vendor not in sys.path:sys.path.insert(0,vendor)
    import h2o
    if h2o.connection() is None:
        # One JVM per worker, separate from the other lineage and its workers.
        with socket.socket() as s:s.bind(('127.0.0.1',0));port=s.getsockname()[1]
        h2o.init(ip='127.0.0.1',port=port,start_h2o=True,name=f'stromal_{os.getpid()}',
                 max_mem_size=cfg['jungclock'].get('max_mem_size','8G'),
                 nthreads=cfg['resources']['threads'],bind_to_localhost=True,
                 ice_root=cfg['jungclock'].get('ice_root'),verbose=False)
        h2o.no_progress()
    return h2o


def transform_age(y):
    y=np.asarray(y,float)
    return np.where(y<=20,np.log1p(y)-np.log(21),(y-20)/21)


def inverse_age(z):
    z=np.asarray(z,float);out=np.empty_like(z);mask=z<=0
    out[mask]=np.exp(z[mask]+np.log(21))-1;out[~mask]=21*z[~mask]+20
    return out


class JungClock:
    protocol='JungClock-TabulaStromal-adapted; native H2O; no TPM/fsva'
    def __init__(self,cfg,seed=42,feature_names=None,checkpoint=None):
        self.cfg=cfg;self.seed=seed;self.feature_names=list(feature_names or [])

    def fit(self,X,y,metadata=None):
        h2o=backend(self.cfg)
        from h2o.estimators.glm import H2OGeneralizedLinearEstimator
        if len(self.feature_names)!=X.shape[1]:raise ValueError('JungClock requires exact feature schema')
        if metadata is None or 'group_id' not in metadata:raise ValueError('Grouped inner CV metadata required')
        from .training import grouped_splits
        n=min(self.cfg['cv']['inner_splits'],metadata.group_id.nunique())
        if n<2:raise ValueError('JungClock needs >=2 independent inner groups')
        folds=np.empty(len(y),int)
        for i,(_,va) in enumerate(grouped_splits(metadata.reset_index(drop=True),n,self.seed)):folds[va]=i
        ont=pd.read_csv(self.cfg['jungclock']['ontology'])
        ont=ont[(ont.ProcessLevel==3)&(ont.ProcessName!='Background genes')]
        # Mapping is supplied explicitly from schema, never split subtype names.
        mapping=self.cfg['jungclock'].get('feature_gene_map',{})
        if not mapping:raise ValueError('Missing explicit Celltype_Gene -> Gene mapping')
        genes=np.array([mapping[f] for f in self.feature_names])
        self.process_manifest=[];columns={};ty=transform_age(y)
        for process,tab in ont.groupby('ProcessName',sort=True):
            ix=np.flatnonzero(np.isin(genes,tab.Symbol.unique()))
            if len(set(genes[ix]))<self.cfg['jungclock'].get('min_genes',10):continue
            z=PCA(n_components=min(len(y)-1,len(ix)),svd_solver='full').fit_transform(X[:,ix])
            frame=h2o.H2OFrame(pd.DataFrame(z,columns=[f'PC{i}' for i in range(z.shape[1])]).assign(target=ty,inner_fold=folds))
            model=H2OGeneralizedLinearEstimator(family='gaussian',alpha=0,lambda_search=True,
                       standardize=True,seed=self.seed,solver='IRLSM',max_iterations=1000)
            model.train(x=[f'PC{i}' for i in range(z.shape[1])],y='target',fold_column='inner_fold',training_frame=frame)
            r2=float(model.r2(xval=True));selected=bool(np.isfinite(r2) and r2>=0)
            key=f'process_{len(self.process_manifest)}'
            if selected:columns[key]=model.predict(frame).as_data_frame().iloc[:,0].to_numpy()
            self.process_manifest.append(dict(process=process,key=key,inner_R2=r2,eligible_for_stack=selected,
                n_pairs=len(ix),n_genes=len(set(genes[ix])),indices=ix.tolist()))
            h2o.remove(model);h2o.remove(frame)
        if not columns:raise RuntimeError('No train-only process passes inner CV R2>=0; no Ridge fallback')
        frame=h2o.H2OFrame(pd.DataFrame(columns).assign(target=ty,inner_fold=folds))
        selector=H2OGeneralizedLinearEstimator(family='gaussian',alpha=1,lambda_search=True,standardize=True,seed=self.seed+1,solver='IRLSM',max_iterations=1000)
        selector.train(x=list(columns),y='target',fold_column='inner_fold',training_frame=frame)
        active={k for k,v in selector.coef().items() if k!='Intercept' and v!=0}
        for item in self.process_manifest:item['selected']=item['key'] in active
        self.indices=np.array(sorted({j for r in self.process_manifest if r['selected'] for j in r['indices']}),int)
        h2o.remove(selector);h2o.remove(frame)
        if not len(self.indices):raise RuntimeError('Process-selection Lasso selected no process; no fallback')
        frame=h2o.H2OFrame(pd.DataFrame(X[:,self.indices],columns=[f'f{i}' for i in self.indices]).assign(target=ty,inner_fold=folds,Weight=1.))
        self.final_solver='IRLSM' if len(self.indices)<=self.cfg['jungclock'].get('max_irlsm_features',4000) else 'COORDINATE_DESCENT'
        self.model=H2OGeneralizedLinearEstimator(family='gaussian',link='identity',alpha=0,
            lambda_search=True,standardize=True,seed=self.seed,solver=self.final_solver,max_iterations=1000)
        self.model.train(x=[f'f{i}' for i in self.indices],y='target',fold_column='inner_fold',weights_column='Weight',training_frame=frame)
        h2o.remove(frame);return self

    def predict(self,X,metadata=None):
        h2o=backend(self.cfg);frame=h2o.H2OFrame(pd.DataFrame(X[:,self.indices],columns=[f'f{i}' for i in self.indices]))
        result=self.model.predict(frame);pred=inverse_age(result.as_data_frame().iloc[:,0].to_numpy())
        h2o.remove(result);h2o.remove(frame);return pred

    def save(self,path):
        path=Path(path).resolve();path.mkdir(parents=True,exist_ok=True);h2o=backend(self.cfg)
        native=Path(h2o.save_model(self.model,path=str(path),force=True))
        state=dict(native_model=str(native),indices=self.indices.tolist(),feature_names=self.feature_names,
                   seed=self.seed,protocol=self.protocol,final_solver=self.final_solver,process_manifest=self.process_manifest)
        write_json(path/'jungclock.json',state)
        write_csv(path/'process_selection.csv',pd.DataFrame(self.process_manifest))
        write_csv(path/'process_feature_manifest.csv',pd.DataFrame([dict(process=r['process'],feature=self.feature_names[j],selected=r['selected']) for r in self.process_manifest for j in r['indices']]))
        return [native,path/'jungclock.json',path/'process_selection.csv',path/'process_feature_manifest.csv']

    @classmethod
    def load(cls,path,cfg):
        import json
        state=json.loads((Path(path)/'jungclock.json').read_text());obj=cls(cfg,state['seed'],state['feature_names'])
        obj.indices=np.array(state['indices'],int);obj.process_manifest=state['process_manifest'];obj.final_solver=state.get('final_solver','IRLSM')
        obj.model=backend(cfg).load_model(state['native_model']);return obj
