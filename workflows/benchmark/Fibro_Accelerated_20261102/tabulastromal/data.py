"""Read-only input inspection and chunked sparse sample/subtype mean expression."""
from pathlib import Path
import json, logging, re, platform, sys, importlib.metadata, shutil
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import psutil
from .storage import write_csv, write_json, write_text, dump, atomic_path, file_info
LOG = logging.getLogger(__name__)
NULLS = {'nan','none','unknown','unassigned','','NA'}


def canonical_sample(s):
    return re.sub(r'(_query)+$', '', str(s))


def metadata(obs, cfg):
    f=cfg['fields']; o=obs.copy()
    age_key=next((k for k in [f['age']]+f.get('age_fallbacks',[]) if k in o),None)
    if age_key is None: raise ValueError('No configured age target field is available')
    for key in ['sample','gse','celltype']:
        if f[key] not in o: raise ValueError(f'Missing required obs field {f[key]}')
    out=pd.DataFrame(index=o.index)
    out['sample']=o[f['sample']].astype(str).map(canonical_sample)
    out['GSE']=o[f['gse']].astype(str)
    out['GSM']=o[f['gsm']].astype(str) if f['gsm'] in o else ''
    out['sample_id']=out['GSE']+'|'+out['sample']
    # Keep biological sample rows; use conservative GSM groups when donor identity is absent.
    dk=f.get('donor')
    out['group_id']=(out['GSE']+'|donor:'+o[dk].astype(str) if dk else
                     out['GSE']+'|GSM:'+out['GSM'].where(~out['GSM'].isin(NULLS),out['sample']))
    out['Age']=pd.to_numeric(o[age_key],errors='coerce')
    out['Celltype']=o[f['celltype']].astype(str)
    out['valid']=np.isfinite(out.Age)&~out.Celltype.str.lower().isin(NULLS)&~out['sample'].str.lower().isin(NULLS)
    for c in ['tissue','sex','disease_global','Disease_Status','Disease_Global_Category']:
        if c in o: out[c]=o[c].astype(str)
    if out.loc[out.valid,'Age'].lt(0).any() or out.loc[out.valid,'Age'].gt(130).any():
        raise ValueError('Age outside [0,130]; confirm units in config before continuing')
    conflict=out.loc[out.valid].groupby('sample_id').Age.nunique()
    if (conflict>1).any(): raise ValueError(f'Conflicting ages for samples: {conflict[conflict>1].index.tolist()}')
    return out


def inventory(paths,cfg):
    items=[dict(dataset='train',role='development',method='stored_log1p',expression='X',path=paths['train']),
           dict(dataset='combined',role='external_combined',method='stored_log1p',expression='X',path=paths['external_combined'])]
    for p in sorted(Path(paths['external_gse_dir']).rglob('*.h5ad')):
        method='scanorama' if 'scanorama' in p.name else 'bbknn' if 'bbknn' in p.name else 'other'
        if method != cfg['external']['primary_method'] and not cfg['external']['include_alternate_method']: continue
        items.append(dict(dataset=p.parent.name+'_'+method,role='external_per_gse',method=method,
                          expression=cfg['external'].get(method+'_expression','X'),path=str(p)))
    if len({i['dataset'] for i in items})!=len(items): raise ValueError('Nonunique GSE dataset IDs')
    return items


def matrix_source(a, expression):
    if expression=='raw':
        if a.raw is None: raise ValueError('Configured raw expression is absent; refusing scaled X fallback')
        return a.raw.X, a.raw.var_names
    if expression=='X': return a.X,a.var_names
    if expression in a.layers: return a.layers[expression],a.var_names
    raise ValueError(f'Unknown expression source: {expression}')


def audit(ctx):
    root=ctx.root/'audit'; items=inventory(ctx.paths,ctx.cfg)
    rows=[]; sample_tables={}; obs_tables=[];excluded_samples=[]
    for it in items:
        p=Path(it['path']); LOG.info('AUDIT %s %s',it['dataset'],p)
        if not p.is_file(): raise FileNotFoundError(p)
        a=ad.read_h5ad(p,backed='r')
        try:
            m=metadata(a.obs,ctx.cfg); x,genes=matrix_source(a,it['expression'])
            excluded=m[~m.valid].groupby(['GSE','sample','GSM'],dropna=False).size().reset_index(name='n_cells')
            for r in excluded.itertuples():
                excluded_samples.append(dict(dataset=it['dataset'],GSE=r.GSE,sample=r.sample,GSM=r.GSM,n_cells=r.n_cells,
                    reason='missing/non-numeric age (e.g. Adult) or invalid celltype/sample; no age imputation'))
            if x.shape[0]!=len(m): raise ValueError(f'Expression/obs rows differ for {p}')
            if not genes.is_unique: raise ValueError(f'Duplicate gene IDs in {p}; explicit mapping required')
            probe=x[:min(256,len(m))]; vals=probe.data if sp.issparse(probe) else np.asarray(probe).ravel()
            mn=float(vals.min()) if len(vals) else 0; mx=float(vals.max()) if len(vals) else 0
            if mn< -1e-6 and it['expression']=='raw': raise ValueError(f'Negative raw expression: {p}')
            tab=m[m.valid].groupby('sample_id',sort=True).agg({**{c:'first' for c in m if c not in ['sample_id','valid','Celltype']},'valid':'size'}).rename(columns={'valid':'n_cells'}).reset_index()
            sample_tables[it['dataset']]=tab
            ct=m[m.valid].groupby(['sample_id','Celltype']).size().rename('n_cells').reset_index()
            write_csv(root/(it['dataset']+'_sample_subtype_counts.csv'),ct)
            write_csv(root/(it['dataset']+'_samples.csv'),tab)
            record=dict(**{**it,**file_info(p)},n_cells=a.n_obs,n_genes=x.shape[1],n_samples=len(tab),
                        n_groups=tab.group_id.nunique(),n_GSE=tab.GSE.nunique(),n_celltypes=m.loc[m.valid,'Celltype'].nunique(),
                        age_min=tab.Age.min(),age_max=tab.Age.max(),age_missing=int((~np.isfinite(m.Age)).sum()),
                        expression_type=type(x).__name__,probe_min=mn,probe_max=mx,
                        integer_fraction=float(np.mean(np.isclose(vals,np.round(vals)))) if len(vals) else np.nan,
                        excluded_cells=int((~m.valid).sum()),layers=';'.join(a.layers),has_raw=a.raw is not None,
                        log1p_uns='log1p' in a.uns,duplicate_genes=int(genes.duplicated().sum()))
            rows.append(record)
            age_key=next(k for k in [ctx.cfg['fields']['age']]+ctx.cfg['fields'].get('age_fallbacks',[]) if k in a.obs)
            obs_tables.append(dict(dataset=it['dataset'],obs_columns=list(a.obs),var_columns=list(a.var),
                                   age_key=age_key,age_dtype=str(a.obs[age_key].dtype),
                                   selected_fields=ctx.cfg['fields'],age_examples=a.obs[age_key].drop_duplicates().head(25).tolist()))
        finally: a.file.close()
    train=sample_tables['train']; combined=sample_tables['combined']; leaks=[]
    for ds,tab in sample_tables.items():
        if ds=='train': continue
        for r in tab.itertuples():
            gsm_match=bool(r.GSM not in NULLS and r.GSM in set(train.GSM))
            sample_match=r.sample_id in set(train.sample_id)
            donor_match=ctx.cfg['fields'].get('donor') and r.group_id in set(train.group_id)
            leaks.append(dict(dataset=ds,sample_id=r.sample_id,group_id=r.group_id,GSM=r.GSM,
                              overlaps_training=bool(gsm_match or sample_match or donor_match),
                              overlaps_combined=r.sample_id in set(combined.sample_id) or (r.GSM not in NULLS and r.GSM in set(combined.GSM)),
                              action='exclude_external' if gsm_match or sample_match or donor_match else 'keep',
                              reason='training_sample/GSM/donor_overlap' if gsm_match or sample_match or donor_match else ''))
    write_csv(root/'input_summary.csv',pd.DataFrame(rows));write_csv(root/'leakage_audit.csv',pd.DataFrame(leaks))
    write_csv(root/'excluded_samples.csv',pd.DataFrame(excluded_samples,columns=['dataset','GSE','sample','GSM','n_cells','reason']))
    write_json(root/'obs_schema.json',obs_tables);write_json(root/'inventory.json',items)
    write_csv(ctx.root/'external/gse_manifest.csv',pd.DataFrame(rows))
    resource=dict(cpu_count=psutil.cpu_count(),affinity=len(psutil.Process().cpu_affinity()),
                  ram_available_gb=psutil.virtual_memory().available/2**30,disk_free_gb=shutil.disk_usage(ctx.root).free/2**30,
                  threads=ctx.cfg['resources']['threads'],benchmark_combinations=len(ctx.cfg['features']['K'])*ctx.cfg['cv']['n_splits']*len(ctx.cfg['models']['names']),
                  adaptation_upper_bound_per_dataset=len(ctx.cfg['features']['K'])*ctx.cfg['cv']['n_splits']*len(ctx.cfg['models']['names'])*len(ctx.cfg['external']['fractions'])*ctx.cfg['external']['repeats'])
    write_json(root/'resource_estimate.json',resource)
    refs=[]
    for key,path in ctx.cfg['reference_paths'].items():
        try:
            p=Path(path); text=p.read_text()
            if p.suffix=='.ipynb': text='\n'.join(''.join(c['source']) for c in json.loads(text)['cells'] if c['cell_type']=='code')
            write_text(root/('reference_'+key+'.source.txt'),text);refs.append(f'{key}: read {path} ({len(text)} source characters; mixed-language notebook cells are saved as text)')
        except (OSError,ValueError) as e: refs.append(f'{key}: unavailable: {e}')
    versions=[]
    for pkg in ['anndata','numpy','pandas','scipy','scikit-learn','torch','xgboost','PyYAML','matplotlib','joblib','psutil','nbformat']:
        try: versions.append(pkg+'=='+importlib.metadata.version(pkg))
        except importlib.metadata.PackageNotFoundError: versions.append(pkg+' MISSING')
    write_text(root/'environment.txt','\n'.join([sys.executable,sys.version,platform.platform()]+versions)+'\n')
    import yaml
    write_text(root/'config.yaml',yaml.safe_dump(ctx.cfg,sort_keys=False))
    report=f'''# Input audit: {ctx.lineage}

Real inputs: {len(items)} files; development sample rows: {len(train)}; conservative independent groups: {train.group_id.nunique()}.
Fields: {ctx.cfg['fields']}. No verified donor field: {ctx.cfg['fields'].get('donor') is None}.
Rows are biological sample IDs GSE|sample, stripping only the known trailing _query marker.
Without verified donor IDs we conservatively group the same GSE/GSM, including _1/_2 samples, to prevent cross-fold leakage; independent donor identity is unverified.
Age is the stored numeric Age, interpreted as years (range {train.Age.min()}–{train.Age.max()}); fractional neonatal ages are retained. Upstream source unit conversion is not independently verified.
External exclusions: {sum(r['overlaps_training'] for r in leaks)} sample representations. A and B overlap: {sum(r['overlaps_combined'] for r in leaks if r['dataset']!='combined')} representations; they are not independent replication cohorts.
Celltype retains pan/tissue suffixes exactly. Myfibro enrichment is described, never forcibly balanced at feature selection.

## Audited reference behavior and changes

- Legacy 10-model pseudobulk averages stored log expression per sample/subtype, zero-fills absent subtype, then variance/MI/scaling; keep that mean and absence policy, move all fitted transformations inside training folds.
- Legacy regression_metrics MAD = mean(abs(y_true - mean(y_true))); preserve it as MAD/Age_MAD. fig3 instead labels median prediction error MAD: here that error is explicitly MedAE.
- Legacy estimator references can alias across folds; create a fresh model and save independent artifacts per combination.
- Legacy IDEAL-Age pools zero-expanded cells over all cells, indirectly using subtype proportions, builds a train/external feature union, and fits stacking on OOF targets before reporting those same targets. These behaviors conflict with this task.
- New IDEAL-Age variant uses the same K sample/subtype schema, deterministic per-subtype cell reservoirs and equal subtype pooling before gated attention. Pseudobulk averages four transformed regressors; Ensemble is a preregistered 50:50 blend, with no outer-validation-trained stacking. These are composition-free comparable variants, not exact legacy IDEAL-Age performance replications. Subtype presence still reflects observability.
- bbknn per-GSE X is scaled/HVG-restricted; use its nonnegative raw log expression. Scanorama X is a separately labeled alternate expression representation; never pool it with bbknn samples.
- EN_finetune refits on external subsets and evaluates the full set, including fitted samples, with optional error/negative-prediction filtering. Replace with disjoint adaptation/holdout and preserve all predictions. Linear/deep models use refit or fine_tune; trees/nonincremental predictors use explicitly labeled slope/intercept calibration.
- Stored training preprocessing upstream of these h5ads (gene filtering, integration, transferred subtype labels) cannot be undone. Fold-level guarantees start from the supplied data; no claim that upstream integration was fold-blind.

## References actually inspected

'''+'\n'.join(refs)+'\n\nCPU/RAM/disk estimates are in resource_estimate.json. Wall time depends on measured model fits; no fabricated completion estimate. Sequential lineage execution is default.\n'
    write_text(root/'audit_report.md',report)
    return [root/'input_summary.csv',root/'leakage_audit.csv',root/'excluded_samples.csv',root/'inventory.json',root/'audit_report.md',root/'environment.txt']


def sparse_save(path,x):
    with atomic_path(path) as tmp:
        with tmp.open('wb') as f: sp.save_npz(f,x,compressed=True)


def aggregate(ctx,it,schema=None):
    """One sparse chunk at a time; no whole-cell dense expression expansion."""
    out=ctx.root/'features'/it['dataset'];out.mkdir(parents=True,exist_ok=True)
    marker=out/'complete.json'
    if ctx.cache.valid(marker): return load_dataset(out)
    a=ad.read_h5ad(it['path'],backed='r')
    try:
        m=metadata(a.obs,ctx.cfg); source,gn=matrix_source(a,it['expression']); genes=list(map(str,gn))
        if schema is None:
            ctypes=sorted(m.loc[m.valid,'Celltype'].unique());schema={'celltypes':ctypes,'genes':genes,'separator':'_'}
        ctypes=schema['celltypes']; train_genes=schema['genes'];ng=len(train_genes);nc=len(ctypes)
        gi=pd.Index(train_genes).get_indexer(genes);ci=pd.Index(ctypes).get_indexer(m.Celltype)
        keep=m.valid.to_numpy() & (ci>=0)
        if it['dataset']!='train':
            leaks=pd.read_csv(ctx.root/'audit/leakage_audit.csv'); banned=set(leaks.loc[(leaks.dataset==it['dataset'])&(leaks.action=='exclude_external'),'sample_id'])
            keep &= ~m.sample_id.isin(banned).to_numpy()
        counts=m.loc[keep].groupby('sample_id').size();ids=sorted(counts[counts>=ctx.cfg['features']['min_cells_per_sample']].index)
        si=pd.Index(ids).get_indexer(m.sample_id);keep &= si>=0
        tab=m.loc[keep].groupby('sample_id').first().reindex(ids).drop(columns=['Celltype','valid']).reset_index()
        sums=sp.csr_matrix((len(ids)*nc,ng),dtype=np.float32);ctcounts=np.zeros(len(ids)*nc,np.int64)
        rng=np.random.default_rng(ctx.cfg['seed']);reservoir=[]
        groups={}
        for ix in np.flatnonzero(keep): groups.setdefault(int(si[ix]*nc+ci[ix]),[]).append(ix)
        for group,index in groups.items():
            if len(index)<ctx.cfg['features']['min_cells_per_subtype']: keep[index]=False;continue
            reservoir.extend(rng.choice(index,min(len(index),ctx.cfg['features']['cells_per_subtype']),replace=False).tolist())
        reservoir=np.array(sorted(reservoir),int);res_index={int(row):i for i,row in enumerate(reservoir)}
        bagchunks=[];sample_rows=[];ct_rows=[]
        chunk=ctx.cfg['features']['chunk_cells'];validgenes=gi>=0
        for start in range(0,len(m),chunk):
            stop=min(start+chunk,len(m)); local=np.flatnonzero(keep[start:stop])
            if not len(local):continue
            x=sp.csr_matrix(source[start:stop],dtype=np.float32)[local]
            if not np.isfinite(x.data).all(): raise ValueError('Nonfinite input expression')
            # Map available genes only; training schema never changes with external inputs.
            co=x[:,validgenes].tocoo();mapped=sp.csr_matrix((co.data,(co.row,gi[validgenes][co.col])),shape=(len(local),ng))
            globalrows=local+start;gr=si[globalrows]*nc+ci[globalrows]
            membership=sp.csr_matrix((np.ones(len(gr)),(gr,np.arange(len(gr)))),shape=(len(ids)*nc,len(gr)))
            sums += (membership@mapped).astype(np.float32);np.add.at(ctcounts,gr,1)
            selected=np.array([j for j,row in enumerate(globalrows) if int(row) in res_index],int)
            if len(selected):
                co=mapped[selected].tocoo();chosen=globalrows[selected]
                expanded=sp.csr_matrix((co.data,(co.row,ci[chosen][co.row]*ng+co.col)),shape=(len(chosen),nc*ng))
                bagchunks.append(expanded);sample_rows.extend(si[chosen]);ct_rows.extend(ci[chosen])
            if start%(chunk*100)==0:LOG.info('%s aggregation %d/%d cells',it['dataset'],stop,len(m))
        means=sp.diags(1/np.maximum(ctcounts,1))@sums;co=means.tocoo()
        X=sp.csr_matrix((co.data,(co.row//nc,(co.row%nc)*ng+co.col)),shape=(len(ids),nc*ng),dtype=np.float32)
        bags=sp.vstack(bagchunks,format='csr') if bagchunks else sp.csr_matrix((0,nc*ng),dtype=np.float32)
        feature_names=[f'{ct}_{g}' for ct in ctypes for g in train_genes]
        if len(set(feature_names))!=len(feature_names):raise ValueError('Ambiguous Celltype_Gene names')
        sparse_save(out/'pseudobulk.npz',X);sparse_save(out/'cell_reservoir.npz',bags)
        write_csv(out/'samples.csv',tab);write_json(out/'schema.json',{**schema,'features':feature_names})
        dump(out/'bag_index.joblib',dict(sample=np.asarray(sample_rows,int),celltype=np.asarray(ct_rows,int)))
        write_csv(out/'celltype_coverage.csv',pd.DataFrame(dict(Celltype=ctypes,n_cells=ctcounts.reshape(len(ids),nc).sum(0),
                     n_samples=(ctcounts.reshape(len(ids),nc)>0).sum(0))))
        write_json(out/'alignment.json',dict(**it,n_samples=len(ids),n_features=X.shape[1],matched_genes=int(validgenes.sum()),
                    missing_genes=sorted(set(train_genes)-set(genes)),unknown_celltypes=sorted(set(m.Celltype)-set(ctypes)),
                    missing_policy='zero for unobserved gene/subtype; source expression retained',reservoir_cells=bags.shape[0]))
        outputs=[out/f for f in ['pseudobulk.npz','cell_reservoir.npz','samples.csv','schema.json','bag_index.joblib','alignment.json','celltype_coverage.csv']]
        ctx.cache.mark(marker,outputs)
    finally:a.file.close()
    return load_dataset(out)


def load_dataset(path):
    import joblib
    p=Path(path)
    return dict(path=p,X=sp.load_npz(p/'pseudobulk.npz'),cells=sp.load_npz(p/'cell_reservoir.npz'),
                meta=pd.read_csv(p/'samples.csv',keep_default_na=False),schema=json.loads((p/'schema.json').read_text()),
                bag_index=joblib.load(p/'bag_index.joblib'))


def features(ctx):
    items=json.loads((ctx.root/'audit/inventory.json').read_text())
    train=aggregate(ctx,items[0]);schema=train['schema']; del train
    for it in items[1:]:
        dataset=aggregate(ctx,it,schema);del dataset
    return [ctx.root/'features'/it['dataset']/'complete.json' for it in items]
