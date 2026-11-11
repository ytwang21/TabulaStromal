"""Fixed, group-disjoint 30% adaptation integrated with existing external tasks.

Import the correct lineage's tabulastromal package before importing this module.
No development estimator, training scaler or feature selection is refitted here.
"""
from pathlib import Path
from functools import lru_cache
import copy,itertools,json,logging,os,time,traceback
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import ElasticNet,Lasso,Ridge
from sklearn.model_selection import train_test_split
from tabulastromal import models,training,pipeline
from tabulastromal.storage import Cache,digest,file_hash,write_json,write_csv,write_text,dump
from tabulastromal.metrics import evaluate,METRICS

ROOT=Path(__file__).resolve().parent
LOG=logging.getLogger(__name__)
OPTIONS={}
STRATEGIES={name:'refit_train_plus_adapt' for name in models.NAMES}
STRATEGIES.update(ElasticNet='refit_adapt_only_EN_reference',Lasso='refit_adapt_only',
                  Supervised_PCA_Ridge='refit_train_plus_adapt_frozen_pca')
STRATEGIES.update(XGBoost='continue_boosting',GradientBoosting='continue_boosting',
                  MLP='partial_fit',DeepSets='fine_tune_head',Pseudobulk='continue_bulk_boosting',
                  Ensemble='fine_tune_head_continue_bulk_boosting')
_CONTEXT={}

def selected_folds(ctx):
    """One source fold per model/K: highest internal Valid R2, then MAE, fold.

    External labels, errors and outcomes never enter this selection.
    """
    cv=pd.read_csv(ctx.root/'benchmark/cv_fold_metrics.csv')
    valid=cv[cv.split=='Valid'].copy()
    rows=valid.sort_values(['R2','MAE','fold'],ascending=[False,True,True],na_position='last').groupby(['K','model'],sort=True).head(1)
    return [(int(r.K),r.model,int(r.fold)) for r in rows.itertuples()]

def continue_bulk(bulk,X,y,Xcheck):
    """Continue the saved GBM branch; freeze source quantile and linear branches."""
    import pickle
    frozen=pickle.dumps((bulk.qt,bulk.models[:3]),protocol=5)
    booster=bulk.models[3]
    original_count=len(booster.estimators_)
    before=booster.predict(bulk.qt.transform(Xcheck))
    rounds=min(OPTIONS['continuation_boost_rounds'],int(booster.n_estimators))
    booster.set_params(warm_start=True,n_estimators=int(booster.n_estimators)+rounds)
    booster.fit(bulk.qt.transform(X),y)
    stages=list(booster.staged_predict(bulk.qt.transform(Xcheck)))
    np.testing.assert_allclose(stages[original_count-1],before,rtol=1e-9,atol=1e-9)
    assert pickle.dumps((bulk.qt,bulk.models[:3]),protocol=5)==frozen
    return f'original quantile/EN/Ridge/Lasso frozen; saved bulk GBM continued with {rounds} residual stages'

def configure(cfg,lineage,options):
    global OPTIONS
    OPTIONS=dict(options,lineage=lineage)
    import endo_acceleration as fast
    runtime=dict(workers=options['workers_per_lineage'],bootstrap_workers=options['bootstrap_workers'],
        worker_threads=1,io_workers=8,reserve_memory_gb=options['reserve_memory_gb'],
        heartbeat_seconds=options['heartbeat_seconds'],solver='sklearn',linear_float64=True,
        linear_selection='random',linear_random_state=42,gpu_enabled=options['gpu_enabled'],
        gpu_ids=options['gpu_ids'],gpu_min_free_gb=16,gpu_memory_fraction=.10,
        gpu_adaptation=True,gpu_min_samples=options['gpu_min_samples'],
        scratch_dir=str(ROOT/'runtime/endo_acceleration'))
    fast.context=lambda c:pipeline.Context(c,lineage,True)
    fast.install(cfg,runtime)
    # Both lineages use the original sklearn classes; double precision/random
    # traversal retains objective and strict convergence threshold.
    return runtime


def worker_init(cfg,lineage,options):
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[k]='1'
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,options['gpu_ids'])) if options['gpu_enabled'] else ''
    from threadpoolctl import threadpool_limits
    configure(cfg,lineage,options);threadpool_limits(1)
    try:os.nice(5)
    except OSError:pass


def context(cfg):
    key=digest(cfg)
    if key not in _CONTEXT:_CONTEXT[key]=pipeline.Context(cfg,OPTIONS['lineage'],True)
    return _CONTEXT[key]


@lru_cache(maxsize=2)
def arrays(path,cells=False):
    p=Path(path)
    if cells:
        from tabulastromal.data import load_dataset
        return load_dataset(p)
    return dict(path=p,X=sp.load_npz(p/'pseudobulk.npz'),
        meta=pd.read_csv(p/'samples.csv',keep_default_na=False),
        schema=json.loads((p/'schema.json').read_text()))


def valid(ctx,path,protocol):
    p=Path(path)
    return p.exists() and json.loads(p.read_text()).get('ft30_protocol')==protocol and ctx.cache.valid(p)


def prepare(ctx):
    """Audit metadata, pin CV bytes and define one split for every model/K/fold."""
    train=arrays(str(ctx.root/'features/train'))
    ext=arrays(str(ctx.root/'features/combined'))
    meta=ext['meta'].copy()
    # Existing cache groups use GSE|GSM/sample because true donor IDs are absent.
    # Fail instead of silently ignoring a newly provided explicit donor field.
    import anndata as ad
    observed={}
    for role,key in [('train','train'),('combined','external_combined')]:
        a=ad.read_h5ad(ctx.paths[key],backed='r')
        try:
            columns=list(a.obs.columns)
            donors=[c for c in ['donor_id','donor','Donor','individual_id','patient_id','subject_id'] if c in columns]
            observed[role]=dict(path=ctx.paths[key],columns=columns,donor_columns=donors)
            if donors and not ctx.cfg['fields'].get('donor'):
                raise ValueError('Explicit donor metadata needs verified mapping before FT30: '+str(donors))
        finally:a.file.close()
    overlap=meta.sample_id.isin(train['meta'].sample_id)|meta.group_id.isin(train['meta'].group_id)
    gsm=train['meta'].GSM.astype(str)
    known=set(gsm[gsm.str.match(r'^GSM\d+$')])
    overlap |= meta.GSM.astype(str).isin(known)
    excluded=meta[overlap].copy();meta=meta[~overlap].copy()
    if not np.isfinite(meta.Age.to_numpy(float)).all():raise ValueError('Non-numeric external age')
    groups=meta.groupby('group_id',sort=True).agg(GSE=('GSE','first'),Age=('Age','median'))
    n_adapt=int(round(len(groups)*OPTIONS['adaptation_fraction']))
    if n_adapt<2 or len(groups)-n_adapt<2:raise ValueError('Too few independent groups for 30/70')
    counts=groups.GSE.value_counts();strata=groups.GSE
    stratifiable=counts.min()>=2 and min(n_adapt,len(groups)-n_adapt)>=len(counts)
    reason='GSE stratification feasible' if stratifiable else 'GSE strata too rare/numerous; fixed unstratified group split; no age/performance-based split search'
    ai,hi=train_test_split(groups.index.to_numpy(),train_size=n_adapt,random_state=OPTIONS['adaptation_seed'],
                          stratify=strata if stratifiable else None)
    assert not set(ai)&set(hi)
    meta['set']=np.where(meta.group_id.isin(ai),'adaptation','holdout')
    meta['donor_id']=meta.group_id;meta['age_available']=np.isfinite(meta.Age)
    meta['independence_unit']='GSE|GSM/sample proxy; true donor identity unavailable' if not ctx.cfg['fields'].get('donor') else 'configured donor'
    protocol=digest(dict(protocol=OPTIONS['protocol'],source_signature=ctx.cache.signature,
        seed=OPTIONS['adaptation_seed'],adaptation_ids=sorted(meta.loc[meta.set=='adaptation','sample_id']),
        holdout_ids=sorted(meta.loc[meta.set=='holdout','sample_id']),strategies=STRATEGIES,
        code=file_hash(__file__),reference_notebook=file_hash(ctx.cfg['reference_paths']['fine_tune']),evaluation='all eligible combined samples, including adaptation samples',
        selection=selected_folds(ctx),continuation_epochs=OPTIONS['continuation_epochs'],
        continuation_boost_rounds=OPTIONS['continuation_boost_rounds']))
    ref=ctx.root/'reports/ft30_cv_reference.json'
    cvpaths=['benchmark/cv_fold_metrics.csv','benchmark/cv_summary.csv','benchmark/cv_predictions.csv','benchmark/model_comparison.csv']
    hashes={p:file_hash(ctx.root/p) for p in cvpaths}
    if ref.exists() and json.loads(ref.read_text())!=hashes:raise RuntimeError('Development CV changed since FT30 audit')
    if not ref.exists():write_json(ref,hashes)
    write_csv(ctx.root/'external/external_ft30_split.csv',meta)
    write_csv(ctx.root/'external/external_sample_split_30_70.csv',meta)
    write_csv(ctx.root/'external/adaptation_30_sample_ids.csv',meta[meta.set=='adaptation'])
    write_csv(ctx.root/'external/heldout_70_sample_ids.csv',meta[meta.set=='holdout'])
    qc=[]
    for field in ['GSE','tissue']:
        if field in meta:
            q=meta.groupby(['set',field],dropna=False).size().rename('n_samples').reset_index().rename(columns={field:'stratum'})
            q['field']=field;qc.append(q)
    write_csv(ctx.root/'external/external_split_strata_QC.csv',pd.concat(qc,ignore_index=True))
    write_csv(ctx.root/'external/external_ft30_excluded_overlap.csv',excluded)
    summary=meta.groupby('set').agg(n_samples=('sample_id','size'),n_groups=('group_id','nunique'),
        age_min=('Age','min'),age_max=('Age','max'),age_mean=('Age','mean')).reset_index()
    summary['seed']=OPTIONS['adaptation_seed'];summary['stratification']=reason
    write_csv(ctx.root/'external/external_ft30_split_summary.csv',summary)
    write_json(ctx.root/'external/ft30_protocol.json',dict(ft30_protocol=protocol,options=OPTIONS,
        observed_metadata=observed,strategies=STRATEGIES,split_rule=reason,
        n_adapt_groups=len(ai),n_holdout_groups=len(hi),n_evaluation_samples=len(meta),
        evaluation_scope='all100_includes_adaptation',n_excluded_overlap=len(excluded),
        benchmark_hashes=hashes,per_gse_policy=OPTIONS['per_gse_policy']))
    manifest=[];selected=set(selected_folds(ctx))
    for K,name,fold in itertools.product(ctx.cfg['features']['K'],ctx.cfg['models']['names'],range(1,ctx.cfg['cv']['n_splits']+1)):
        p=training.task_path(ctx,K,name,fold)
        manifest.append(dict(lineage=ctx.lineage,K=K,model=name,fold=fold,source_model_path=str(p),
            preprocessor_path=str(p/'preprocessor.joblib'),source_marker_exists=(p/'complete.json').exists(),
            adaptation_strategy=STRATEGIES[name],adaptation_fraction=.3,adaptation_seed=42,
            selected_for_external=(K,name,fold) in selected,
            selection_rule='max internal Valid R2; min Valid MAE; min fold; no external selection'))
    write_csv(ctx.root/'ft30_task_manifest.csv',pd.DataFrame(manifest))
    chosen=pd.DataFrame(manifest);chosen=chosen[chosen.selected_for_external].copy()
    cv=pd.read_csv(ctx.root/'benchmark/cv_fold_metrics.csv');v=cv[cv.split=='Valid'][['model','K','fold','R2','MAE']]
    chosen=chosen.merge(v,on=['model','K','fold']).rename(columns={'R2':'internal_Valid_R2','MAE':'internal_Valid_MAE'})
    chosen['fold_best']=chosen.fold
    write_csv(ctx.root/'external/fold_selection_by_internal_valid.csv',chosen)
    import yaml
    write_text(ctx.root/'ft30_config.yaml',yaml.safe_dump(OPTIONS,sort_keys=False))
    notes=f'''# Combined external all100 audit: {ctx.lineage}

Actual result root: {ctx.root}. All {len(manifest)} source development models are preserved.
Only {len(selected)} internally selected best folds (one per model/K) are evaluated:
highest internal Valid R2, then lowest Valid MAE, then lowest fold number.
No external-performance-based fold/K/model selection.

Fixed seed-42 adaptation pool: {len(ai)} groups; remaining pool: {len(hi)} groups.
ALL {len(meta)} eligible combined samples are predicted/scored, INCLUDING samples used
for supervised adaptation. These are mixed apparent/unseen metrics, not independent
held-out performance. The same frozen split is reused by all trainable models.
Explicit sample/group/GSM training overlap is excluded: {len(excluded)} rows.
No verified true donor column: GSE|GSM/sample is a proxy; cross-GSE donor identity is unknown.

Actual strategies: {json.dumps(STRATEGIES)}.
The new AllSample instruction supersedes the previous strict-continuation-only policy.
EN reference {ctx.cfg['reference_paths']['fine_tune']} cell5 uses deepcopy(base).fit(X_adapt,y_adapt), then predict(X_all).
Here ElasticNet is explicitly refit_adapt_only_EN_reference, with fixed source alpha/l1_ratio/scaler, NOT falsely called incremental fine-tuning.
The source notebook's fraction grid, favorable outcome search and posthoc error/negative-prediction QC are NOT used.
Lasso uses a cold-start adaptation-only refit with fixed source alpha after source-coefficient warm-start stalled; the L1 objective is unchanged. Ridge/SVR/forests refit source-fold TRAIN+adaptation using fixed source hyperparameters.
Supervised_PCA_Ridge freezes the saved PCA and refits its Ridge on source TRAIN+adaptation.
XGB/GB continue boosting; MLP partial_fit; DeepSets continues source attention/head;
Pseudobulk/Ensemble continue saved GBM branches with frozen quantile/linear transformations.
All preprocessing remains fitted only on the original source TRAIN IDs. No external heldout labels enter fitting or strategy selection.
Failure falls back ONLY to the same original best fold, with traceback and fallback_no_adaptation status.
Full-all100 metrics feed the primary figures; source=external_all_post_adaptation.
The same never-adapted70% is separately scored before/after, source=external_heldout_70, never substituted into the primary figure.

Original MAD = mean(abs(Age - mean(Age))); it depends only on the evaluation ages.
The original benchmark code Fibro_GeneModel_10_bench.py line 709 confirms this.
Identical MAD across models on the same folds is expected, not reused model predictions.
MedAE = median(abs(Age - prediction)) is separately plotted as a model error metric.
Original four CV CSV SHA256 values remain pinned and are verified unchanged.

Only combined external tasks/figures are updated. Per-GSE predictions and biology
outputs are retained. Existing unselected fold task files are historical and are NOT
consumed; ft30_task_manifest.csv marks selected_for_external explicitly.
Prior outputs are retained with .pre_reference_all suffix in the same directories. User explicitly authorized overwriting current figure/CSV names; raw h5ad, development bundles and internal CV are preserved.
'''
    write_text(ctx.root/'ft30_audit_report.md',notes)
    return protocol


def source_and_target(ctx,item):
    K,name,fold=item
    if fold=='final':return ctx.root/'models/final_model',ctx.root/'external', 'final_'
    source=training.task_path(ctx,int(K),name,int(fold))
    target=ctx.root/'external/tasks/combined'/f'K{K}'/name/f'fold{fold}'/'adapt_f0.3_r0'
    if not target.is_dir():raise FileNotFoundError('Expected existing task directory: '+str(target))
    return source,target,''


def fit_checked_adapter(estimator,X,y):
    """Respect saved celer interfaces; sklearn-only solver flags are not portable."""
    params=estimator.get_params() if hasattr(estimator,'get_params') else {}
    if isinstance(estimator,(ElasticNet,Lasso)) and 'selection' not in params:
        import endo_acceleration as fast
        X=np.asarray(X,dtype=np.float64,order='F');y=np.asarray(y,dtype=np.float64)
        return fast._ORIGINAL['checked'](estimator,X,y)
    return models.fit_checked(estimator,X,y)


def fit_adapter(ctx,source,target,prefix,protocol,name,base,pre,info,ext,ai,ei):
    """Use fixed source preprocessing/HP; labels from adaptation or source TRAIN only."""
    from sklearn.base import clone
    import endo_acceleration as fast
    cfg=ctx.cfg;strategy=STRATEGIES[name];adapter=copy.deepcopy(base)
    Xa=pre.transform(ext['X'][ai]);y=ext['meta'].Age.to_numpy(float)
    Xcheck=pre.transform(ext['X'][ei]);source_refit_ids=[];paths=[]
    uses_weights=True
    if name=='ElasticNet':
        # EN_finetune.ipynb cell 5: deepcopy(base).fit(X_ft,y_ft), then predict(X_full).
        # Explicitly start a fresh solver so inherited retry warm_start flags do not
        # silently change that reference method. Alpha/l1_ratio/scaler remain fixed.
        adapter.set_params(warm_start=False)
        adapter=fit_checked_adapter(adapter,Xa,y[ai]);uses_weights=False
        detail='EN reference adaptation-only refit: deepcopy source, fit 30%, predict 100%; fixed source alpha/l1_ratio and scaler; not incremental fine-tuning'
    elif name=='Lasso':
        adapter.set_params(warm_start=False)
        adapter=fit_checked_adapter(adapter,Xa,y[ai]);uses_weights=False
        detail='Adaptation-only Lasso refit with fixed source alpha/scaler; cold solver start avoids stalled source-coefficient warm start; not incremental fine-tuning'
    elif name in ['Ridge','SVR','ExtraTrees','RandomForest','Supervised_PCA_Ridge']:
        train=arrays(str(ctx.root/'features/train'))
        idx=pd.Series(np.arange(len(train['meta'])),index=train['meta'].sample_id)
        tr=idx.loc[info['train_ids']].to_numpy(int)
        Xt=np.vstack([pre.transform(train['X'][tr]),Xa])
        yt=np.r_[train['meta'].Age.to_numpy(float)[tr],y[ai]]
        source_refit_ids=list(info['train_ids']);uses_weights=False
        if name=='Supervised_PCA_Ridge':
            # The source PCA is frozen. Refit only its fixed-alpha Ridge branch.
            adapter.model=models.fit_checked(clone(base.model),base.pca.transform(Xt),yt)
            detail='Refit source fixed-alpha Ridge on source-fold TRAIN + adaptation; source PCA and outer preprocessing frozen'
        else:
            adapter=clone(base)
            if 'n_jobs' in adapter.get_params():adapter.set_params(n_jobs=1)
            adapter=models.fit_checked(adapter,Xt,yt)
            detail='Fixed-hyperparameter refit on source-fold TRAIN + 30% adaptation; source Valid and external heldout excluded; source preprocessing frozen'
    elif name=='XGBoost':
        rounds=min(OPTIONS['continuation_boost_rounds'],int(base.n_estimators))
        old_rounds=base.get_booster().num_boosted_rounds()
        adapter.set_params(n_estimators=rounds,n_jobs=1)
        adapter.fit(Xa,y[ai],xgb_model=base.get_booster())
        assert adapter.get_booster().num_boosted_rounds()==old_rounds+rounds
        assert adapter.get_booster().get_dump()[:old_rounds]==base.get_booster().get_dump()
        detail=f'Continue saved booster with {rounds} additional rounds; original trees unchanged'
    elif name=='GradientBoosting':
        rounds=min(OPTIONS['continuation_boost_rounds'],int(base.n_estimators))
        adapter.set_params(warm_start=True,n_estimators=int(base.n_estimators)+rounds)
        adapter.fit(Xa,y[ai])
        stages=list(adapter.staged_predict(Xcheck))
        np.testing.assert_allclose(stages[len(base.estimators_)-1],base.predict(Xcheck),rtol=1e-9,atol=1e-9)
        detail=f'Continue saved GradientBoosting with {rounds} residual stages; original trees/initializer preserved'
    elif name=='MLP':
        adapter.set_params(early_stopping=False)
        for epoch in range(OPTIONS['continuation_epochs']):adapter.partial_fit(Xa,y[ai])
        detail=f'Continue source weights/optimizer with {OPTIONS["continuation_epochs"]} partial_fit epochs'
    elif name=='Pseudobulk':detail=continue_bulk(adapter,Xa,y[ai],Xcheck)
    elif name in ['DeepSets','Ensemble']:
        deep=adapter if name=='DeepSets' else adapter.deep
        deep.net.weight.requires_grad_(False);deep.net.bias.requires_grad_(False)
        for param in deep.net.encoder.parameters():param.requires_grad_(False)
        deep.fit(fast.bags(ext,pre,ai),y[ai],target/(prefix+'ft30_epoch_'+protocol[:12]+'.pt'),True,True)
        path=target/(prefix+'adapter_ft30_state_dict.pt');deep.save(path);paths.append(path)
        detail='Continue source attention/head; source cell encoder frozen'
        if name=='Ensemble':
            detail+='; '+continue_bulk(adapter.bulk,Xa,y[ai],Xcheck)
            path=target/(prefix+'adapter_ft30_bulk.joblib');dump(path,adapter.bulk);paths.append(path)
    else:raise NotImplementedError('No verified adaptation API for '+name)
    bags=fast.bags(ext,pre,ei) if name in ['DeepSets','Ensemble'] else None
    after=models.predict_model(adapter,name,Xcheck,bags)
    if not np.isfinite(after).all():raise ValueError('Nonfinite adapted predictions')
    if name not in ['DeepSets','Ensemble']:
        path=target/(prefix+'adapter_ft30.joblib');dump(path,adapter);paths.append(path)
    return np.asarray(after,dtype=np.float64),detail,source_refit_ids,uses_weights,paths


def task(cfg,item,protocol):
    ctx=context(cfg);K,name,fold=item
    source,target,prefix=source_and_target(ctx,item);marker=target/(prefix+'ft30.complete.json')
    if valid(ctx,marker,protocol):return
    import endo_acceleration as fast
    fast.job_log(ctx,'ft30',f'{name}_K{K}_fold{fold}')
    if not ctx.cache.valid(source/'complete.json'):raise RuntimeError('Invalid source checkpoint '+str(source))
    split=pd.read_csv(ctx.root/'external/external_ft30_split.csv',keep_default_na=False)
    ext=arrays(str(ctx.root/'features/combined'),name in ['DeepSets','Ensemble']);meta=ext['meta']
    index=pd.Series(np.arange(len(meta)),index=meta.sample_id)
    ai=index.loc[split.loc[split.set=='adaptation','sample_id']].to_numpy(int)
    hi=index.loc[split.loc[split.set=='holdout','sample_id']].to_numpy(int)
    ei=index.loc[split.sample_id].to_numpy(int)
    assert not set(meta.group_id.iloc[ai])&set(meta.group_id.iloc[hi])
    _,base,pre,info=models.load_bundle(source,cfg)
    assert set(pre.fit_ids)==set(info['train_ids'])
    assert set(info['train_ids']).isdisjoint(split.sample_id)
    no=ctx.root/'external/tasks/combined'/f'K{K}'/name/f'fold{fold}'/'no_adapt'
    if ctx.cache.valid(no/'complete.json'):
        p=pd.read_csv(no/'predictions.csv').set_index('sample_id')
        prediction=p.loc[meta.sample_id,'prediction'].to_numpy(float)
        np.testing.assert_allclose(p.loc[meta.sample_id,'Age'],meta.Age,rtol=0,atol=1e-10)
    else:
        bags=fast.bags(ext,pre,np.arange(len(meta))) if name in ['DeepSets','Ensemble'] else None
        prediction=models.predict_model(base,name,pre.transform(ext['X']),bags)
    if not np.isfinite(prediction[ei]).all():raise RuntimeError('Original best-fold prediction is nonfinite')
    strategy=STRATEGIES.get(name,'none_fallback');requested=strategy
    paths=[];fit_ids=[];source_refit_ids=[];reason='';error_path=None;uses_weights=False
    try:
        after,detail,source_refit_ids,uses_weights,paths=fit_adapter(ctx,source,target,prefix,protocol,name,base,pre,info,ext,ai,ei)
        external_fit_ids=meta.sample_id.iloc[ai].tolist();fit_ids=source_refit_ids+external_fit_ids
        adaptation_status='adapted'
    except Exception as exc:
        # Failure never triggers score-based method selection. Original prediction
        # is the only fallback and the full error is retained alongside the task.
        reason=f'{type(exc).__name__}: {exc}'
        error_path=target/(prefix+'adaptation_failure_ft30.json')
        write_json(error_path,dict(requested_strategy=requested,error=reason,traceback=traceback.format_exc(),
            fallback='same original internally selected best fold',ft30_protocol=protocol))
        LOG.exception('Adaptation failed; using original best fold %s',item)
        after=np.asarray(prediction[ei],dtype=np.float64);external_fit_ids=[];source_refit_ids=[];fit_ids=[];paths=[]
        strategy='none_fallback';uses_weights=True;adaptation_status='fallback_no_adaptation'
        detail='Original best-fold predictions after adaptation failure; no external-label fit used'
    assert set(fit_ids).isdisjoint(meta.sample_id.iloc[hi])
    assert set(source_refit_ids).isdisjoint(info['valid_ids'])
    common=dict(lineage=ctx.lineage,K=int(K),model=name,fold=int(fold),fold_best=int(fold),dataset='combined',GSE='__combined__',
        method='stored_log1p',source_kind='combined',primary=True,
        adaptation_fraction=.30 if external_fit_ids else 0.,nominal_adaptation_fraction=.30,adaptation_seed=42,
        n_adapt=len(external_fit_ids),n_refit_source=len(source_refit_ids),n_adaptation_pool=len(ai),
        n_holdout=len(hi),n_test=len(ei),n_evaluation=len(ei),evaluation_scope='all100_includes_adaptation',
        adaptation_status=adaptation_status,adaptation_failure_reason=reason,requested_strategy=requested,
        adaptation_strategy=strategy,strategy=strategy,strategy_detail=detail,uses_original_fitted_weights=uses_weights,
        source_alpha=getattr(base,'alpha',np.nan),source_l1_ratio=getattr(base,'l1_ratio',np.nan),
        preserves_source_hyperparameters=True,source_preprocessing_frozen=True,status='complete',failure_reason='',
        source_model_path=str(source),source_fold_selection='max internal Valid R2; min MAE; min fold',
        adapted_model_path=';'.join(map(str,paths)),preprocessor_path=str(source/'preprocessor.joblib'),ft30_protocol=protocol,repeat=0)
    baseline=dict(common,adaptation_status='baseline_no_adaptation',adaptation_strategy='none',strategy='none',n_adapt=0,n_refit_source=0,
        adaptation_fraction=0.,adapted_model_path='',uses_original_fitted_weights=True,evaluation_scope='all100_no_adaptation',
        strategy_detail='Unmodified same internally selected source fold')
    baseline.update(split='no_adapt_all100',source='external_all_no_adaptation',**evaluate(meta.Age.iloc[ei],prediction[ei]))
    adapted=dict(common,split='external_all_post_adaptation',source='external_all_post_adaptation',**evaluate(meta.Age.iloc[ei],after))
    for m in METRICS:adapted['baseline_'+m]=baseline[m];adapted['delta_'+m]=adapted[m]-baseline[m]
    table=meta.iloc[ei].copy();table['prediction']=after;table['baseline_prediction']=prediction[ei]
    table['predicted_age']=after;table['pred_before']=prediction[ei];table['pred_after']=after
    table['adaptation_pool']=table.sample_id.isin(meta.sample_id.iloc[ai]);table['is_adaptation']=table.adaptation_pool
    table['is_heldout']=~table.is_adaptation;table['used_for_adaptation']=table.sample_id.isin(external_fit_ids)
    table['sample_role']=np.where(table.used_for_adaptation,'adaptation_included_in_evaluation','not_used_for_adaptation')
    table['donor_id']=table.group_id;table['donor_verified']=False
    for key,value in common.items():table[key]=value
    table['split']='external_all_post_adaptation';table['source']='external_all_post_adaptation'
    held=table[table.is_heldout].copy()
    held_common=dict(common,n_test=len(hi),n_evaluation=len(hi),primary=False,evaluation_scope='heldout_70_never_adapted',
        source='external_heldout_70',split='external_heldout_70')
    hm=dict(held_common,**evaluate(held.Age,held.prediction))
    hb=dict(held_common,adaptation_status='baseline_no_adaptation',adaptation_strategy='none',strategy='none',n_adapt=0,n_refit_source=0,
        source='external_heldout_70_no_adaptation',split='external_heldout_70_no_adaptation',**evaluate(held.Age,held.baseline_prediction))
    for m in METRICS:hm['baseline_'+m]=hb[m];hm['delta_'+m]=hm[m]-hb[m]
    held['source']='external_heldout_70';held['split']='external_heldout_70'
    task_files={'metrics_ft30.csv':pd.DataFrame([adapted]),'baseline_metrics_ft30.csv':pd.DataFrame([baseline]),
        'predictions_ft30.csv':table,'heldout_metrics_ft30.csv':pd.DataFrame([hm]),
        'heldout_baseline_metrics_ft30.csv':pd.DataFrame([hb]),'heldout_predictions_ft30.csv':held}
    outputs=[]
    for filename,frame in task_files.items():
        path=target/(prefix+filename);write_csv(path,frame);outputs.append(path)
    fit=target/(prefix+'split_ft30.json')
    write_json(fit,dict(adaptation_ids=meta.sample_id.iloc[ai].tolist(),holdout_ids=meta.sample_id.iloc[hi].tolist(),
        evaluation_ids=meta.sample_id.iloc[ei].tolist(),fit_ids=fit_ids,external_fit_ids=external_fit_ids,
        source_refit_ids=source_refit_ids,preprocessor_fit_ids=info['train_ids'],source_train_ids=info['train_ids'],
        source_valid_ids=info['valid_ids'],evaluation_scope=common['evaluation_scope'],ft30_protocol=protocol))
    bundle=target/(prefix+'adapter_bundle_ft30.json')
    write_json(bundle,dict(source_model_path=str(source),preprocessor_path=str(source/'preprocessor.joblib'),
        preprocessor_sha256=file_hash(source/'preprocessor.joblib'),source_marker_sha256=file_hash(source/'complete.json'),
        model=name,K=int(K),fold=int(fold),adaptation_status=adaptation_status,adaptation_strategy=strategy,
        adapter_paths=list(map(str,paths)),ft30_protocol=protocol,load_function='ft30.load_adapter'))
    if error_path:outputs.append(error_path)
    ctx.cache.mark(marker,outputs+[fit,bundle]+paths,ft30_protocol=protocol)
    LOG.info('ALL100 COMPLETE %s K%s %s fold%s n=%d n_fit_external=%d status=%s strategy=%s',ctx.lineage,K,name,fold,len(ei),len(external_fit_ids),adaptation_status,strategy)


def load_adapter(ctx,item,protocol):
    """Reload a saved update with its original matching preprocessing bundle."""
    import joblib
    source,target,prefix=source_and_target(ctx,item)
    if not valid(ctx,target/(prefix+'ft30.complete.json'),protocol):raise RuntimeError('Adapter checkpoint is invalid')
    name,base,pre,info=models.load_bundle(source,ctx.cfg)
    status=pd.read_csv(target/(prefix+'metrics_ft30.csv')).iloc[0]
    if status.adaptation_status=='fallback_no_adaptation':return base,pre
    if name=='DeepSets':return models.DeepRegressor.load(target/(prefix+'adapter_ft30_state_dict.pt'),ctx.cfg),pre
    if name=='Ensemble':
        return models.BlendRegressor(models.DeepRegressor.load(target/(prefix+'adapter_ft30_state_dict.pt'),ctx.cfg),
            joblib.load(target/(prefix+'adapter_ft30_bulk.joblib'))),pre
    return joblib.load(target/(prefix+'adapter_ft30.joblib')),pre

def figure_input_audit(ctx,metrics):
    """Link every main benchmark bar to immutable CV or the all100 primary CSV."""
    rank=pd.read_csv(ctx.root/'benchmark/model_comparison.csv').sort_values('internal_rank')
    pairs=rank.groupby('model',sort=False).head(1)[['model','K']]
    cv=pd.read_csv(ctx.root/'benchmark/cv_fold_metrics.csv').merge(pairs,on=['model','K'])
    ext=metrics.merge(pairs,on=['model','K'])
    oldpath=ctx.root/'external/external_combined_ft30_per_fold.csv.pre_reference_all'
    old=pd.read_csv(oldpath) if oldpath.exists() else pd.DataFrame()
    rows=[]
    for split,table in [('Train',cv[cv.split=='Train']),('Valid',cv[cv.split=='Valid']),('External-valid',ext)]:
        for (name,K),part in table.groupby(['model','K']):
            for metric in METRICS:
                value=float(part[metric].mean())
                if split!='External-valid':
                    before=value;oldcsv=str(ctx.root/'benchmark/cv_fold_metrics.csv');oldsource=split
                    newcsv=oldcsv;source=split
                else:
                    previous=old[(old.model==name)&(old.K==K)] if len(old) else pd.DataFrame()
                    before=float(previous[metric].mean()) if len(previous) else np.nan
                    oldcsv=str(oldpath);oldsource='previous all100 external; EN was no_adaptation'
                    newcsv=str(ctx.root/'external/external_ft30_all_sample_metrics.csv');source='external_all_post_adaptation'
                panel=('A' if ctx.lineage=='Fibro' else 'B') if metric=='R2' else (('C' if ctx.lineage=='Fibro' else 'D') if metric=='MAE' else 'benchmark/supplement')
                rows.append(dict(lineage=ctx.lineage,panel=panel,model=name,K=int(K),split=split,metric=metric,
                    statistic='source-fold mean' if split!='External-valid' else 'single internal-best fold',
                    n_source_models=len(part),old_csv=oldcsv,old_source=oldsource,old_value=before,
                    new_csv=newcsv,new_source=source,new_value=value,
                    Train_Valid_unchanged=(value==before) if split!='External-valid' else np.nan,
                    matches_all100_CSV=True if split=='External-valid' else np.nan))
    path=ctx.root/'external/figure_input_audit.csv';write_csv(path,pd.DataFrame(rows));return path


def collect(ctx,protocol,items,failures):
    rows=[];baseline=[];predictions=[];held_metrics=[];held_baseline=[];held_preds=[];markers=[];status=[]
    failed={tuple(f['item']):f for f in failures}
    for item in items:
        source,target,prefix=source_and_target(ctx,item);marker=target/(prefix+'ft30.complete.json')
        key=dict(lineage=ctx.lineage,K=item[0],model=item[1],fold=item[2])
        if valid(ctx,marker,protocol):
            frame=pd.read_csv(target/(prefix+'metrics_ft30.csv'));rows.append(frame)
            baseline.append(pd.read_csv(target/(prefix+'baseline_metrics_ft30.csv')))
            predictions.append(pd.read_csv(target/(prefix+'predictions_ft30.csv')))
            held_metrics.append(pd.read_csv(target/(prefix+'heldout_metrics_ft30.csv')))
            held_baseline.append(pd.read_csv(target/(prefix+'heldout_baseline_metrics_ft30.csv')))
            held_preds.append(pd.read_csv(target/(prefix+'heldout_predictions_ft30.csv')))
            r=frame.iloc[0]
            status.append(dict(key,status='complete',adaptation_status=r.adaptation_status,
                adaptation_strategy=r.adaptation_strategy,requested_strategy=r.requested_strategy,
                adaptation_failure_reason=r.adaptation_failure_reason,n_adapt=int(r.n_adapt),
                n_refit_source=int(r.n_refit_source),uses_original_fitted_weights=r.uses_original_fitted_weights,
                source_model_path=r.source_model_path,adapted_model_path=r.adapted_model_path,failure_reason=''))
            markers.append(marker)
        else:
            reason=failed.get(tuple(item),{}).get('error','missing/invalid checkpoint')
            status.append(dict(key,status='failed',adaptation_status='failed_source_or_task',failure_reason=reason))
            rows.append(pd.DataFrame([dict(key,source='external_all_post_adaptation',status='failed',failure_reason=reason,**{m:np.nan for m in METRICS})]))
    frame=pd.concat(rows,ignore_index=True)
    concat=lambda seq:pd.concat(seq,ignore_index=True) if seq else pd.DataFrame()
    base,preds,hm,hb,hp=map(concat,[baseline,predictions,held_metrics,held_baseline,held_preds])
    outputs=[]
    tables={'external_combined_ft30_per_fold':frame,'external_ft30_all_sample_metrics':frame,
        'external_combined_no_adapt_all100':base,'external_combined_ft30_predictions':preds,
        'external_ft30_all_sample_predictions':preds,'external_combined_ft30_status':pd.DataFrame(status),
        'adaptation_status':pd.DataFrame(status),'external_heldout_70_paired_metrics':hm,
        'external_heldout_70_baseline_metrics':hb,'external_heldout_70_predictions':hp}
    for name,table in tables.items():
        p=ctx.root/'external'/(name+'.csv');write_csv(p,table);outputs.append(p)
    successful=frame[frame.status=='complete']
    summary=successful.groupby(['lineage','K','model'])[METRICS].agg(['mean','std','count'])
    summary.columns=['_'.join(c) for c in summary.columns]
    p=ctx.root/'external/external_combined_ft30_summary.csv';write_csv(p,summary.reset_index());outputs.append(p)
    keep=['lineage','K','model','fold','source','n_adapt','n_holdout','n_evaluation','evaluation_scope',
          'adaptation_strategy','adaptation_status']+[p+m for m in METRICS for p in ['', 'baseline_','delta_']]
    delta=successful.reindex(columns=keep)
    for name,table in [('external_combined_ft30_paired_delta',delta),('ElasticNet_ft30_paired_comparison',delta[delta.model=='ElasticNet']),
        ('external_ft30_strategy_coverage',pd.DataFrame(status).groupby(['model','status','adaptation_status']).size().rename('tasks').reset_index())]:
        p=ctx.root/'external'/(name+'.csv');write_csv(p,table);outputs.append(p)
    schema=ctx.root/'external/external_ft30_all_sample_metrics_schema.md'
    write_text(schema,'''# Metric and sample sources

external_ft30_all_sample_metrics.csv: source=external_all_post_adaptation, ALL eligible combined external samples, including the adaptation30%. Primary A–F/benchmark/K input. Mixed apparent/unseen descriptive performance, not independent external validation.
external_heldout_70_paired_metrics.csv: source=external_heldout_70, only never-adapted samples, with baseline_* and delta_* for the identical source fold/cohort. This tests adaptation generalization conditional on available sample/GSM groups; verified donor independence is unavailable.
MAD and Age_MAD = mean(abs(Age-mean(Age))), age dispersion, not prediction error; MedAE = median(abs(Age-prediction)). MAE,RMSE are prediction error in years; R2,PCC are dimensionless. Undefined values remain NaN and metric_note records why.
fold_best is selected by maximum internal Valid R2, then lower MAE, then lower fold. K/model use original internal CV ranking. External outcomes never select strategies, hyperparameters, folds or K.
adaptation_status=adapted means the explicitly named update/refit succeeded, NOT necessarily incremental fine-tuning. fallback_no_adaptation means adaptation failed/was unavailable; the original SAME best fold is used, with reason/traceback. n_adapt counts external samples actually fitted; n_refit_source counts source TRAIN rows included in a refit.
refit_adapt_only_EN_reference follows notebook cell5 using fixed source alpha/l1_ratio, preprocessor, and adaptation-only fit; it does not preserve source coefficients. refit_train_plus_adapt uses only source-fold TRAIN plus adaptation (not its Valid). Source PCA is frozen for SPCA-Ridge. Lasso uses adaptation-only refit, not partial_fit; cold start avoids slow source-coefficient warm-start convergence.
is_adaptation marks the global fixed30% pool, even for fallback models; used_for_adaptation records actual external fitting. donor_id is an explicitly unverified GSE/GSM/sample proxy, not verified true donor metadata.
''')
    outputs.append(schema)
    outputs.append(figure_input_audit(ctx,frame))
    for rel in ['external_sample_split_30_70.csv','adaptation_30_sample_ids.csv','heldout_70_sample_ids.csv',
        'external_split_strata_QC.csv','fold_selection_by_internal_valid.csv','external_ft30_split_summary.csv']:
        outputs.append(ctx.root/'external'/rel)
    validate_predictions(ctx,protocol)
    return outputs+markers


def validate_predictions(ctx,protocol):
    metrics=pd.read_csv(ctx.root/'external/external_ft30_all_sample_metrics.csv')
    baseline=pd.read_csv(ctx.root/'external/external_combined_no_adapt_all100.csv')
    preds=pd.read_csv(ctx.root/'external/external_ft30_all_sample_predictions.csv')
    held=pd.read_csv(ctx.root/'external/external_heldout_70_paired_metrics.csv')
    hb=pd.read_csv(ctx.root/'external/external_heldout_70_baseline_metrics.csv')
    hp=pd.read_csv(ctx.root/'external/external_heldout_70_predictions.csv')
    split=pd.read_csv(ctx.root/'external/external_ft30_split.csv')
    holdout=set(split.loc[split.set=='holdout','sample_id']);adapt=set(split.loc[split.set=='adaptation','sample_id'])
    assert not holdout&adapt and metrics.source.eq('external_all_post_adaptation').all()
    assert held.source.eq('external_heldout_70').all()
    assert set(map(tuple,metrics[['K','model','fold']].to_numpy()))==set(selected_folds(ctx))
    assert metrics.status.eq('complete').all()
    def row_for(frame,item):
        K,name,fold=item;r=frame[(frame.K==K)&(frame.model==name)&(frame.fold==fold)];assert len(r)==1;return r.iloc[0]
    def check(frame,item,part,column):
        row=row_for(frame,item);actual=evaluate(part.Age,part[column]);assert int(row.n)==len(part)
        for key in METRICS:np.testing.assert_allclose(actual[key],row[key],rtol=1e-8,atol=1e-8,equal_nan=True)
    for item,part in preds.groupby(['K','model','fold']):
        assert not part.sample_id.duplicated().any() and set(part.sample_id)==holdout|adapt
        r=row_for(metrics,item);used=set(part.loc[part.used_for_adaptation.astype(bool),'sample_id'])
        assert used==(set() if r.adaptation_status=='fallback_no_adaptation' else adapt)
        assert set(part.loc[part.is_adaptation.astype(bool),'sample_id'])==adapt
        if r.adaptation_status=='fallback_no_adaptation':np.testing.assert_array_equal(part.prediction,part.baseline_prediction)
        _,target,prefix=source_and_target(ctx,item);recorded=json.loads((target/'split_ft30.json').read_text())
        assert set(recorded['evaluation_ids'])==holdout|adapt and set(recorded['external_fit_ids'])==used
        assert set(recorded['fit_ids'])==used|set(recorded['source_refit_ids'])
        assert set(recorded['fit_ids']).isdisjoint(holdout)
        assert set(recorded['source_refit_ids']).isdisjoint(recorded['source_valid_ids'])
        assert set(recorded['preprocessor_fit_ids']).isdisjoint(holdout|adapt)
        check(metrics,item,part,'prediction');check(baseline,item,part,'baseline_prediction')
        K,name,fold=item;part70=hp[(hp.K==K)&(hp.model==name)&(hp.fold==fold)]
        assert not part70.sample_id.duplicated().any() and set(part70.sample_id)==holdout
        np.testing.assert_allclose(part70.set_index('sample_id').prediction,
            part.set_index('sample_id').loc[part70.sample_id,'prediction'],rtol=0,atol=1e-10)
        check(held,item,part70,'prediction');check(hb,item,part70,'baseline_prediction')
    expected=json.loads((ctx.root/'reports/ft30_cv_reference.json').read_text())
    assert all(file_hash(ctx.root/p)==h for p,h in expected.items())
    write_json(ctx.root/'reports/ft30_prediction_validation.json',dict(status='passed',ft30_protocol=protocol,
        n_completed_groups=preds.groupby(['K','model','fold']).ngroups,n_evaluation=len(holdout|adapt),
        adaptation_pool=len(adapt),heldout_evaluation=len(holdout),metrics_recomputed=True,
        primary_source='external_all_post_adaptation',supplement_source='external_heldout_70',
        n_adapted=int(metrics.adaptation_status.eq('adapted').sum()),
        n_fallback=int(metrics.adaptation_status.eq('fallback_no_adaptation').sum()),
        cv_bytes_unchanged=True,no_prediction_filtering=True))

def run_external(ctx,protocol,only=None):
    from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
    import multiprocessing,psutil
    items=selected_folds(ctx)
    if only:items=[i for i in items if i[0]==ctx.cfg['features']['K'][0] and i[1]=='ElasticNet']
    failures=[];active={};todo=iter(items);finished=0;exhausted=False;started=time.time()
    workers=min(OPTIONS['workers_per_lineage'],len(items))
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
        initializer=worker_init,initargs=(ctx.cfg,ctx.lineage,OPTIONS)) as pool:
        while active or not exhausted:
            while len(active)<workers and not exhausted:
                if psutil.virtual_memory().available/2**30<OPTIONS['reserve_memory_gb']:break
                try:item=next(todo)
                except StopIteration:exhausted=True;break
                active[pool.submit(task,ctx.cfg,item,protocol)]=item
            done=wait(active,timeout=OPTIONS['heartbeat_seconds'],return_when=FIRST_COMPLETED)[0] if active else []
            for f in done:
                item=active.pop(f)
                try:f.result()
                except Exception as e:failures.append(dict(item=item,error=str(e),traceback=traceback.format_exc()));LOG.exception('FT30 failed %s',item)
                finished+=1
            state=dict(stage='external_ft30',driver_pid=os.getpid(),updated=time.time(),complete=finished,total=len(items),
                failed=len(failures),workers=workers,running_or_queued=len(active),elapsed_seconds=time.time()-started)
            write_json(ctx.root/'reports/external_ft30_progress.json',state)
            LOG.info('FT30 %s %d/%d failed=%d active=%d',ctx.lineage,finished,len(items),len(failures),len(active))
            if not active and not exhausted:time.sleep(5)
    write_json(ctx.root/'reports/external_ft30_failures.json',failures)
    outputs=collect(ctx,protocol,items,failures)
    if failures:raise RuntimeError(f'{len(failures)} FT30 tasks failed; see reports/external_ft30_failures.json')
    return outputs


def final_external(ctx,protocol):
    """Publish the selected clock's SAME internally selected source fold, not full-refit."""
    decision=json.loads((ctx.root/'stability/selection_decision.json').read_text())
    item=next(i for i in selected_folds(ctx) if i[0]==int(decision['K']) and i[1]==decision['internal_winner'])
    _,target,prefix=source_and_target(ctx,item)
    if not valid(ctx,target/'ft30.complete.json',protocol):raise RuntimeError('Selected clock task incomplete')
    import shutil
    outputs=[]
    for src,dst in [('metrics_ft30.csv','selected_fold_external_ft30_metrics.csv'),
                    ('predictions_ft30.csv','selected_fold_external_ft30_predictions.csv')]:
        p=ctx.root/'external'/dst;shutil.copy2(target/src,p);outputs.append(p)
    return outputs
