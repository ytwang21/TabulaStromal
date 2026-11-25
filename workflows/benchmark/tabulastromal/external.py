"""Fixed-pipeline external evaluation and disjoint, paired adaptation experiments."""
from pathlib import Path
import copy, gc, itertools, json, logging, traceback
import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from .storage import write_csv,write_json,dump,Cache
from .metrics import evaluate,METRICS
from .data import load_dataset
from .models import load_bundle,bags_for,predict_model,fit_model,DeepRegressor,BulkRegressor,fit_checked
from .training import task_path
LOG=logging.getLogger(__name__)



def external_scope(ctx):
    """Reload only the external scheduling policy at phase entry, including live drivers.

    The matching output root/revision prevents a synthetic or alternate run from
    inheriting this project's live policy. No scientific/model parameter is reloaded.
    """
    import yaml
    resources=ctx.cfg['resources']
    config=Path(resources.get('external_policy_config',Path(__file__).resolve().parents[1]/'config.yaml'))
    scope=resources.get('external_scope','all')
    if config.is_file():
        live=yaml.safe_load(config.read_text())
        if (live.get('output_root')==ctx.cfg.get('output_root') and
                live.get('execution_revision')==ctx.cfg.get('execution_revision')):
            scope=live.get('resources',{}).get('external_scope',scope)
    if scope not in ('all','combined'):
        raise ValueError(f'Unknown external_scope: {scope!r}')
    resources['external_scope']=scope
    return scope


def selected_external_items(ctx,items=None):
    items=json.loads((ctx.root/'audit/inventory.json').read_text())[1:] if items is None else items
    scope=external_scope(ctx)
    selected=[it for it in items if scope=='all' or it['dataset']=='combined']
    if scope=='combined' and not selected:
        raise ValueError('combined-only external policy has no combined input')
    return selected


def record_external_policy(ctx,items):
    """Persist actual phase policy without changing retained training signatures."""
    from .storage import file_hash
    root=Path(__file__).resolve().parents[1]
    code={str(p.relative_to(root)):file_hash(p) for p in sorted((root/'tabulastromal').glob('*.py'))}
    code['run_tabulastromal.py']=file_hash(root/'run_tabulastromal.py')
    path=ctx.root/'reports/run_identity.json'
    if path.is_file():
        identity=json.loads(path.read_text())
        identity.setdefault('execution',{})['external_phase_policy']={
            'scope':ctx.cfg['resources']['external_scope'],
            'datasets':[it['dataset'] for it in items],
            'per_gse_deferred':ctx.cfg['resources']['external_scope']=='combined',
            'phase_sources':code,
            'checkpoint_signature_retained':ctx.cache.signature}
        write_json(path,identity)
    LOG.info('EXTERNAL POLICY scope=%s datasets=%s; benchmark checkpoints retained',
             ctx.cfg['resources']['external_scope'],[it['dataset'] for it in items])


def adaptation_split(meta,fraction,seed,min_adapt=2,min_test=2):
    """Group-disjoint split; deterministic age/GSE/tissue-balanced candidate sampling."""
    groups=meta.groupby('group_id').agg(Age=('Age','median'),GSE=('GSE','first'))
    if 'tissue' in meta:groups['tissue']=meta.groupby('group_id').tissue.first()
    n=len(groups);n_adapt=int(round(n*fraction))
    if n_adapt<min_adapt or n-n_adapt<min_test:
        return None,f'fraction gives {n_adapt} adapt and {n-n_adapt} test groups; minimum {min_adapt}/{min_test}'
    rng=np.random.default_rng(seed);best=None;best_score=np.inf
    age=groups.Age.to_numpy();age_sd=max(float(np.std(age)),1)
    for _ in range(32):
        a=rng.choice(n,n_adapt,replace=False);t=np.setdiff1d(np.arange(n),a)
        imbalance=abs(age[a].mean()-age.mean())/age_sd+abs(age[t].mean()-age.mean())/age_sd
        for col in ['GSE']+(['tissue'] if 'tissue' in groups else []):
            full=groups[col].value_counts(normalize=True);aa=groups.iloc[a][col].value_counts(normalize=True).reindex(full.index,fill_value=0)
            imbalance+=float(np.abs(full-aa).sum())
        if imbalance<best_score:best=(a,t);best_score=imbalance
    ag=set(groups.index[best[0]]);tg=set(groups.index[best[1]])
    assert ag.isdisjoint(tg)
    return (np.flatnonzero(meta.group_id.isin(ag)),np.flatnonzero(meta.group_id.isin(tg))),''


def adapt_model(name,model,pre,train,ext,tr,ai,ctx,info,path):
    Xtr=pre.transform(train['X'][tr]);Xa=pre.transform(ext['X'][ai]);ytr=train['meta'].Age.to_numpy(float)[tr];ya=ext['meta'].Age.to_numpy(float)[ai]
    if name in ['ElasticNet','Ridge','Lasso','IdealAge(Pseudobulk)']:
        # Refit with source+adaptation observations, fixed source preprocessing and hyperparameters.
        adapted=fit_model(name,np.vstack([Xtr,Xa]),np.r_[ytr,ya],ctx.cfg,ctx.cfg['seed'],info.get('params',{}))
        dump(path/'adapter.joblib',adapted)
        return adapted,'refit','joint source+adaptation refit; source hyperparameters fixed'
    if name in ['IdealAge(DeepSets)','IdealAge(Ensemble)']:
        adapted=copy.deepcopy(model);deep=adapted if name=='IdealAge(DeepSets)' else adapted.deep
        deep.net.weight.requires_grad_(False);deep.net.bias.requires_grad_(False)
        for p in deep.net.encoder.parameters():p.requires_grad_(False)
        epoch=path/'epochs'/ctx.cache.signature/'adapter_epoch.pt';epoch.parent.mkdir(parents=True,exist_ok=True)
        deep.fit(bags_for(ext,pre,ai),ya,epoch,ctx.resume,fine_tune=True)
        deep.save(path/'adapter_state_dict.pt')
        if name=='IdealAge(Ensemble)':
            adapted.bulk=fit_checked(BulkRegressor(ctx.cfg,ctx.cfg['seed']),np.vstack([Xtr,Xa]),np.r_[ytr,ya]);dump(path/'adapter_bulk.joblib',adapted.bulk)
        return adapted,'fine_tune','freeze cell encoder; update subtype attention+head; Ensemble bulk jointly refit'
    # Nonincremental trees/MLP/PCA/SVR: explicit residual slope/intercept calibration.
    # Fit only on adaptation labels. It is never called fine-tuning.
    pa=predict_model(model,name,Xa);center=float(pa.mean())
    cal=Ridge(alpha=ctx.cfg['external']['calibration_alpha']).fit((pa-center)[:,None],ya-pa)
    adapter=dict(calibrator=cal,prediction_center=center);dump(path/'adapter.joblib',adapter)
    return adapter,'calibration','ridge residual slope/intercept; alpha fixed in config'


def predict_adapted(adapter,strategy,base,name,pre,ext,rows,ctx):
    X=pre.transform(ext['X'][rows]);bags=bags_for(ext,pre,rows) if name in ['IdealAge(DeepSets)','IdealAge(Ensemble)'] else None
    if strategy=='calibration':
        pred=predict_model(base,name,X,bags);return pred+adapter['calibrator'].predict((pred-adapter['prediction_center'])[:,None])
    return predict_model(adapter,name,X,bags)


def external(ctx,items_override=None,collect=True,combinations_override=None):
    if not collect:
        from .storage import worker_logging
        worker_logging(ctx.root,'external')
    items=selected_external_items(ctx,items_override)
    if collect:record_external_policy(ctx,items)
    train=load_dataset(ctx.root/'features/train')
    failures=[]
    if collect and ctx.cfg['resources'].get('external_workers',1)>1:
        from joblib import Parallel,delayed,parallel_config
        with parallel_config(backend='loky',inner_max_num_threads=ctx.cfg['resources']['threads']):
            from .pipeline import task_context
            # Shard the combined cohort by K/fold so dataset filtering retains
            # parallel execution. Every worker owns disjoint checkpoint paths.
            jobs=[(it,[(K,name,fold) for name in ctx.cfg['models']['names']])
                  for it in items for K,fold in itertools.product(
                      ctx.cfg['features']['K'],range(1,ctx.cfg['cv']['n_splits']+1))]
            batches=Parallel(n_jobs=ctx.cfg['resources']['external_workers'],pre_dispatch=ctx.cfg['resources']['external_workers'])(
                delayed(external)(task_context(ctx,i),[it],False,combos)
                for i,(it,combos) in enumerate(jobs))
        failures=[x for batch in batches for x in batch]
        collect_external(ctx)
        write_csv(ctx.root/'external/failed_combinations.csv',pd.DataFrame(failures,columns=['lineage','K','model','fold','error']))
        if failures:raise RuntimeError(f'{len(failures)} external tasks failed; see failed_combinations.csv')
        return [ctx.root/'external'/f for f in ['external_combined_no_adapt.csv','external_per_gse_no_adapt.csv','external_combined_adaptation_all.csv','external_per_gse_adaptation.csv','skipped_gse.csv']]
    for it in items:
        ext=load_dataset(ctx.root/'features'/it['dataset']);meta=ext['meta'];y=meta.Age.to_numpy(float);kind='combined' if it['dataset']=='combined' else 'per_gse'
        if len(meta)==0:
            write_json(ctx.root/'external/tasks'/it['dataset']/'excluded.json',dict(reason='no nonoverlapping eligible samples'));continue
        combinations=combinations_override if combinations_override is not None else itertools.product(ctx.cfg['features']['K'],ctx.cfg['models']['names'],range(1,ctx.cfg['cv']['n_splits']+1))
        for K,name,fold in combinations:
            source=task_path(ctx,K,name,fold)
            if not Cache(ctx.cache.signature,True).valid(source/'complete.json'):raise RuntimeError(f'Incomplete benchmark task: {source}')
            basepath=ctx.root/'external/tasks'/it['dataset']/f'K{K}'/name/f'fold{fold}';basepath.mkdir(parents=True,exist_ok=True)
            _,model,pre,info=load_bundle(source,ctx.cfg)
            tr=np.flatnonzero(train['meta'].sample_id.isin(info['train_ids']))
            baseinfo=dict(lineage=ctx.lineage,K=K,model=name,fold=fold,GSE='__combined__' if kind=='combined' else str(meta.GSE.iloc[0]),
                          dataset=it['dataset'],method=it['method'],source_kind=kind,primary=it['method'] in ['stored_log1p',ctx.cfg['external']['primary_method']])
            nop=basepath/'no_adapt';nop.mkdir(exist_ok=True)
            try:
                if ctx.cache.valid(nop/'complete.json'):pred=pd.read_csv(nop/'predictions.csv').prediction.to_numpy()
                else:
                    LOG.info('EXTERNAL %s K=%d %s fold=%d n=%d',it['dataset'],K,name,fold,len(meta))
                    bags=bags_for(ext,pre,np.arange(len(meta))) if name in ['IdealAge(DeepSets)','IdealAge(Ensemble)'] else None
                    pred=predict_model(model,name,pre.transform(ext['X']),bags)
                    record=dict(**baseinfo,adaptation_fraction=0.,repeat=0,strategy='no_adaptation',n_adapt=0,n_test=len(meta),**evaluate(y,pred))
                    write_csv(nop/'metrics.csv',pd.DataFrame([record]));predtable=meta.copy();predtable['prediction']=pred
                    for key,value in baseinfo.items():predtable[key]=value
                    predtable['adaptation_fraction']=0.;predtable['repeat']=0;predtable['strategy']='no_adaptation'
                    write_csv(nop/'predictions.csv',predtable);ctx.cache.mark(nop/'complete.json',[nop/'metrics.csv',nop/'predictions.csv'])
                for fraction,rep in itertools.product(ctx.cfg['external']['fractions'],range(ctx.cfg['external']['repeats'])):
                    path=basepath/f'adapt_f{fraction:.1f}_r{rep}';marker=path/'complete.json'
                    if ctx.cache.valid(marker):continue
                    path.mkdir(exist_ok=True)
                    skip='GSE has fewer than five independent groups' if kind=='per_gse' and meta.group_id.nunique()<ctx.cfg['external']['min_gse_samples'] else ''
                    split,reason=adaptation_split(meta,fraction,ctx.cfg['seed']+rep,min_adapt=ctx.cfg['external']['min_adapt_groups'],min_test=ctx.cfg['external']['min_test_groups'])
                    if skip or split is None:
                        write_json(path/'skipped.json',dict(**baseinfo,adaptation_fraction=fraction,repeat=rep,reason=skip or reason))
                        ctx.cache.mark(marker,[path/'skipped.json'],outcome='skipped');continue
                    ai,ti=split;assert set(meta.group_id.iloc[ai]).isdisjoint(meta.group_id.iloc[ti])
                    write_json(path/'split.json',dict(adapt_ids=meta.sample_id.iloc[ai].tolist(),test_ids=meta.sample_id.iloc[ti].tolist(),
                        adapt_groups=meta.group_id.iloc[ai].unique().tolist(),test_groups=meta.group_id.iloc[ti].unique().tolist(),
                        n_groups_adapt=meta.group_id.iloc[ai].nunique(),n_groups_test=meta.group_id.iloc[ti].nunique(),
                        actual_fraction=meta.group_id.iloc[ai].nunique()/meta.group_id.nunique(),seed=ctx.cfg['seed']+rep,
                        adapt_age=meta.Age.iloc[ai].describe().to_dict(),test_age=meta.Age.iloc[ti].describe().to_dict(),
                        adapt_GSE=meta.GSE.iloc[ai].value_counts().to_dict(),test_GSE=meta.GSE.iloc[ti].value_counts().to_dict()))
                    adapter,strategy,detail=adapt_model(name,model,pre,train,ext,tr,ai,ctx,info,path)
                    after=predict_adapted(adapter,strategy,model,name,pre,ext,ti,ctx);before=pred[ti]
                    bm=evaluate(y[ti],before);am=evaluate(y[ti],after)
                    record=dict(**baseinfo,adaptation_fraction=fraction,repeat=rep,strategy=strategy,strategy_detail=detail,
                                n_adapt=len(ai),n_test=len(ti),small_sample_note='held-out n<20: unstable; repeat SD is descriptive' if len(ti)<20 else '',**am)
                    record.update({'baseline_'+m:bm[m] for m in METRICS});record.update({'delta_'+m:am[m]-bm[m] for m in METRICS})
                    write_csv(path/'metrics.csv',pd.DataFrame([record]));table=meta.iloc[ti].copy();table['prediction']=after;table['baseline_prediction']=before
                    for key,value in baseinfo.items():table[key]=value
                    table['adaptation_fraction']=fraction;table['repeat']=rep;table['strategy']=strategy
                    write_csv(path/'predictions.csv',table)
                    adapter_files=list(path.glob('adapter*'))
                    ctx.cache.mark(marker,[path/'metrics.csv',path/'predictions.csv',path/'split.json']+adapter_files)
                    del adapter
            except Exception as e:
                LOG.exception('External task failed');fail=dict(**baseinfo,error=str(e),traceback=traceback.format_exc());failures.append(fail);write_json(basepath/'failure.json',fail)
            del model;gc.collect()
        if collect:collect_external(ctx)
        del ext;gc.collect()
    if not collect:return failures
    collect_external(ctx);write_csv(ctx.root/'external/failed_combinations.csv',pd.DataFrame(failures,columns=['lineage','K','model','fold','GSE','dataset','method','source_kind','primary','error','traceback']))
    if failures:raise RuntimeError(f'{len(failures)} external tasks failed; completed tasks retained')
    return [ctx.root/'external'/f for f in ['external_combined_no_adapt.csv','external_per_gse_no_adapt.csv','external_combined_adaptation_all.csv','external_per_gse_adaptation.csv','skipped_gse.csv']]


def collect_external(ctx):
    tables={};preds={};skips=[]
    scope=external_scope(ctx)
    for p in sorted((ctx.root/'external/tasks').glob('*/*/*/*/*/complete.json')):
        if scope=='combined' and p.relative_to(ctx.root/'external/tasks').parts[0]!='combined':continue
        if not ctx.cache.valid(p):continue
        path=p.parent
        if (path/'skipped.json').exists():skips.append(json.loads((path/'skipped.json').read_text()));continue
        if not (path/'metrics.csv').exists():continue
        t=pd.read_csv(path/'metrics.csv');kind=t.source_kind.iloc[0];adapt=path.name!='no_adapt';key=(kind,adapt)
        tables.setdefault(key,[]).append(t);preds.setdefault(key,[]).append(pd.read_csv(path/'predictions.csv'))
    for kind in ['combined','per_gse']:
        no=pd.concat(tables.get((kind,False),[]),ignore_index=True) if tables.get((kind,False)) else pd.DataFrame(columns=['lineage','K','model','fold','GSE','dataset','method']+METRICS)
        ada=pd.concat(tables.get((kind,True),[]),ignore_index=True) if tables.get((kind,True)) else pd.DataFrame(columns=list(no.columns)+['baseline_MAE','delta_MAE'])
        write_csv(ctx.root/f'external/external_{kind}_no_adapt.csv',no)
        adaptfile = ctx.root/('external/external_combined_adaptation_all.csv' if kind=='combined' else 'external/external_per_gse_adaptation.csv')
        write_csv(adaptfile,pd.concat([no,ada],ignore_index=True) if kind=='combined' else ada)
        for adapt in [False,True]:
            tab=pd.concat(preds.get((kind,adapt),[]),ignore_index=True) if preds.get((kind,adapt)) else pd.DataFrame(columns=['sample_id','Age','prediction'])
            write_csv(ctx.root/f'external/external_{kind}_{"adaptation" if adapt else "no_adapt"}_predictions.csv',tab)
    write_csv(ctx.root/'external/skipped_gse.csv',pd.DataFrame(skips,columns=['lineage','K','model','fold','GSE','dataset','method','source_kind','primary','adaptation_fraction','repeat','reason']))
    # External-selected fold is saved strictly as exploratory, never used in final selection.
    nofile=ctx.root/'external/external_combined_no_adapt.csv'
    t=pd.read_csv(nofile)
    if len(t):
        best=t.sort_values(['R2','MAE'],ascending=[False,True]).groupby(['K','model']).head(1).copy();best['selection_label']='exploratory_external_selected_fold'
        write_csv(ctx.root/'external/exploratory_external_selected_fold.csv',best)
