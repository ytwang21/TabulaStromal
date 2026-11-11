"""Grouped outer CV and fully nested, train-only linear-model tuning."""
from pathlib import Path
import gc, itertools, json, logging, time, traceback
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold,StratifiedKFold
from .storage import dump,write_csv,write_json
from .metrics import evaluate,METRICS
from .data import load_dataset
from .models import Preprocessor,rank_features,make_estimator,fit_checked,fit_model,predict_model,bags_for,save_bundle
LOG=logging.getLogger(__name__)


def grouped_splits(meta,n_splits=5,seed=42,age_bins=8):
    grouped=meta.groupby('group_id',sort=True).Age.median()
    if len(grouped)<n_splits:raise ValueError(f'{len(grouped)} independent groups < {n_splits} folds')
    bins=pd.qcut(grouped,q=min(age_bins,max(2,len(grouped)//n_splits)),labels=False,duplicates='drop')
    can_stratify=(bins.nunique()>1 and bins.value_counts().min()>=n_splits)
    split=(StratifiedKFold(n_splits,shuffle=True,random_state=seed).split(grouped,bins) if can_stratify else
           KFold(n_splits,shuffle=True,random_state=seed).split(grouped))
    for tr,va in split:
        train_groups=set(grouped.index[tr]);valid_groups=set(grouped.index[va])
        assert train_groups.isdisjoint(valid_groups)
        yield np.flatnonzero(meta.group_id.isin(train_groups)),np.flatnonzero(meta.group_id.isin(valid_groups))


def cached_rank(ctx,ds,rows,tag):
    path=ctx.root/'checkpoints/preprocessing'/tag/'ranking.joblib';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return joblib.load(path)
    rank=rank_features(ds['X'][rows],ds['meta'].Age.to_numpy(float)[rows],ctx.cfg)
    dump(path,rank);ctx.cache.mark(marker,[path]);return rank


def cached_pre(ctx,ds,rows,K,tag,rank=None):
    path=ctx.root/'checkpoints/preprocessing'/tag/f'K{K}.joblib';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return joblib.load(path)
    rank=cached_rank(ctx,ds,rows,tag) if rank is None else rank
    pre=Preprocessor.fit(ds['X'][rows],ds['meta'].Age.to_numpy(float)[rows],ds['schema']['features'],ds['meta'].sample_id.iloc[rows],K,ctx.cfg,rank)
    dump(path,pre);ctx.cache.mark(marker,[path]);return pre


def tune_linear(ctx,ds,rows,K,name,tag):
    c=ctx.cfg['models']
    if name not in ['ElasticNet','Lasso','Ridge'] or not c['tune_linear']:return {}
    path=ctx.root/'checkpoints/tuning'/tag/f'K{K}_{name}.json';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return json.loads(path.read_text())['params']
    candidates=[dict(alpha=a,**({'l1_ratio':r} if name=='ElasticNet' else {})) for a,r in itertools.product(c['linear_alphas'],c['l1_ratios'] if name=='ElasticNet' else [0])]
    scores=np.zeros(len(candidates));failures=[];meta=ds['meta'].iloc[rows].reset_index(drop=True)
    n=min(ctx.cfg['cv']['inner_splits'],meta.group_id.nunique())
    if n<2:raise ValueError('Not enough groups for nested tuning')
    for inner,(itr,iva) in enumerate(grouped_splits(meta,n,ctx.cfg['seed'])):
        tr,va=rows[itr],rows[iva];pre=cached_pre(ctx,ds,tr,K,tag+f'/inner{inner}')
        assert set(pre.fit_ids).isdisjoint(set(ds['meta'].sample_id.iloc[va]))
        X=pre.transform(ds['X'][tr]);V=pre.transform(ds['X'][va]);y=ds['meta'].Age.to_numpy(float)
        for i,params in enumerate(candidates):
            try:
                model=fit_checked(make_estimator(name,ctx.cfg,ctx.cfg['seed'],params),X,y[tr]);scores[i]+=np.mean(np.abs(model.predict(V)-y[va]))/n
            except Exception as e:
                scores[i]=np.inf;failures.append(dict(inner=inner,params=params,error=str(e)))
    if not np.isfinite(scores).any():raise RuntimeError(f'All {name} inner candidates failed: {failures}')
    params=candidates[int(np.argmin(scores))]
    write_json(path,dict(params=params,candidates=candidates,mean_inner_MAE=scores,failures=failures,
                         tuning_ids=meta.sample_id.tolist(),criterion='inner_valid_MAE; full preprocessing refit per inner train'))
    ctx.cache.mark(marker,[path]);return params


def task_path(ctx,K,name,fold):return ctx.root/'models'/f'K{K}'/name/f'fold{fold}'


def benchmark(ctx):
    ds=load_dataset(ctx.root/'features/train');meta=ds['meta'];y=meta.Age.to_numpy(float)
    splits=list(grouped_splits(meta,ctx.cfg['cv']['n_splits'],ctx.cfg['seed'],ctx.cfg['cv']['age_bins']))
    splitrows=[];failures=[]
    for fold,(tr,va) in enumerate(splits,1):
        for split,ix in [('Train',tr),('Valid',va)]:
            splitrows.extend(dict(fold=fold,split=split,sample_id=meta.sample_id.iloc[i],group_id=meta.group_id.iloc[i]) for i in ix)
        assert set(meta.group_id.iloc[tr]).isdisjoint(meta.group_id.iloc[va])
    write_csv(ctx.root/'benchmark/cv_splits.csv',pd.DataFrame(splitrows))
    for fold,(tr,va) in enumerate(splits,1):
        for K in ctx.cfg['features']['K']:
            pre=cached_pre(ctx,ds,tr,K,f'outer{fold}')
            X=pre.transform(ds['X'][tr]);V=pre.transform(ds['X'][va]);trainbags=validbags=None
            for name in ctx.cfg['models']['names']:
                path=task_path(ctx,K,name,fold);marker=path/'complete.json'
                if ctx.cache.valid(marker):continue
                LOG.info('BENCHMARK %s K=%d %s fold=%d, effective K=%d',ctx.lineage,K,name,fold,len(pre.indices));start=time.time()
                try:
                    params=tune_linear(ctx,ds,tr,K,name,f'outer{fold}')
                    if name in ['DeepSets','Ensemble'] and trainbags is None:
                        trainbags=bags_for(ds,pre,tr);validbags=bags_for(ds,pre,va)
                    # Only resume epoch states from this exact source/config/input signature.
                    epochdir=path/'epochs'/ctx.cache.signature;epochdir.mkdir(parents=True,exist_ok=True)
                    model=fit_model(name,X,y[tr],ctx.cfg,ctx.cfg['seed']+fold,params,trainbags,epochdir/'last.pt',ctx.resume)
                    predtr=predict_model(model,name,X,trainbags);predva=predict_model(model,name,V,validbags)
                    info=dict(lineage=ctx.lineage,K=K,fold=fold,seed=ctx.cfg['seed']+fold,params=params,
                              train_ids=meta.sample_id.iloc[tr].tolist(),valid_ids=meta.sample_id.iloc[va].tolist(),
                              train_groups=meta.group_id.iloc[tr].unique().tolist(),valid_groups=meta.group_id.iloc[va].unique().tolist(),
                              training_seconds=time.time()-start,warnings=getattr(model,'fit_warnings_',[]))
                    outputs=save_bundle(path,name,model,pre,info);result=[];preds=[]
                    for split,ix,p in [('Train',tr,predtr),('Valid',va,predva)]:
                        result.append(dict(lineage=ctx.lineage,K=K,model=name,fold=fold,split=split,
                                           effective_K=len(pre.indices),training_seconds=info['training_seconds'],**evaluate(y[ix],p)))
                        table=meta.iloc[ix].copy();table['prediction']=p;table['lineage']=ctx.lineage;table['K']=K;table['model']=name;table['fold']=fold;table['split']=split;preds.append(table)
                    write_csv(path/'metrics.csv',pd.DataFrame(result));write_csv(path/'predictions.csv',pd.concat(preds))
                    ctx.cache.mark(marker,outputs+[path/'metrics.csv',path/'predictions.csv'])
                    del model;gc.collect()
                except Exception as e:
                    LOG.exception('Benchmark failed: %s K=%d fold=%d',name,K,fold)
                    error=dict(K=K,model=name,fold=fold,error=str(e),traceback=traceback.format_exc())
                    failures.append(error);write_json(path/'failure.json',error)
            collect_benchmark(ctx)
    collect_benchmark(ctx)
    write_csv(ctx.root/'benchmark/failed_combinations.csv',pd.DataFrame(failures,columns=['K','model','fold','error','traceback']))
    expected=len(ctx.cfg['features']['K'])*len(ctx.cfg['models']['names'])*len(splits)
    complete=sum(ctx.cache.valid(task_path(ctx,K,n,f)/'complete.json') for K,n,f in itertools.product(ctx.cfg['features']['K'],ctx.cfg['models']['names'],range(1,len(splits)+1)))
    write_json(ctx.root/'benchmark/completeness.json',dict(expected=expected,complete=complete,failed=len(failures)))
    if complete!=expected:raise RuntimeError(f'Benchmark partial: {complete}/{expected}; inspect failed_combinations.csv; rerun --resume')
    return [ctx.root/'benchmark'/s for s in ['cv_fold_metrics.csv','cv_summary.csv','cv_predictions.csv','model_comparison.csv','completeness.json']]


def collect_benchmark(ctx):
    metrics=[];preds=[]
    for p in sorted((ctx.root/'models').glob('K*/*/fold*/complete.json')):
        if ctx.cache.valid(p):
            metrics.append(pd.read_csv(p.parent/'metrics.csv'));preds.append(pd.read_csv(p.parent/'predictions.csv'))
    if not metrics:return
    frame=pd.concat(metrics,ignore_index=True);pred=pd.concat(preds,ignore_index=True)
    write_csv(ctx.root/'benchmark/cv_fold_metrics.csv',frame);write_csv(ctx.root/'benchmark/cv_predictions.csv',pred)
    summary=frame.groupby(['lineage','K','model','split'])[METRICS].agg(['mean','std','count']);summary.columns=['_'.join(c) for c in summary.columns];summary=summary.reset_index()
    write_csv(ctx.root/'benchmark/cv_summary.csv',summary);write_csv(ctx.root/'benchmark/feature_k_comparison.csv',summary)
    valid=summary[summary.split=='Valid'].copy();train=summary[summary.split=='Train'][['K','model','R2_mean']].rename(columns={'R2_mean':'Train_R2_mean'})
    valid=valid.merge(train,on=['K','model']);valid['generalization_gap']=valid.Train_R2_mean-valid.R2_mean
    # Ranking uses internal validation only; gap and variability are prespecified tie breakers.
    valid=valid.sort_values(['R2_mean','MAE_mean','R2_std','generalization_gap'],ascending=[False,True,True,True]);valid['internal_rank']=np.arange(1,len(valid)+1)
    write_csv(ctx.root/'benchmark/model_comparison.csv',valid)
    # Tissue and GSE sensitivity uses held-out OOF rows, not training fit quality.
    checks=[]
    for field in ['GSE','tissue']:
        if field not in pred:continue
        for (K,name,group),sub in pred[pred.split=='Valid'].groupby(['K','model',field]):
            checks.append(dict(K=K,model=name,stratification=field,stratum=group,**evaluate(sub.Age,sub.prediction)))
    write_csv(ctx.root/'benchmark/cv_gse_tissue_sensitivity.csv',pd.DataFrame(checks))
