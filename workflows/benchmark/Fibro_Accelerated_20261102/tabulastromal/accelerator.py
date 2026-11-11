"""Fibro-only acceleration; isolated source snapshot and scientifically checked migration."""
from pathlib import Path
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from functools import lru_cache
import copy,gc,hashlib,itertools,json,logging,multiprocessing,os,shutil,time,traceback,warnings
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import QuantileTransformer
from sklearn.linear_model import Ridge
from sklearn.ensemble import GradientBoostingRegressor
from celer import ElasticNet as CelerEN,Lasso as CelerLasso
from threadpoolctl import threadpool_limits
from . import models,training,external,stability,pipeline
from .storage import Cache,write_json,write_csv,write_text,dump,digest,process_lock,file_hash
from .data import load_dataset
from .metrics import evaluate
LOG=logging.getLogger(__name__)
_ORIGINAL={k:getattr(m,n) for k,m,n in [('fit',models,'fit_model'),('predict',models,'predict_model'),('factory',models,'make_estimator'),('checked',models,'fit_checked'),('bootstrap',stability,'bootstrap_en'),('external_collect',external,'collect_external')]}
_CFG=None
_CACHED_PRE=training.cached_pre


def solver(name,cfg,seed,params=None):
    params=params or {};a=cfg['accelerated']
    if name in ['ElasticNet','Lasso']:
        kw=dict(alpha=params.get('alpha',.1),tol=a['solver_tolerance'],max_iter=a['celer_max_iter'],max_epochs=a['celer_max_epochs'])
        if name=='ElasticNet':kw['l1_ratio']=params.get('l1_ratio',.5)
        return (CelerEN if name=='ElasticNet' else CelerLasso)(**kw)
    return _ORIGINAL['factory'](name,cfg,seed,params)


def checked(model,X,y):
    if not isinstance(model,(CelerEN,CelerLasso)):return _ORIGINAL['checked'](model,X,y)
    notes=[];initial=model.max_iter
    for factor in [1,5,25]:
        model.set_params(max_iter=initial*factor,warm_start=factor>1)
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always');model.fit(np.asfortranarray(X,dtype=np.float64),np.asarray(y,dtype=np.float64))
        notes.extend(str(w.message) for w in ws)
        # Do not silently accept an unconverged candidate merely because a warning was absent.
        gap=float(np.max(np.asarray(getattr(model,'dual_gap_',np.inf))))
        tolerance=model.tol*max(float(np.var(y)),1e-15)
        if np.isfinite(gap) and gap<=tolerance*1.05+1e-9:
            model.fit_warnings_=notes;model.solver_backend_='celer-0.7.4';return model
    raise RuntimeError(f'Celer convergence not attained: dual_gap={gap}, allowed={tolerance}; warnings={notes[-2:]}')


def bulk_fit(self,X,y):
    self.qt=QuantileTransformer(n_quantiles=min(100,len(X)),output_distribution='normal',random_state=self.seed)
    Z=self.qt.fit_transform(X)
    self.models=[solver('ElasticNet',self.cfg,self.seed,dict(alpha=.1,l1_ratio=.5)),Ridge(alpha=10),
                 solver('Lasso',self.cfg,self.seed,dict(alpha=.1)),GradientBoostingRegressor(n_estimators=min(200,self.cfg['models']['boost_estimators']),max_depth=3,learning_rate=.05,subsample=.8,random_state=self.seed)]
    for model in self.models:checked(model,Z,y)
    return self


@lru_cache(maxsize=3)
def dataset(path):return load_dataset(Path(path))


_BAG_CACHE={}
def bags(ds,pre,rows):
    key=(str(ds['path']),hashlib.sha256(pre.indices.tobytes()+pre.scaler.scale_.tobytes()).hexdigest())
    if key not in _BAG_CACHE:
        if len(_BAG_CACHE)>=2:_BAG_CACHE.pop(next(iter(_BAG_CACHE)))
        cells=ds['cells'][:,pre.indices].multiply(1/np.maximum(pre.scaler.scale_,1e-12)).tocsr();ix=ds['bag_index']
        _BAG_CACHE[key]=[(cells[np.flatnonzero(ix['sample']==i)],ix['celltype'][ix['sample']==i]) for i in range(len(ds['meta']))]
    return [_BAG_CACHE[key][int(i)] for i in rows]


@contextmanager
def gpu_lease(cfg):
    import fcntl
    path=Path(cfg['accelerated']['scratch_dir'])/'gpu.lock';path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        yield
        fcntl.flock(lock,fcntl.LOCK_UN)


def fit(name,X,y,cfg,seed,params=None,bags=None,checkpoint=None,resume=True):
    if name not in ['DeepSets','Ensemble'] or not cfg['accelerated']['gpu_enabled']:
        return _ORIGINAL['fit'](name,X,y,cfg,seed,params,bags,checkpoint,resume)
    with gpu_lease(cfg):
        import torch
        if not torch.cuda.is_available():
            LOG.warning('CUDA unavailable in this execution environment; using CPU for %s',name)
            return _ORIGINAL['fit'](name,X,y,cfg,seed,params,bags,checkpoint,resume)
        torch.cuda.set_device(0);torch.cuda.set_per_process_memory_fraction(cfg['accelerated']['gpu_memory_fraction'],0)
        gpu_cfg=copy.deepcopy(cfg);gpu_cfg['resources']['device']='cuda'
        LOG.info('GPU TRAIN %s on physical GPU %s, memory fraction %.2f',name,cfg['accelerated']['gpu_ids'],cfg['accelerated']['gpu_memory_fraction'])
        try:
            model=_ORIGINAL['fit'](name,X,y,gpu_cfg,seed,params,bags,checkpoint,resume)
            deep=model if name=='DeepSets' else model.deep;deep.net.to('cpu');deep.device='cpu';deep.cfg=cfg
            torch.cuda.empty_cache();return model
        except RuntimeError as e:
            if not any(s in str(e).lower() for s in ['cuda','out of memory','driver']):raise
            LOG.exception('GPU runtime failed; resume the same task on CPU')
            torch.cuda.empty_cache()
            return _ORIGINAL['fit'](name,X,y,cfg,seed,params,bags,checkpoint,resume)


def install(cfg):
    global _CFG
    _CFG=cfg
    models.make_estimator=training.make_estimator=solver
    models.fit_checked=training.fit_checked=external.fit_checked=checked
    models.BulkRegressor.fit=bulk_fit
    models.fit_model=training.fit_model=stability.fit_model=external.fit_model=fit
    models.bags_for=training.bags_for=external.bags_for=stability.bags_for=bags
    training.tune_linear=stability.tune_linear=tune_linear
    training.cached_pre=stability.cached_pre=cached_pre
    stability.bootstrap_en=bootstrap_parallel


def worker_init(cfg,thread_override=None):
    threads=thread_override or cfg['accelerated']['worker_threads']
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]=str(threads)
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,cfg['accelerated']['gpu_ids'])) if cfg['accelerated']['gpu_enabled'] else ''
    try:os.nice(5)
    except OSError:LOG.warning('Could not lower worker CPU priority')
    install(cfg);threadpool_limits(limits=threads)


def context(cfg):return pipeline.Context(cfg,'Fibro',True)


def job_log(ctx,kind,key):
    path=ctx.root/'logs/workers'/kind/(str(key)+'.log');path.parent.mkdir(parents=True,exist_ok=True)
    handler=logging.FileHandler(path);handler.setFormatter(logging.Formatter('%(asctime)s | %(levelname)s | %(message)s'))
    logger=logging.getLogger();logger.setLevel(logging.INFO)
    for old in logger.handlers[:]:old.close();logger.removeHandler(old)
    logger.addHandler(handler);logging.captureWarnings(True)
    return path


def run_jobs(ctx,kind,fn,items,workers):
    """A bounded spawn pool; heartbeat is written even during a long native-code solver."""
    if not items:return []
    import psutil
    failures=[];completed=0;root=ctx.root/'reports';active={}
    workers=max(1,min(workers,len(items)));start=time.time()
    thread_override=1 if kind.startswith('bootstrap') else None
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),initializer=worker_init,initargs=(ctx.cfg,thread_override)) as pool:
        todo=iter(items)
        def submit():
            try:item=next(todo)
            except StopIteration:return False
            future=pool.submit(fn,ctx.cfg,item);active[future]=item;return True
        for _ in range(workers*2):
            if not submit():break
        while active:
            done,_=wait(active,timeout=ctx.cfg['accelerated']['heartbeat_seconds'],return_when=FIRST_COMPLETED)
            for future in done:
                item=active.pop(future)
                try:future.result()
                except BaseException as e:
                    error=dict(item=item,error=str(e),traceback=traceback.format_exc());failures.append(error)
                    LOG.error('%s task failed %s: %s',kind,item,e)
                completed+=1
                if psutil.virtual_memory().available/2**30<ctx.cfg['accelerated']['reserve_memory_gb']:
                    LOG.warning('Low available RAM; drain running work before submitting more')
                else:submit()
            if not active:
                for _ in range(workers):
                    if not submit():break
            state=dict(stage=kind,driver_pid=os.getpid(),updated=time.time(),complete=completed,total=len(items),failed=len(failures),running_or_queued=len(active),elapsed_seconds=time.time()-start,workers=workers)
            write_json(root/f'{kind}_progress.json',state);LOG.info('PROGRESS %s %d/%d failed=%d active=%d',kind,completed,len(items),len(failures),len(active))
    write_json(root/f'{kind}_failures.json',failures)
    if failures:raise RuntimeError(f'{kind}: {len(failures)} tasks failed; see {root}/{kind}_failures.json. Successful units remain resumable.')
    return failures


def cached_pre(ctx,ds,rows,K,tag,rank=None):
    # Full-training and final-refit use the identical training rows and inner split seed.
    tag=tag.replace('full_train_final','full_train')
    return _CACHED_PRE(ctx,ds,rows,K,tag,rank)


def prepare_fold(cfg,item):
    fold,inner=item;ctx=context(cfg);ds=dataset(str(ctx.root/'features/train'))
    if fold=='full':tr=np.arange(len(ds['meta']));tag='full_train'
    else:
        tr,va=list(training.grouped_splits(ds['meta'],cfg['cv']['n_splits'],cfg['seed'],cfg['cv']['age_bins']))[int(fold)-1];tag=f'outer{fold}'
    if inner is not None:
        inner_n=min(cfg['cv']['inner_splits'],ds['meta'].iloc[tr].group_id.nunique())
        itr,iva=list(training.grouped_splits(ds['meta'].iloc[tr].reset_index(drop=True),inner_n,cfg['seed']))[inner]
        tr=tr[itr];tag+=f'/inner{inner}'
    job_log(ctx,'prepare',tag.replace('/','_'))
    for K in cfg['features']['K']:training.cached_pre(ctx,ds,tr,K,tag)


def tune_linear(ctx,ds,rows,K,name,tag):
    c=ctx.cfg['models']
    if name not in ['ElasticNet','Lasso','Ridge'] or not c['tune_linear']:return {}
    path=ctx.root/'checkpoints/tuning'/tag/f'K{K}_{name}.json';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return json.loads(path.read_text())['params']
    candidates=[dict(alpha=a,**({'l1_ratio':r} if name=='ElasticNet' else {})) for a,r in itertools.product(c['linear_alphas'],c['l1_ratios'] if name=='ElasticNet' else [0])]
    scores=np.zeros(len(candidates));failures=[];meta=ds['meta'].iloc[rows].reset_index(drop=True);y=ds['meta'].Age.to_numpy(float)
    n=min(ctx.cfg['cv']['inner_splits'],meta.group_id.nunique())
    for inner,(itr,iva) in enumerate(training.grouped_splits(meta,n,ctx.cfg['seed'])):
        tr,va=rows[itr],rows[iva];pre=training.cached_pre(ctx,ds,tr,K,tag+f'/inner{inner}')
        assert set(pre.fit_ids).isdisjoint(set(ds['meta'].sample_id.iloc[va]))
        X=pre.transform(ds['X'][tr]);V=pre.transform(ds['X'][va])
        for i,params in enumerate(candidates):
            task=path.parent/(path.stem+f'_inner{inner}_candidate{i}.json');done=task.with_suffix('.complete.json')
            if ctx.cache.valid(done):r=json.loads(task.read_text())
            else:
                start=time.time();LOG.info('TUNE %s K%d inner%d candidate%d %s START',name,K,inner,i,params)
                try:
                    model=checked(solver(name,ctx.cfg,ctx.cfg['seed'],params),X,y[tr]);mae=float(np.mean(np.abs(model.predict(V)-y[va])))
                    r=dict(MAE=mae,params=params,inner=inner,solver='celer' if name!='Ridge' else 'sklearn',seconds=time.time()-start)
                    write_json(task,r);ctx.cache.mark(done,[task])
                except Exception as e:r=dict(MAE=np.inf,error=str(e));failures.append(dict(inner=inner,params=params,error=str(e)))
            scores[i]+=r['MAE']/n
    if not np.isfinite(scores).any():raise RuntimeError(f'All nested {name} candidates failed: {failures}')
    params=candidates[int(np.argmin(scores))]
    write_json(path,dict(params=params,candidates=candidates,mean_inner_MAE=scores,failures=failures,tuning_ids=meta.sample_id.tolist(),criterion='inner_valid_MAE; preprocessing fits each inner train; per-candidate checkpoints'))
    ctx.cache.mark(marker,[path]);return params


def benchmark_group(cfg,item):
    fold,K=map(int,item);ctx=context(cfg);job_log(ctx,'benchmark',f'fold{fold}_K{K}')
    ds=dataset(str(ctx.root/'features/train'));meta=ds['meta'];y=meta.Age.to_numpy(float)
    tr,va=list(training.grouped_splits(meta,cfg['cv']['n_splits'],cfg['seed'],cfg['cv']['age_bins']))[fold-1]
    pre=training.cached_pre(ctx,ds,tr,K,f'outer{fold}');X=pre.transform(ds['X'][tr]);V=pre.transform(ds['X'][va]);tb=vb=None
    for name in cfg['models']['names']:
        path=training.task_path(ctx,K,name,fold);marker=path/'complete.json'
        if ctx.cache.valid(marker):continue
        start=time.time();LOG.info('BENCHMARK Fibro K=%d %s fold=%d START',K,name,fold)
        params=tune_linear(ctx,ds,tr,K,name,f'outer{fold}')
        if name in ['DeepSets','Ensemble'] and tb is None:tb=bags(ds,pre,tr);vb=bags(ds,pre,va)
        epoch=path/'epochs'/ctx.cache.signature/'last.pt';epoch.parent.mkdir(parents=True,exist_ok=True)
        model=fit(name,X,y[tr],cfg,cfg['seed']+fold,params,tb,epoch,True)
        predtr=models.predict_model(model,name,X,tb);predva=models.predict_model(model,name,V,vb)
        info=dict(lineage='Fibro',K=K,fold=fold,seed=cfg['seed']+fold,params=params,train_ids=meta.sample_id.iloc[tr].tolist(),valid_ids=meta.sample_id.iloc[va].tolist(),
            train_groups=meta.group_id.iloc[tr].unique().tolist(),valid_groups=meta.group_id.iloc[va].unique().tolist(),training_seconds=time.time()-start,
            warnings=getattr(model,'fit_warnings_',[]),engine='accelerated',solver_backend=getattr(model,'solver_backend_','source estimator'))
        outputs=models.save_bundle(path,name,model,pre,info);metrics=[];predictions=[]
        for split,rows,p in [('Train',tr,predtr),('Valid',va,predva)]:
            metrics.append(dict(lineage='Fibro',K=K,model=name,fold=fold,split=split,effective_K=len(pre.indices),training_seconds=info['training_seconds'],**evaluate(y[rows],p)))
            table=meta.iloc[rows].copy();table['prediction']=p;table['lineage']='Fibro';table['K']=K;table['model']=name;table['fold']=fold;table['split']=split;predictions.append(table)
        write_csv(path/'metrics.csv',pd.DataFrame(metrics));write_csv(path/'predictions.csv',pd.concat(predictions));ctx.cache.mark(marker,outputs+[path/'metrics.csv',path/'predictions.csv'])
        write_json(path/'solver_diagnostics.json',dict(backend=getattr(model,'solver_backend_','source estimator'),seconds=time.time()-start,warnings=getattr(model,'fit_warnings_',[])))
        LOG.info('BENCHMARK %s K%d fold%d COMPLETE %.2fs',name,K,fold,time.time()-start);del model;gc.collect()


def benchmark(ctx):
    ds=dataset(str(ctx.root/'features/train'));rows=[]
    for fold,(tr,va) in enumerate(training.grouped_splits(ds['meta'],ctx.cfg['cv']['n_splits'],ctx.cfg['seed'],ctx.cfg['cv']['age_bins']),1):
        for split,ix in [('Train',tr),('Valid',va)]:rows.extend(dict(fold=fold,split=split,sample_id=ds['meta'].sample_id.iloc[i],group_id=ds['meta'].group_id.iloc[i]) for i in ix)
    write_csv(ctx.root/'benchmark/cv_splits.csv',pd.DataFrame(rows))
    preparations=list(itertools.product(list(range(1,ctx.cfg['cv']['n_splits']+1))+['full'],[None]+list(range(ctx.cfg['cv']['inner_splits']))))
    run_jobs(ctx,'prepare',prepare_fold,preparations,ctx.cfg['accelerated']['workers'])
    tasks=list(itertools.product(range(1,ctx.cfg['cv']['n_splits']+1),ctx.cfg['features']['K']))
    try:run_jobs(ctx,'benchmark',benchmark_group,tasks,ctx.cfg['accelerated']['workers'])
    finally:training.collect_benchmark(ctx)
    expected=len(tasks)*len(ctx.cfg['models']['names']);complete=sum(ctx.cache.valid(training.task_path(ctx,K,n,f)/'complete.json') for f,K in tasks for n in ctx.cfg['models']['names'])
    write_json(ctx.root/'benchmark/completeness.json',dict(expected=expected,complete=complete,failed=expected-complete))
    write_csv(ctx.root/'benchmark/failed_combinations.csv',pd.DataFrame(columns=['K','model','fold','error','traceback']))
    if complete!=expected:raise RuntimeError(f'Only {complete}/{expected} models complete')
    return [ctx.root/'benchmark'/f for f in ['cv_fold_metrics.csv','cv_summary.csv','cv_predictions.csv','model_comparison.csv','completeness.json']]


def adapt_fast(name,base,pre,train,ext,tr,ai,ctx,info,path):
    cached=ctx._external_arrays;Xtr,ytr,Xe,ye=cached;Xa,ya=Xe[ai],ye[ai]
    if name in ['ElasticNet','Ridge','Lasso','Pseudobulk']:
        model=fit(name,np.vstack([Xtr,Xa]),np.r_[ytr,ya],ctx.cfg,ctx.cfg['seed'],info.get('params',{}));dump(path/'adapter.joblib',model)
        return model,'refit','joint source+adaptation refit; fixed source preprocessing/hyperparameters; accelerated solver'
    if name in ['DeepSets','Ensemble']:
        model=copy.deepcopy(base);deep=model if name=='DeepSets' else model.deep
        deep.net.weight.requires_grad_(False);deep.net.bias.requires_grad_(False)
        for p in deep.net.encoder.parameters():p.requires_grad_(False)
        # Small adaptation panels have higher GPU launch overhead; use parallel CPU workers.
        deep.cfg=ctx.cfg;deep.device='cpu';deep.net.to('cpu')
        epoch=path/'epochs'/ctx.cache.signature/'adapter_epoch.pt';epoch.parent.mkdir(parents=True,exist_ok=True)
        deep.fit(bags(ext,pre,ai),ya,epoch,True,fine_tune=True);deep.save(path/'adapter_state_dict.pt')
        if name=='Ensemble':model.bulk=checked(models.BulkRegressor(ctx.cfg,ctx.cfg['seed']),np.vstack([Xtr,Xa]),np.r_[ytr,ya]);dump(path/'adapter_bulk.joblib',model.bulk)
        return model,'fine_tune','freeze cell encoder; tune attention/head; Ensemble bulk joint refit; source hyperparameters fixed'
    pa=models.predict_model(base,name,Xa);center=float(pa.mean())
    model=dict(calibrator=Ridge(alpha=ctx.cfg['external']['calibration_alpha']).fit((pa-center)[:,None],ya-pa),prediction_center=center)
    dump(path/'adapter.joblib',model);return model,'calibration','ridge residual slope/intercept; fixed alpha'


def external_group(cfg,item):
    dsid,K,name,fold=item;K=int(K);fold=int(fold);ctx=context(cfg);job_log(ctx,'external',f'{dsid}_K{K}_{name}_fold{fold}')
    it=next(r for r in json.loads((ctx.root/'audit/inventory.json').read_text()) if r['dataset']==dsid)
    ext=dataset(str(ctx.root/'features'/dsid));train=dataset(str(ctx.root/'features/train'));meta=ext['meta'];y=meta.Age.to_numpy(float)
    path=ctx.root/'external/tasks'/dsid/f'K{K}'/name/f'fold{fold}';groupmarker=path/'group_complete.json'
    if ctx.cache.valid(groupmarker):return
    path.mkdir(parents=True,exist_ok=True)
    if not len(meta):
        write_json(path/'excluded.json',dict(dataset=dsid,reason='no nonoverlapping numeric-age samples; see input audit'))
        ctx.cache.mark(groupmarker,[path/'excluded.json']);return
    source=training.task_path(ctx,K,name,fold)
    if not ctx.cache.valid(source/'complete.json'):raise RuntimeError(f'Source task incomplete: {source}')
    _,base,pre,info=models.load_bundle(source,cfg);tr=np.flatnonzero(train['meta'].sample_id.isin(info['train_ids']))
    Xe=pre.transform(ext['X']);Xtr=pre.transform(train['X'][tr]);ctx._external_arrays=(Xtr,train['meta'].Age.to_numpy(float)[tr],Xe,y)
    kind='combined' if dsid=='combined' else 'per_gse'
    bi=dict(lineage='Fibro',K=K,model=name,fold=fold,GSE='__combined__' if kind=='combined' else str(meta.GSE.iloc[0]),dataset=dsid,method=it['method'],source_kind=kind,primary=it['method'] in ['stored_log1p',cfg['external']['primary_method']])
    nop=path/'no_adapt';nop.mkdir(exist_ok=True);marker=nop/'complete.json'
    if ctx.cache.valid(marker):prediction=pd.read_csv(nop/'predictions.csv').prediction.to_numpy()
    else:
        eb=bags(ext,pre,np.arange(len(meta))) if name in ['DeepSets','Ensemble'] else None
        prediction=models.predict_model(base,name,Xe,eb)
        write_csv(nop/'metrics.csv',pd.DataFrame([dict(**bi,adaptation_fraction=0.,repeat=0,strategy='no_adaptation',n_adapt=0,n_test=len(meta),**evaluate(y,prediction))]))
        table=meta.copy();table['prediction']=prediction
        for key,value in bi.items():table[key]=value
        table['adaptation_fraction']=0.;table['repeat']=0;table['strategy']='no_adaptation'
        write_csv(nop/'predictions.csv',table);ctx.cache.mark(marker,[nop/'metrics.csv',nop/'predictions.csv'])
    outputs=[marker];skips=[]
    for fraction,rep in itertools.product(cfg['external']['fractions'],range(cfg['external']['repeats'])):
        ap=path/f'adapt_f{fraction:.1f}_r{rep}';am=ap/'complete.json'
        if ctx.cache.valid(am):outputs.append(am);continue
        too_small=kind=='per_gse' and meta.group_id.nunique()<cfg['external']['min_gse_samples']
        split,reason=external.adaptation_split(meta,fraction,cfg['seed']+rep,cfg['external']['min_adapt_groups'],cfg['external']['min_test_groups'])
        if too_small or split is None:
            skips.append(dict(**bi,adaptation_fraction=fraction,repeat=rep,reason='GSE has fewer than five independent groups' if too_small else reason));continue
        ai,ti=split;assert set(meta.group_id.iloc[ai]).isdisjoint(meta.group_id.iloc[ti]);ap.mkdir(exist_ok=True)
        LOG.info('ADAPT K%d %s fold%d %s f%.1f repeat%d n_adapt%d n_test%d',K,name,fold,dsid,fraction,rep,len(ai),len(ti))
        write_json(ap/'split.json',dict(adapt_ids=meta.sample_id.iloc[ai].tolist(),test_ids=meta.sample_id.iloc[ti].tolist(),adapt_groups=meta.group_id.iloc[ai].unique().tolist(),test_groups=meta.group_id.iloc[ti].unique().tolist(),
           n_groups_adapt=meta.group_id.iloc[ai].nunique(),n_groups_test=meta.group_id.iloc[ti].nunique(),actual_fraction=meta.group_id.iloc[ai].nunique()/meta.group_id.nunique(),seed=cfg['seed']+rep,
           adapt_age=meta.Age.iloc[ai].describe().to_dict(),test_age=meta.Age.iloc[ti].describe().to_dict(),adapt_GSE=meta.GSE.iloc[ai].value_counts().to_dict(),test_GSE=meta.GSE.iloc[ti].value_counts().to_dict()))
        adapter,strategy,detail=adapt_fast(name,base,pre,train,ext,tr,ai,ctx,info,ap)
        if strategy=='calibration':after=prediction[ti]+adapter['calibrator'].predict((prediction[ti]-adapter['prediction_center'])[:,None])
        else:after=models.predict_model(adapter,name,Xe[ti],bags(ext,pre,ti) if name in ['DeepSets','Ensemble'] else None)
        before=evaluate(y[ti],prediction[ti]);after_metrics=evaluate(y[ti],after)
        record=dict(**bi,adaptation_fraction=fraction,repeat=rep,strategy=strategy,strategy_detail=detail,n_adapt=len(ai),n_test=len(ti),small_sample_note='n<20: unstable; repeat SD descriptive' if len(ti)<20 else '',**after_metrics)
        record.update({'baseline_'+m:before[m] for m in before if m in ['R2','MAE','MAD','Age_MAD','MedAE','RMSE','PCC']})
        record.update({'delta_'+m:after_metrics[m]-before[m] for m in ['R2','MAE','MAD','Age_MAD','MedAE','RMSE','PCC']})
        write_csv(ap/'metrics.csv',pd.DataFrame([record]));table=meta.iloc[ti].copy();table['prediction']=after;table['baseline_prediction']=prediction[ti]
        for key,value in bi.items():table[key]=value
        table['adaptation_fraction']=fraction;table['repeat']=rep;table['strategy']=strategy;write_csv(ap/'predictions.csv',table)
        ctx.cache.mark(am,[ap/'metrics.csv',ap/'predictions.csv',ap/'split.json']+list(ap.glob('adapter*')));outputs.append(am);del adapter
    if skips:write_csv(path/'skipped_fractions.csv',pd.DataFrame(skips));outputs.append(path/'skipped_fractions.csv')
    ctx.cache.mark(groupmarker,outputs)


def external_parallel(ctx):
    items=json.loads((ctx.root/'audit/inventory.json').read_text())[1:]
    tasks=list(itertools.product([it['dataset'] for it in items],ctx.cfg['features']['K'],ctx.cfg['models']['names'],range(1,ctx.cfg['cv']['n_splits']+1)))
    # Heavy models first spread across all workers; no external outcome determines task selection.
    order={'Pseudobulk':0,'Ensemble':1,'ElasticNet':2,'Lasso':3,'DeepSets':4}
    tasks.sort(key=lambda t:(order.get(t[2],5),t[0],t[1],t[3]))
    try:run_jobs(ctx,'external',external_group,tasks,ctx.cfg['accelerated']['workers'])
    finally:
        _ORIGINAL['external_collect'](ctx)
        skips=[pd.read_csv(p) for p in (ctx.root/'external/tasks').glob('*/*/*/*/skipped_fractions.csv')]
        old=pd.read_csv(ctx.root/'external/skipped_gse.csv')
        if skips:write_csv(ctx.root/'external/skipped_gse.csv',pd.concat([old]+skips,ignore_index=True))
    write_csv(ctx.root/'external/failed_combinations.csv',pd.DataFrame(columns=['lineage','K','model','fold','GSE','error']))
    return [ctx.root/'external'/f for f in ['external_combined_no_adapt.csv','external_per_gse_no_adapt.csv','external_combined_adaptation_all.csv','external_per_gse_adaptation.csv','skipped_gse.csv']]+list((ctx.root/'external/tasks').glob('*/*/*/*/group_complete.json'))


def bootstrap_one(cfg,item):
    tag,rows,K,params,b=item;ctx=context(cfg);job_log(ctx,'bootstrap',f'{tag}_{b:04d}');ds=dataset(str(ctx.root/'features/train'));rows=np.asarray(rows,int)
    path=ctx.root/'stability'/tag/'bootstrap'/f'boot_{b:04d}.joblib';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return
    meta=ds['meta'].iloc[rows];groups=meta.group_id.unique();rng=np.random.default_rng(cfg['seed']+b);draw=rng.choice(groups,len(groups),replace=True)
    ix=np.concatenate([rows[np.flatnonzero(meta.group_id.to_numpy()==g)] for g in draw]);LOG.info('BOOTSTRAP %s %d started n=%d',tag,b,len(ix))
    rank=models.rank_features(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],cfg)
    pre=models.Preprocessor.fit(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],ds['schema']['features'],ds['meta'].sample_id.iloc[ix],K,cfg,rank)
    model=fit('ElasticNet',pre.transform(ds['X'][ix]),ds['meta'].Age.to_numpy(float)[ix],cfg,cfg['seed']+b,params)
    record=dict(indices=pre.indices,coefficient=np.asarray(model.coef_,np.float32),input_coefficient=np.asarray(model.coef_/pre.scaler.scale_,np.float32),bootstrap=b,draw_groups=draw,fit_ids=pre.fit_ids,params=params,solver='celer')
    dump(path,record);ctx.cache.mark(marker,[path]);LOG.info('BOOTSTRAP %s %d complete',tag,b)


def bootstrap_parallel(ctx,ds,rows,K,params,tag):
    tasks=[(tag,np.asarray(rows).tolist(),K,params,b) for b in range(ctx.cfg['stability']['n_bootstrap'])]
    run_jobs(ctx,'bootstrap_'+tag,bootstrap_one,tasks,ctx.cfg['accelerated']['bootstrap_workers'])
    return _ORIGINAL['bootstrap'](ctx,ds,rows,K,params,tag)


def migrate(ctx):
    """Validate source scientific settings/hashes, then import immutable artifacts only."""
    target=ctx.root/'reports/migration.json'
    if ctx.cache.valid(target.with_suffix('.complete.json')):return json.loads(target.read_text())
    base=Path(ctx.cfg['accelerated']['source_root']);source=base/'Fibro';identity=json.loads((source/'reports/run_identity.json').read_text());old=identity['config'];sc=Cache(identity['signature'],True)
    for key in ['seed','fields','features','cv','models','external','stability','enrichment']:
        if old[key]!=ctx.cfg[key]:raise ValueError(f'Cannot migrate changed scientific setting {key}')
    if old['lineages']['Fibro']!=ctx.cfg['lineages']['Fibro']:raise ValueError('Fibro input paths changed')
    for item in identity['inputs']:
        p=Path(item['path'])
        if not p.exists() or p.stat().st_size!=item.get('size') or p.stat().st_mtime_ns!=item.get('mtime_ns'):raise ValueError(f'Source input fingerprint changed: {p}')
    report=dict(source_root=str(source),source_signature=identity['signature'],new_signature=ctx.cache.signature,models=[],preprocessing=0,tuning=0,
        policy='same inputs, sampling, normalization, splits, K and hyperparameter grids; Celer solves the same convex objective; imported estimators preserve their original solver provenance')
    def copy_file(p,dest):
        dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():return
        # Hard links preserve data bytes without duplicating large cache files; all writes are atomic replacement.
        try:os.link(p,dest)
        except OSError:shutil.copy2(p,dest)
    if not sc.valid(source/'checkpoints/stage_audit.json'):raise ValueError('Source audit checkpoint invalid')
    for p in (source/'audit').iterdir():
        if p.is_file():copy_file(p,ctx.root/'audit'/p.name)
    import yaml
    write_text(ctx.root/'audit/config_accelerated.yaml',yaml.safe_dump(ctx.cfg,sort_keys=False))
    audit_outputs=[ctx.root/'audit'/f for f in ['input_summary.csv','leakage_audit.csv','excluded_samples.csv','inventory.json','audit_report.md','environment.txt']]
    ctx.cache.mark(ctx.root/'checkpoints/stage_audit.json',audit_outputs)
    if not sc.valid(source/'checkpoints/stage_features.json'):raise ValueError('Source feature checkpoint invalid')
    feature_markers=[]
    for old_marker in sorted((source/'features').glob('*/complete.json')):
        if not sc.valid(old_marker):raise ValueError(f'Invalid feature cache {old_marker}')
        dest=ctx.root/'features'/old_marker.parent.name;outputs=[]
        for p in old_marker.parent.iterdir():
            if p.is_file() and p.name!='complete.json':copy_file(p,dest/p.name);outputs.append(dest/p.name)
        marker=dest/'complete.json';ctx.cache.mark(marker,outputs);feature_markers.append(marker)
    ctx.cache.mark(ctx.root/'checkpoints/stage_features.json',feature_markers)
    for sub in ['preprocessing','tuning']:
        for old_marker in (source/'checkpoints'/sub).rglob('*.complete.json'):
            if not sc.valid(old_marker):continue
            record=json.loads(old_marker.read_text());outputs=[]
            for f in record['outputs']:
                p=Path(f['path']);dest=ctx.root/p.relative_to(source);copy_file(p,dest);outputs.append(dest)
            ctx.cache.mark(ctx.root/old_marker.relative_to(source),outputs);report[sub]+=1
    for marker in sorted((source/'models').glob('K*/*/fold*/complete.json')):
        if not sc.valid(marker):continue
        dest=ctx.root/marker.parent.relative_to(source);outputs=[]
        record=json.loads(marker.read_text())
        for f in record['outputs']:
            p=Path(f['path']);target_file=dest/p.name;copy_file(p,target_file);outputs.append(target_file)
        write_json(dest/'import_provenance.json',dict(source=str(marker.parent),source_signature=identity['signature'],preserved_original_model=True))
        ctx.cache.mark(dest/'complete.json',outputs+[dest/'import_provenance.json']);report['models'].append(str(dest.relative_to(ctx.root)))
    report['imported_model_count']=len(report['models']);write_json(target,report);ctx.cache.mark(target.with_suffix('.complete.json'),[target])
    training.collect_benchmark(ctx);return report


def stage(ctx,name):
    if name=='benchmark':return benchmark(ctx)
    if name=='external':return external_parallel(ctx)
    mapping={'audit':('data','audit'),'features':('data','features'),'stability':('stability','stability'),'enrichment':('biology','enrichment'),'figures':('figures','figures')}
    import importlib
    module,fn=mapping[name];return getattr(importlib.import_module('tabulastromal.'+module),fn)(ctx)
