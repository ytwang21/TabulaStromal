"""Endo runtime acceleration adapter; original pipeline architecture remains intact."""
from pathlib import Path
from functools import lru_cache
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED,ThreadPoolExecutor
from collections import deque
import copy,csv,gc,hashlib,itertools,json,logging,multiprocessing,os,shutil,tempfile,time,traceback,warnings
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import QuantileTransformer
from sklearn.linear_model import Ridge
from sklearn.linear_model import ElasticNet, Lasso
from sklearn.ensemble import GradientBoostingRegressor
from threadpoolctl import threadpool_limits
from tabulastromal import models,training,external,stability,pipeline
from tabulastromal.storage import Cache,write_json,write_csv,write_text,dump,file_hash,atomic_path
from tabulastromal.data import load_dataset
from tabulastromal.metrics import evaluate,METRICS
LOG=logging.getLogger(__name__)
RUNTIME={}
_BAG_CACHE={}
_BOOTSTRAP_CONTEXT={}
_ORIGINAL={k:getattr(m,n) for k,m,n in [('fit',models,'fit_model'),('factory',models,'make_estimator'),('checked',models,'fit_checked'),('bootstrap',stability,'bootstrap_en'),('deep_fit',models.DeepRegressor,'fit'),('collect',external,'collect_external'),('stage',pipeline,'run_stage'),('tune',training,'tune_linear'),('bulk_fit',models.BulkRegressor,'fit')]}

@lru_cache(maxsize=3)
def dataset(path):return load_dataset(Path(path))

def context(cfg):return pipeline.Context(cfg,'Endo' if 'Endo' in cfg['lineages'] else 'Fibro',True)


def worker_init(cfg, options, thread_override=None):
    threads=thread_override or options['worker_threads']
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key]=str(threads)
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,options['gpu_ids'])) if options['gpu_enabled'] else ''
    try:os.nice(5)
    except OSError:pass
    install(cfg, options)
    threadpool_limits(limits=threads)


@contextmanager
def gpu_slot():
    """Use a free GPU lease or fall back to CPU; never queue CPU tasks on a GPU lock."""
    import fcntl
    if not RUNTIME['gpu_enabled']:
        yield None
        return
    import torch
    if not torch.cuda.is_available():
        yield None
        return
    root=Path(RUNTIME['scratch_dir']);root.mkdir(parents=True,exist_ok=True)
    handle=None;choice=None
    for logical,physical in enumerate(RUNTIME['gpu_ids']):
        candidate=(root/f'gpu_{physical}.lock').open('a+')
        try:fcntl.flock(candidate,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            candidate.close();continue
        try:
            free,_=torch.cuda.mem_get_info(logical)
            if free/2**30 < RUNTIME['gpu_min_free_gb']:
                candidate.close();continue
            handle,choice=candidate,logical
            break
        except RuntimeError:
            candidate.close()
    try:yield choice
    finally:
        if handle is not None:handle.close()


def deep_fit(self, bags, y, checkpoint=None, resume=True, fine_tune=False):
    def cpu_fit():
        original_cfg=self.cfg
        self.cfg=copy.deepcopy(original_cfg)
        self.cfg['resources']['threads']=RUNTIME['worker_threads']
        try:return _ORIGINAL['deep_fit'](self,bags,y,checkpoint,resume,fine_tune)
        finally:self.cfg=original_cfg
    use_gpu=RUNTIME['gpu_enabled'] and (not fine_tune or
        (RUNTIME['gpu_adaptation'] and len(y)>=RUNTIME['gpu_min_samples']))
    if not use_gpu:
        return cpu_fit()
    import torch
    with gpu_slot() as device:
        if device is None:
            return cpu_fit()
        old_cfg,old_device=self.cfg,self.device
        torch.cuda.set_device(device)
        torch.cuda.set_per_process_memory_fraction(RUNTIME['gpu_memory_fraction'],device)
        self.device=f'cuda:{device}'
        self.cfg=copy.deepcopy(old_cfg)
        self.cfg['resources']['device']=self.device
        self.cfg['resources']['threads']=RUNTIME['worker_threads']
        if hasattr(self,'net'):self.net.to(self.device)
        LOG.info('GPU DeepSets physical=%s samples=%d fine_tune=%s',RUNTIME['gpu_ids'][device],len(y),fine_tune)
        try:
            return _ORIGINAL['deep_fit'](self,bags,y,checkpoint,resume,fine_tune)
        finally:
            if hasattr(self,'net'):self.net.to(old_device)
            self.cfg,self.device=old_cfg,old_device
            torch.cuda.empty_cache()


def fit(*args,**kwargs):
    return _ORIGINAL['fit'](*args,**kwargs)


def fit_checked(model, X, y):
    """Keep the objective/tolerance/budget; improve unfinished coordinate descent.

    Float32 coordinate descent can stall above the dual-gap tolerance even after
    millions of iterations. Casting the already preprocessed values changes only
    solver arithmetic, not feature selection, scaling or the fitting objective.
    Completed model/adapter checkpoints are never passed through this function.
    """
    if isinstance(model, (ElasticNet, Lasso)) and RUNTIME.get('linear_float64', True):
        model.set_params(selection=RUNTIME.get('linear_selection', 'random'),
                         random_state=RUNTIME.get('linear_random_state', 42))
        X = X.astype(np.float64) if sp.issparse(X) else np.asarray(X, dtype=np.float64, order='F')
        y = np.asarray(y, dtype=np.float64)
        if getattr(model, 'warm_start', False) and hasattr(model, 'coef_'):
            model.coef_ = np.asarray(model.coef_, dtype=np.float64)
            model.intercept_ = np.asarray(model.intercept_, dtype=np.float64)[()]
    result = _ORIGINAL['checked'](model, X, y)
    if isinstance(model, (ElasticNet, Lasso)):
        result.numerical_execution_policy_ = dict(dtype=str(X.dtype), selection=model.selection,
            random_state=model.random_state, objective='unchanged', tolerance=model.tol)
    return result


def prepare_tuning(cfg,item):
    tag,K,inner,tr,va=item
    ctx=context(cfg);job_log(ctx,'tuning',f'{tag.replace("/","_")}_K{K}_inner{inner}_prepare')
    ds=dataset(str(ctx.root/'features/train'))
    training.cached_pre(ctx,ds,np.asarray(tr,int),K,tag+f'/inner{inner}')


def tuning_candidate(cfg,item):
    tag,K,name,inner,index,params,tr,va=item
    ctx=context(cfg);job_log(ctx,'tuning',f'{tag.replace("/","_")}_K{K}_{name}_{inner}_{index}')
    path=ctx.root/'checkpoints/tuning'/tag/f'K{K}_{name}_inner{inner}_candidate{index}.json'
    marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return
    ds=dataset(str(ctx.root/'features/train'));tr=np.asarray(tr,int);va=np.asarray(va,int)
    pre=training.cached_pre(ctx,ds,tr,K,tag+f'/inner{inner}')
    assert set(pre.fit_ids).isdisjoint(set(ds['meta'].sample_id.iloc[va]))
    y=ds['meta'].Age.to_numpy(float)
    try:
        estimator=models.make_estimator(name,cfg,cfg['seed'],params)
        estimator=models.fit_checked(estimator,pre.transform(ds['X'][tr]),y[tr])
        record=dict(MAE=float(np.mean(abs(estimator.predict(pre.transform(ds['X'][va]))-y[va]))),
                    params=params,inner=inner,solver='original_sklearn')
    except Exception as exc:
        record=dict(MAE=float('inf'),params=params,inner=inner,error=str(exc),solver='original_sklearn')
    write_json(path,record);ctx.cache.mark(marker,[path])


def tune_linear(ctx,ds,rows,K,name,tag):
    c=ctx.cfg['models']
    if name not in ['ElasticNet','Lasso','Ridge'] or not c['tune_linear']:return {}
    path=ctx.root/'checkpoints/tuning'/tag/f'K{K}_{name}.json';marker=path.with_suffix('.complete.json')
    if ctx.cache.valid(marker):return json.loads(path.read_text())['params']
    candidates=[dict(alpha=a,**({'l1_ratio':r} if name=='ElasticNet' else {}))
                for a,r in itertools.product(c['linear_alphas'],c['l1_ratios'] if name=='ElasticNet' else [0])]
    meta=ds['meta'].iloc[rows].reset_index(drop=True)
    n=min(ctx.cfg['cv']['inner_splits'],meta.group_id.nunique())
    if n<2:raise ValueError('Not enough groups for nested tuning')
    splits=[(inner,rows[tr],rows[va]) for inner,(tr,va) in enumerate(training.grouped_splits(meta,n,ctx.cfg['seed']))]
    run_jobs(ctx,'tuning_prepare_'+tag.replace('/','_'),prepare_tuning,
             [(tag,K,i,tr.tolist(),va.tolist()) for i,tr,va in splits],min(n,RUNTIME['workers']))
    tasks=[(tag,K,name,i,j,params,tr.tolist(),va.tolist())
           for i,tr,va in splits for j,params in enumerate(candidates)]
    run_jobs(ctx,'tuning_'+tag.replace('/','_'),tuning_candidate,tasks,RUNTIME['workers'])
    scores=np.zeros(len(candidates));failures=[]
    for i,tr,va in splits:
        for j,params in enumerate(candidates):
            task=path.parent/f'K{K}_{name}_inner{i}_candidate{j}.json'
            record=json.loads(task.read_text());scores[j]+=record['MAE']/n
            if 'error' in record:failures.append(dict(inner=i,params=params,error=record['error']))
    if not np.isfinite(scores).any():raise RuntimeError(f'All {name} inner candidates failed: {failures}')
    params=candidates[int(np.argmin(scores))]
    write_json(path,dict(params=params,candidates=candidates,mean_inner_MAE=scores,failures=failures,
               tuning_ids=meta.sample_id.tolist(),criterion='inner_valid_MAE; original preprocessing and sklearn estimator, parallel execution'))
    ctx.cache.mark(marker,[path]);return params


def install(cfg, options):
    global RUNTIME
    RUNTIME=dict(options)
    if RUNTIME['solver']!='sklearn':
        raise ValueError('This execution adapter preserves the original sklearn solvers')
    models.fit_checked=training.fit_checked=external.fit_checked=fit_checked
    models.DeepRegressor.fit=deep_fit
    models.bags_for=training.bags_for=external.bags_for=stability.bags_for=bags
    stability.bootstrap_en=bootstrap_parallel
    training.tune_linear=stability.tune_linear=tune_linear
    def stage(ctx,name):
        if name=='external':
            import endo_external_parallel
            return endo_external_parallel.external_consistent(ctx)
        return _ORIGINAL['stage'](ctx,name)
    pipeline.run_stage=stage


def collect_streaming(ctx):
    """Same CSV products, with bounded memory, concurrent checksum I/O and heartbeats."""
    start=time.time();last=start;count=0
    markers=sorted((ctx.root/'external/tasks').glob('*/*/*/*/*/complete.json'))
    dimensions=['lineage','K','model','fold','GSE','dataset','method','source_kind','primary']
    no_columns=dimensions+['adaptation_fraction','repeat','strategy','n_adapt','n_test']+list(evaluate([0,1],[0,1]))
    ada_columns=dimensions+['adaptation_fraction','repeat','strategy','strategy_detail','n_adapt','n_test','small_sample_note']+list(evaluate([0,1],[0,1]))+['baseline_'+m for m in METRICS]+['delta_'+m for m in METRICS]
    all_columns=no_columns+[c for c in ada_columns if c not in no_columns]
    meta_columns=[]
    for item in json.loads((ctx.root/'audit/inventory.json').read_text())[1:]:
        with (ctx.root/'features'/item['dataset']/'samples.csv').open(newline='') as handle:
            for column in next(csv.reader(handle)):
                if column not in meta_columns:meta_columns.append(column)
    def unique(cols):return list(dict.fromkeys(cols))
    pred_columns=unique(meta_columns+['prediction']+dimensions+['adaptation_fraction','repeat','strategy'])
    adapt_pred_columns=unique(meta_columns+['prediction','baseline_prediction']+dimensions+['adaptation_fraction','repeat','strategy'])
    skip_columns=dimensions+['adaptation_fraction','repeat','reason']
    temp_root=Path(tempfile.mkdtemp(prefix='.collect-',dir=ctx.root/'external'))
    handles={};writers={};paths={};row_counts={}
    def target(key,name,columns):
        path=temp_root/name
        handle=path.open('w',newline='')
        writer=csv.DictWriter(handle,fieldnames=columns,extrasaction='raise')
        writer.writeheader();handles[key]=handle;writers[key]=writer;paths[key]=path;row_counts[key]=0
    for kind in ['combined','per_gse']:
        target((kind,'no'),f'external_{kind}_no_adapt.csv',no_columns)
        target((kind,'ada'),f'external_{kind}_adaptation.csv',ada_columns)
        target((kind,'np'),f'external_{kind}_no_adapt_predictions.csv',pred_columns)
        target((kind,'ap'),f'external_{kind}_adaptation_predictions.csv',adapt_pred_columns)
    target(('skip',),'skipped_gse.csv',skip_columns)
    def validate(marker):
        if not ctx.cache.valid(marker):
            raise RuntimeError('Invalid external checkpoint during collection: '+str(marker))
        return marker.parent
    def progress(force=False):
        nonlocal last
        now=time.time()
        if force or now-last>=RUNTIME['heartbeat_seconds']:
            state=dict(stage='external_collect',driver_pid=os.getpid(),updated=now,complete=count,
                       total=len(markers),failed=0,running_or_queued=0,workers=RUNTIME['io_workers'],
                       elapsed_seconds=now-start,rows={str(k):v for k,v in row_counts.items()})
            write_json(ctx.root/'reports/external_collect_progress.json',state)
            LOG.info('COLLECT external %d/%d checkpoints %.1fs',count,len(markers),now-start)
            last=now
    try:
        with ThreadPoolExecutor(max_workers=RUNTIME['io_workers']) as pool:
            pending=deque();todo=iter(markers)
            for _ in range(RUNTIME['io_workers']*2):
                try:pending.append(pool.submit(validate,next(todo)))
                except StopIteration:break
            while pending:
                future=pending.popleft()
                while not future.done():
                    wait([future],timeout=RUNTIME['heartbeat_seconds']);progress()
                path=future.result()
                try:pending.append(pool.submit(validate,next(todo)))
                except StopIteration:pass
                if (path/'skipped.json').exists():
                    writers[('skip',)].writerow(json.loads((path/'skipped.json').read_text()))
                    row_counts[('skip',)]+=1
                elif (path/'metrics.csv').exists():
                    with (path/'metrics.csv').open(newline='') as handle:
                        metric_rows=list(csv.DictReader(handle))
                    adapt=path.name!='no_adapt'
                    for row in metric_rows:
                        kind=row['source_kind'];key=(kind,'ada' if adapt else 'no')
                        writers[key].writerow(row);row_counts[key]+=1
                    if metric_rows:
                        key=(metric_rows[0]['source_kind'],'ap' if adapt else 'np')
                        with (path/'predictions.csv').open(newline='') as handle:
                            for row in csv.DictReader(handle):
                                writers[key].writerow(row);row_counts[key]+=1
                count+=1;progress()
        for handle in handles.values():handle.close()
        # Preserve the original combined table's baseline-then-adapted row ordering.
        combined=temp_root/'external_combined_adaptation_all.csv'
        with combined.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=all_columns)
            writer.writeheader()
            for key in [('combined','no'),('combined','ada')]:
                with paths[key].open(newline='') as source:
                    writer.writerows(csv.DictReader(source))
        for key,path in paths.items():
            if key==('combined','ada'):continue
            os.replace(path,ctx.root/'external'/path.name)
        os.replace(combined,ctx.root/'external'/combined.name)
        t=pd.read_csv(ctx.root/'external/external_combined_no_adapt.csv')
        if len(t):
            best=t.sort_values(['R2','MAE'],ascending=[False,True]).groupby(['K','model']).head(1).copy()
            best['selection_label']='exploratory_external_selected_fold'
            write_csv(ctx.root/'external/exploratory_external_selected_fold.csv',best)
        progress(True)
    finally:
        for handle in handles.values():
            if not handle.closed:handle.close()
        shutil.rmtree(temp_root,ignore_errors=True)




def bags(ds,pre,rows):
    key=(str(ds['path']),hashlib.sha256(pre.indices.tobytes()+pre.scaler.scale_.tobytes()).hexdigest())
    if key not in _BAG_CACHE:
        if len(_BAG_CACHE)>=2:_BAG_CACHE.pop(next(iter(_BAG_CACHE)))
        cells=ds['cells'][:,pre.indices].multiply(1/np.maximum(pre.scaler.scale_,1e-12)).tocsr();ix=ds['bag_index']
        _BAG_CACHE[key]=[(cells[np.flatnonzero(ix['sample']==i)],ix['celltype'][ix['sample']==i]) for i in range(len(ds['meta']))]
    return [_BAG_CACHE[key][int(i)] for i in rows]

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
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),initializer=worker_init,initargs=(ctx.cfg,RUNTIME,thread_override)) as pool:
        todo=iter(items)
        def submit():
            try:item=next(todo)
            except StopIteration:return False
            future=pool.submit(fn,ctx.cfg,item);active[future]=item;return True
        for _ in range(workers*2):
            if not submit():break
        while active:
            done,_=wait(active,timeout=RUNTIME['heartbeat_seconds'],return_when=FIRST_COMPLETED)
            for future in done:
                item=active.pop(future)
                try:future.result()
                except BaseException as e:
                    error=dict(item=item,error=str(e),traceback=traceback.format_exc());failures.append(error)
                    LOG.error('%s task failed %s: %s',kind,item,e)
                completed+=1
                if psutil.virtual_memory().available/2**30<RUNTIME['reserve_memory_gb']:
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


@contextmanager
def bootstrap_claim(path):
    """Serialize one bootstrap across the main pipeline and its cache precomputer."""
    import fcntl
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a+') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(handle,fcntl.LOCK_UN)


@lru_cache(maxsize=1)
def bootstrap_arrays(path):
    # Bootstrap needs no per-cell reservoir: retain exactly the same X/meta/schema.
    root=Path(path)
    return dict(path=root,X=sp.load_npz(root/'pseudobulk.npz'),
                meta=pd.read_csv(root/'samples.csv',keep_default_na=False),
                schema=json.loads((root/'schema.json').read_text()))


def bootstrap_one(cfg,item):
    key=json.dumps(cfg,sort_keys=True)
    if key not in _BOOTSTRAP_CONTEXT:_BOOTSTRAP_CONTEXT[key]=context(cfg)
    ctx=_BOOTSTRAP_CONTEXT[key]
    tag,rows,K,params,b=item
    path=ctx.root/'stability'/tag/'bootstrap'/f'boot_{b:04d}.joblib';marker=path.with_suffix('.complete.json')
    with bootstrap_claim(path.with_suffix('.lock')):
        if ctx.cache.valid(marker):return
        _bootstrap_fit(ctx,tag,rows,K,params,b,path,marker)


def _bootstrap_fit(ctx,tag,rows,K,params,b,path,marker):
    cfg=ctx.cfg;job_log(ctx,'bootstrap',f'{tag}_{b:04d}')
    ds=bootstrap_arrays(str(ctx.root/'features/train'));rows=np.asarray(rows,int)
    meta=ds['meta'].iloc[rows];groups=meta.group_id.unique();rng=np.random.default_rng(cfg['seed']+b);draw=rng.choice(groups,len(groups),replace=True)
    ix=np.concatenate([rows[np.flatnonzero(meta.group_id.to_numpy()==g)] for g in draw]);LOG.info('BOOTSTRAP %s %d started n=%d',tag,b,len(ix))
    rank=models.rank_features(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],cfg)
    pre=models.Preprocessor.fit(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],ds['schema']['features'],ds['meta'].sample_id.iloc[ix],K,cfg,rank)
    model=fit('ElasticNet',pre.transform(ds['X'][ix]),ds['meta'].Age.to_numpy(float)[ix],cfg,cfg['seed']+b,params)
    record=dict(indices=pre.indices,coefficient=np.asarray(model.coef_,np.float32),input_coefficient=np.asarray(model.coef_/pre.scaler.scale_,np.float32),bootstrap=b,draw_groups=draw,fit_ids=pre.fit_ids,params=params,solver='sklearn')
    dump(path,record);ctx.cache.mark(marker,[path]);LOG.info('BOOTSTRAP %s %d complete',tag,b)

def bootstrap_parallel(ctx,ds,rows,K,params,tag):
    tasks=[(tag,np.asarray(rows).tolist(),K,params,b) for b in range(ctx.cfg['stability']['n_bootstrap'])]
    run_jobs(ctx,'bootstrap_'+tag,bootstrap_one,tasks,RUNTIME['bootstrap_workers'])
    return _ORIGINAL['bootstrap'](ctx,ds,rows,K,params,tag)
