#!/usr/bin/env python
"""Evaluate every saved Model x K x fold with the existing FT30 -> all100 method."""
from pathlib import Path
import argparse,json,os,sys,subprocess,time,itertools,traceback
ROOT=Path(__file__).resolve().parent

def package_path(lineage):
    package=ROOT/'Fibro_Accelerated_20261102' if lineage=='Fibro' else ROOT
    sys.path[:0]=[str(package/'vendor'),str(package),str(ROOT)]
if os.environ.get('TABULASTROMAL_FT30_LINEAGE'):package_path(os.environ['TABULASTROMAL_FT30_LINEAGE'])

def all_items(cfg):
    return [(int(K),name,int(fold)) for fold,name,K in itertools.product(
        range(1,cfg['cv']['n_splits']+1),cfg['models']['names'],cfg['features']['K'])]


def execute_task(cfg,item,protocol):
    """Same fitted algorithm as ft30.task; correct metadata for non-best folds."""
    import ft30,pandas as pd
    from tabulastromal.storage import write_json,write_csv
    ctx=ft30.context(cfg);ft30.task(cfg,item,protocol)
    best={(K,name):fold for K,name,fold in ft30.selected_folds(ctx)}
    K,name,fold=item
    if best[K,name]==fold:return
    _,target,_=ft30.source_and_target(ctx,item);marker=target/'ft30.complete.json'
    state=json.loads(marker.read_text())
    r=pd.read_csv(target/'metrics_ft30.csv').iloc[0]
    if r.get('source_fold_selection')=='all_benchmark_folds':return
    outputs=[Path(p['path']) for p in state['outputs']]
    for p in outputs:
        if p.suffix=='.csv':
            t=pd.read_csv(p)
            if 'fold_best' in t:
                t['fold_best']=best[K,name];t['source_fold_selection']='all_benchmark_folds'
                t['is_internal_best_fold']=False;write_csv(p,t)
        elif p.name=='split_ft30.json':
            d=json.loads(p.read_text());d.update(evaluation_fold=fold,internal_best_fold=best[K,name],
                evaluation_fold_scope='all_benchmark_folds');write_json(p,d)
    ctx.cache.mark(marker,outputs,ft30_protocol=protocol,evaluation_fold_scope='all_benchmark_folds')


def collect(ctx,protocol,items,failures):
    import ft30,numpy as np,pandas as pd
    from tabulastromal.storage import write_csv,write_json,file_hash
    from tabulastromal.metrics import evaluate,METRICS
    mapping={'metrics_ft30.csv':'metrics','predictions_ft30.csv':'predictions',
        'baseline_metrics_ft30.csv':'baseline_metrics','heldout_metrics_ft30.csv':'heldout_metrics',
        'heldout_predictions_ft30.csv':'heldout_predictions','heldout_baseline_metrics_ft30.csv':'heldout_baseline_metrics'}
    tables={name:[] for name in mapping.values()};markers=[];status=[]
    best={(K,name):fold for K,name,fold in ft30.selected_folds(ctx)}
    for item in items:
        K,name,fold=item;_,target,_=ft30.source_and_target(ctx,item);marker=target/'ft30.complete.json'
        key=dict(lineage=ctx.lineage,K=K,model=name,fold=fold)
        if not ft30.valid(ctx,marker,protocol):
            status.append(dict(key,status='failed',adaptation_status='invalid_task_checkpoint'));continue
        for src,dst in mapping.items():
            t=pd.read_csv(target/src)
            t['fold_best']=best[K,name];t['is_internal_best_fold']=fold==best[K,name]
            t['evaluation_fold_scope']='all_benchmark_folds';tables[dst].append(t)
        r=tables['metrics'][-1].iloc[0]
        status.append(dict(key,status='complete',adaptation_status=r.adaptation_status,
            adaptation_strategy=r.adaptation_strategy,adaptation_failure_reason=r.adaptation_failure_reason,
            n_adapt=r.n_adapt,n_evaluation=r.n_evaluation,source_model_path=r.source_model_path))
        markers.append(marker)
    outputs=[]
    for name,parts in tables.items():
        table=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame()
        tables[name]=table;p=ctx.root/'external'/f'external_ft30_all_folds_{name}.csv'
        write_csv(p,table);outputs.append(p)
        if name in ['metrics','predictions','heldout_metrics']:
            for fold,part in table.groupby('fold'):
                p=ctx.root/'external'/f'external_ft30_fold_{int(fold)}_{name}.csv';write_csv(p,part);outputs.append(p)
    p=ctx.root/'external/external_ft30_all_folds_status.csv';write_csv(p,pd.DataFrame(status));outputs.append(p)
    if len(markers)!=len(items) or failures:raise RuntimeError('Incomplete all-fold tasks; inspect reports/ft30_all_folds_failures.json')
    meta=pd.read_csv(ctx.root/'external/external_sample_split_30_70.csv');allids=set(meta.sample_id)
    adapt=set(meta.loc[meta.set=='adaptation','sample_id']);hold=set(meta.loc[meta.set=='holdout','sample_id'])
    assert not adapt&hold
    for key,part in tables['predictions'].groupby(['K','model','fold']):
        assert set(part.sample_id)==allids and not part.sample_id.duplicated().any()
        K,name,fold=key
        rows=tables['metrics'];r=rows[(rows.K==K)&(rows.model==name)&(rows.fold==fold)].iloc[0]
        expected=evaluate(part.Age,part.prediction)
        for m in METRICS:np.testing.assert_allclose(expected[m],r[m],rtol=1e-8,atol=1e-8,equal_nan=True)
        assert set(part.loc[part.used_for_adaptation,'sample_id'])==(adapt if r.adaptation_status=='adapted' else set())
        _,target,_=ft30.source_and_target(ctx,key);ids=json.loads((target/'split_ft30.json').read_text())
        assert not set(ids['fit_ids'])&hold and not set(ids['preprocessor_fit_ids'])&allids
        assert not set(ids['source_refit_ids'])&set(ids['source_valid_ids'])
    for key,part in tables['heldout_predictions'].groupby(['K','model','fold']):
        assert set(part.sample_id)==hold and not part.sample_id.duplicated().any()
        K,name,fold=key;r=tables['heldout_metrics'];r=r[(r.K==K)&(r.model==name)&(r.fold==fold)].iloc[0]
        for m,value in evaluate(part.Age,part.prediction).items():
            if m in METRICS:np.testing.assert_allclose(value,r[m],rtol=1e-8,atol=1e-8,equal_nan=True)
    pins=json.loads((ctx.root/'reports/ft30_cv_reference.json').read_text());assert all(file_hash(ctx.root/p)==sha for p,sha in pins.items())
    cv=pd.read_csv(ctx.root/'benchmark/cv_fold_metrics.csv');v=cv[cv.split=='Valid']
    selected=v.sort_values(['R2','MAE','K'],ascending=[False,True,True],na_position='last').groupby(['fold','model'],sort=True).head(1).copy()
    selected['K_selection']='max same-fold INTERNAL Valid R2; min MAE; min K'
    p=ctx.root/'external/ft30_fold_internal_K_selection.csv';write_csv(p,selected);outputs.append(p)
    write_json(ctx.root/'reports/ft30_all_folds_validation.json',dict(status='passed',
        expected=len(items),completed=len(markers),n_external=len(allids),n_adaptation=len(adapt),n_heldout=len(hold),
        metrics_recomputed=True,train_valid_unchanged=True,external_not_used_for_K_selection=True))
    return outputs+markers


def run_lineage(args,options):
    lineage=args.lineage;os.environ['TABULASTROMAL_FT30_LINEAGE']=lineage;package_path(lineage)
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,options['gpu_ids'])) if options['gpu_enabled'] else ''
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/endo_matplotlib'))
    import yaml,ft30,logging,psutil,pandas as pd
    from tabulastromal import pipeline
    from tabulastromal.storage import process_lock,write_json,write_csv,file_hash,digest
    cfgpath=Path(args.config) if args.config else ROOT/options['lineages'][lineage]['config']
    cfg=yaml.safe_load(cfgpath.read_text());out=Path(cfg['output_root'])/lineage;items=all_items(cfg)
    prior=json.loads((out/'reports/ft30_driver.json').read_text());identity=json.loads((out/'reports/run_identity.json').read_text())
    record=json.loads((out/'external/ft30_protocol.json').read_text());protocol=record['ft30_protocol']
    if prior['source']['ft30.py']!=file_hash(ROOT/'ft30.py'):raise RuntimeError('Adaptation core changed; original FT30 protocol must be reconciled before reuse')
    assert record['options']['adaptation_fraction']==.3 and record['options']['adaptation_seed']==42
    assert all(record['options'][k]==options[k] for k in ['continuation_epochs','continuation_boost_rounds'])
    ft30.configure(cfg,lineage,options)
    with process_lock(out/'checkpoints/pipeline.lock'):
        ctx=pipeline.Context(cfg,lineage,True)
        assert ctx.cache.signature==identity['signature']
        existing=[i for i in items if ft30.valid(ctx,ft30.source_and_target(ctx,i)[1]/'ft30.complete.json',protocol)]
        pending=[i for i in items if i not in existing]
        if args.dry_run:
            for i in items:
                source,_,_=ft30.source_and_target(ctx,i)
                assert ctx.cache.valid(source/'complete.json'),str(source)
            print(json.dumps(dict(lineage=lineage,expected=len(items),reusable=len(existing),pending=len(pending),
                Ks=cfg['features']['K'],folds=cfg['cv']['n_splits'],output=str(out),source_bundles_valid=True)));return 0
        if args.smoke_test:
            candidate=next((i for i in pending if i[1]=='ElasticNet'),next(i for i in items if i[1]=='ElasticNet'))
            execute_task(cfg,candidate,protocol)
            target=ft30.source_and_target(ctx,candidate)[1]
            pred=pd.read_csv(target/'predictions_ft30.csv');split=pd.read_csv(out/'external/external_sample_split_30_70.csv')
            assert set(pred.sample_id)==set(split.sample_id) and not pred.sample_id.duplicated().any()
            assert pred.used_for_adaptation.sum()==sum(split.set=='adaptation')
            assert pd.read_csv(target/'metrics_ft30.csv').iloc[0].adaptation_status=='adapted'
            write_json(out/'reports/ft30_all_folds_smoke.json',dict(status='passed',item=candidate,
                n_external=len(pred),source_fold=candidate[2],task_protocol=protocol,non_best_fold_supported=True))
            return 0
        pipeline.logger_for(out);log=logging.getLogger('all_folds')
        manifest=[dict(lineage=lineage,K=K,model=name,fold=fold,source_model_path=str(ft30.source_and_target(ctx,(K,name,fold))[0]),
            reusable=(K,name,fold) in existing,task_protocol=protocol) for K,name,fold in items]
        write_csv(out/'external/ft30_all_folds_task_manifest.csv',pd.DataFrame(manifest))
        scope=digest(dict(task_protocol=protocol,items=items,runner=file_hash(__file__)))
        write_json(out/'reports/ft30_all_folds_driver.json',dict(pid=os.getpid(),lineage=lineage,task_protocol=protocol,
            scope_signature=scope,source_signature=ctx.cache.signature,runner_sha256=file_hash(__file__),
            adaptation_core_sha256=file_hash(ROOT/'ft30.py'),expected=len(items),reused=len(existing),pending=len(pending)))
        from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
        import multiprocessing
        failures=[];finished=len(existing);started=time.time();active={};todo=iter(pending);exhausted=False
        workers=min(options['workers_per_lineage'],max(1,len(pending)))
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
            initializer=ft30.worker_init,initargs=(cfg,lineage,options)) as pool:
            while active or not exhausted:
                while len(active)<workers and not exhausted:
                    if psutil.virtual_memory().available/2**30<options['reserve_memory_gb']:break
                    try:item=next(todo)
                    except StopIteration:exhausted=True;break
                    active[pool.submit(execute_task,cfg,item,protocol)]=item
                done=wait(active,timeout=options['heartbeat_seconds'],return_when=FIRST_COMPLETED)[0] if active else []
                for future in done:
                    item=active.pop(future)
                    try:future.result()
                    except Exception as exc:failures.append(dict(item=item,error=str(exc),traceback=traceback.format_exc()));log.exception('All-fold task failed %s',item)
                    finished+=1
                write_json(out/'reports/ft30_all_folds_progress.json',dict(expected=len(items),finished=finished,failed=len(failures),
                    reused=len(existing),active=len(active),workers=workers,pid=os.getpid(),updated=time.time(),elapsed=time.time()-started))
                log.info('ALL FOLDS %s %d/%d failed=%d active=%d reused=%d',lineage,finished,len(items),len(failures),len(active),len(existing))
                if not active and not exhausted:time.sleep(5)
        write_json(out/'reports/ft30_all_folds_failures.json',failures)
        outputs=collect(ctx,protocol,items,failures)
        marker=out/'checkpoints/stage_external_all_folds.json';ctx.cache.mark(marker,outputs,ft30_protocol=protocol,scope_signature=scope)
        write_json(out/'reports/ft30_all_folds_complete.json',dict(status='complete',expected=len(items),task_protocol=protocol,
            source_signature=ctx.cache.signature,scope_signature=scope,completed_at=time.time(),runner_sha256=file_hash(__file__)))
        return 0


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--lineage',choices=['both','Fibro','Endo'],default='both')
    p.add_argument('--stage',choices=['all','external','figures'],default='all');p.add_argument('--resume',action='store_true')
    p.add_argument('--workers',type=int);p.add_argument('--dry-run',action='store_true');p.add_argument('--config');p.add_argument('--smoke-test',action='store_true')
    p.add_argument('--ft30-config',default=str(ROOT/'ft30_config.yaml'));p.add_argument('--dpi',type=int,default=350)
    args=p.parse_args();import yaml
    options=yaml.safe_load(Path(args.ft30_config).read_text())
    if args.workers:options['workers_per_lineage']=args.workers
    if args.lineage=='both':
        if args.config:p.error('Single-lineage --config required for custom scientific config')
        if args.stage!='figures':
            children=[]
            for lineage in ['Fibro','Endo']:
                command=[sys.executable,'-u',str(Path(__file__).resolve()),'--lineage',lineage,'--stage','external',
                    '--resume','--ft30-config',str(Path(args.ft30_config).resolve())]
                if args.workers:command+=['--workers',str(args.workers)]
                if args.dry_run:command+=['--dry-run']
                if args.smoke_test:command+=['--smoke-test']
                with (ROOT/'logs'/f'{lineage}_ft30_all_folds.log').open('ab') as log:
                    children.append((lineage,subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)))
            codes={lineage:child.wait() for lineage,child in children};print(json.dumps(codes),flush=True)
            if any(codes.values()):return 2
        if args.stage in ['all','figures'] and not args.dry_run and not args.smoke_test:
            subprocess.run([sys.executable,'-u',str(ROOT/'Stromal_Combined_Figure/combine_stromal_figure.py'),
                '--all-folds','--dpi',str(args.dpi)],cwd=ROOT,check=True)
        return 0
    if args.stage=='figures':p.error('Use --lineage both for combined figures')
    return run_lineage(args,options)
if __name__=='__main__':raise SystemExit(main())
