#!/usr/bin/env python
"""Precompute original bootstrap caches alongside the formal lineage pipeline.

No stage is marked complete here. The original pipeline retains responsibility
for aggregation, stable-feature decisions and final fits.
Per-bootstrap claims prevent duplicate fits when the main pipeline catches up.
"""
from pathlib import Path
import argparse,json,os,sys,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=str(ROOT/'config.yaml'))
    p.add_argument('--runtime-config',default=str(ROOT/'endo_runtime.yaml'))
    p.add_argument('--workers',type=int,default=32)
    p.add_argument('--outer-only',action='store_true',help='Legacy mode: precompute only outer-fold caches')
    args=p.parse_args()
    if args.workers<1:p.error('workers must be positive')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key]='1'
    import yaml,numpy as np,pandas as pd
    from tabulastromal import pipeline,training
    from tabulastromal.storage import process_lock,write_json,file_hash
    import endo_acceleration as fast
    cfg=yaml.safe_load(Path(args.config).read_text())
    options=yaml.safe_load(Path(args.runtime_config).read_text())
    options.update(workers=args.workers,bootstrap_workers=args.workers,worker_threads=1,gpu_enabled=False)
    out=Path(cfg['output_root'])/'Endo'
    # Separate cache-only launcher lock; never take the main pipeline lock.
    with process_lock(out/'checkpoints/bootstrap_precompute.lock'):
        previous=json.loads((out/'reports/run_identity.json').read_text())['signature']
        fast.install(cfg,options);ctx=fast.context(cfg)
        if ctx.cache.signature!=previous:raise RuntimeError('Scientific signature mismatch')
        if not ctx.cache.valid(out/'checkpoints/stage_benchmark.json'):
            raise RuntimeError('Requires completed current-signature benchmark')
        ds=fast.bootstrap_arrays(str(out/'features/train'))
        ranking=pd.read_csv(out/'benchmark/model_comparison.csv')
        K=int(ranking.loc[ranking.model=='ElasticNet'].iloc[0].K)
        specs=[]
        for fold,(tr,va) in enumerate(training.grouped_splits(ds['meta'],cfg['cv']['n_splits'],cfg['seed'],cfg['cv']['age_bins']),1):
            source=training.task_path(ctx,K,'ElasticNet',fold)
            if not ctx.cache.valid(source/'complete.json'):raise RuntimeError(str(source))
            params=json.loads((source/'model.json').read_text())['params']
            specs.append((f'outer{fold}',tr.tolist(),K,params))
        (out/'logs/precompute/logs').mkdir(parents=True,exist_ok=True)
        pipeline.logger_for(out/'logs/precompute')
        # Work from the opposite end to complement the formal pipeline's queue.
        items=[(*spec,b) for spec in reversed(specs) for b in reversed(range(cfg['stability']['n_bootstrap']))]
        report=out/'reports/bootstrap_precompute_execution.json'
        state=dict(status='running',pid=os.getpid(),started=time.time(),workers=args.workers,
            K=K,total=len(items)+(0 if args.outer_only else cfg['stability']['n_bootstrap']),signature=ctx.cache.signature,
            helper_sha256=file_hash(__file__),adapter_sha256=file_hash(ROOT/'endo_acceleration.py'),
            policy='Original training-only EN tuning, group resampling, MI, scaling, parameters and fit; reverse bootstrap queues; wait for formal full-train tuning rather than duplicate it; per-unit flock; no stage completion or external selection')
        write_json(report,state)
        try:
            fast.run_jobs(ctx,'bootstrap_precompute_outer',fast.bootstrap_one,items,args.workers)
            if not args.outer_only:
                path=out/'checkpoints/tuning/full_train'/f'K{K}_ElasticNet.json'
                marker=path.with_suffix('.complete.json')
                def main_running():
                    import fcntl
                    p=out/'checkpoints/pipeline.lock'
                    if not p.exists():return False
                    with p.open('r') as handle:
                        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
                        except BlockingIOError:return True
                        fcntl.flock(handle,fcntl.LOCK_UN);return False
                while not ctx.cache.valid(marker) and main_running():
                    state.update(status='waiting_for_main_full_train_tuning',updated=time.time())
                    write_json(report,state)
                    __import__('logging').getLogger(__name__).info('Waiting for formal full_train tuning; avoiding duplicate MI/candidate fits')
                    time.sleep(options['heartbeat_seconds'])
                rows=np.arange(len(ds['meta']))
                params=json.loads(path.read_text())['params'] if ctx.cache.valid(marker) else fast.tune_linear(ctx,ds,rows,K,'ElasticNet','full_train')
                state.update(status='running_full_train',updated=time.time());write_json(report,state)
                full=[('full_train',rows.tolist(),K,params,b) for b in reversed(range(cfg['stability']['n_bootstrap']))]
                fast.run_jobs(ctx,'bootstrap_precompute_full_train',fast.bootstrap_one,full,args.workers)
        except BaseException as exc:
            state.update(status='failed',finished=time.time(),error=str(exc));write_json(report,state);raise
        state.update(status='complete',finished=time.time());write_json(report,state)
    return 0

if __name__=='__main__':raise SystemExit(main())
