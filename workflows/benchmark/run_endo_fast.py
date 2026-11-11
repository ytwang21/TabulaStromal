#!/usr/bin/env python
"""Endothelial execution adapter: original seven-stage pipeline, original outputs."""
from pathlib import Path
import argparse,json,os,sys
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT)]


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lineage',choices=['Endo'],default='Endo')
    p.add_argument('--config',default=str(ROOT/'config.yaml'))
    p.add_argument('--runtime-config',default=str(ROOT/'endo_runtime.yaml'))
    p.add_argument('--stage',choices=['all','audit','features','benchmark','external','stability','enrichment','figures'],default='all')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--workers',type=int)
    p.add_argument('--bootstrap-workers',type=int)
    p.add_argument('--worker-threads',type=int)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--stop-after',choices=['audit','features','benchmark','external','stability','enrichment','figures'])
    args=p.parse_args(argv)
    import yaml
    cfg=yaml.safe_load(Path(args.config).read_text())
    runtime=yaml.safe_load(Path(args.runtime_config).read_text())
    for key in ['workers','bootstrap_workers','worker_threads']:
        value=getattr(args,key)
        if value is not None:runtime[key]=value
    if min(runtime[k] for k in ['workers','bootstrap_workers','worker_threads','io_workers'])<1:
        p.error('worker/thread counts must be positive')
    if args.dry_run:
        print(json.dumps(dict(lineage='Endo',output=str(Path(cfg['output_root'])/'Endo'),
            stages=['audit','features','benchmark','external','stability','enrichment','figures'],
            scientific_configuration_unchanged=True,runtime=runtime,
            resume='Existing matching SHA256 checkpoints; no migration or retraining of completed source models'),indent=2))
        return 0
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key]=str(runtime['worker_threads'])
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,runtime['gpu_ids'])) if runtime['gpu_enabled'] else ''
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/endo_matplotlib'))
    os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'runtime/endo_numba'))
    os.environ['PYTHONPATH']=os.pathsep.join([str(ROOT),os.environ.get('PYTHONPATH','')])
    from tabulastromal import pipeline
    from tabulastromal.storage import write_json,process_lock,file_hash
    import endo_acceleration
    endo_acceleration.install(cfg,runtime)
    root=Path(cfg['output_root'])/'Endo'
    # The same pipeline lock excludes both the legacy runner and another fast runner.
    # The driver lock protects preflight/provenance before original run_lineage takes over.
    with process_lock(root/'checkpoints/fast_driver.lock'):
        from contextlib import ExitStack
        with ExitStack() as check:
            check.enter_context(process_lock(root/'checkpoints/pipeline.lock'))
            identity_path=root/'reports/run_identity.json'
            previous=json.loads(identity_path.read_text()) if identity_path.exists() else None
            ctx=pipeline.Context(cfg,'Endo',True)
            if previous and ctx.cache.signature!=previous['signature']:
                write_json(identity_path,previous)
                raise RuntimeError('Scientific source/config/input signature changed; refusing to relabel existing checkpoints')
            write_json(root/'reports/endo_acceleration_identity.json',dict(pid=os.getpid(),runtime=runtime,
                preserved_signature=ctx.cache.signature,config_path=str(Path(args.config).resolve()),
                executor={name:file_hash(ROOT/name) for name in ['run_endo_fast.py','endo_acceleration.py','endo_external_parallel.py']},
                policy='Original stage graph, task schema and saved estimator checkpoints; runtime-only acceleration',
                solver_note='Original sklearn estimators/objective/tolerance/iteration budgets preserved; unfinished EN/Lasso use float64 and fixed-seed random coordinate traversal; numerical/coefficient differences remain explicit'))
        try:
            pipeline.run_lineage(cfg,'Endo',args.stage,True,args.stop_after)
        except KeyboardInterrupt:return 130
        except Exception:
            import logging
            logging.exception('Endo interrupted/incomplete; completed checkpoints retained')
            return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
