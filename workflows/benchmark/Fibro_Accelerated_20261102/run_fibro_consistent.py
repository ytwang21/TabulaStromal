#!/usr/bin/env python
"""Complete resumable Fibro pipeline, Endo-compatible external task artifacts."""
from pathlib import Path
import argparse
import json
import os
import sys
ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'vendor'), str(ROOT)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(ROOT/'config_fibro_fast.yaml'))
    parser.add_argument('--stage', choices=['all','audit','features','benchmark','external','stability','enrichment','figures'], default='all')
    parser.add_argument('--resume', action='store_true', help='Validated completed work is always reused')
    parser.add_argument('--workers', type=int, default=64)
    parser.add_argument('--bootstrap-workers', type=int, default=64)
    parser.add_argument('--worker-threads', type=int, default=1)
    parser.add_argument('--reserve-memory-gb', type=float, default=192)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--prepare-only', action='store_true', help='Validate/migrate existing inputs and caches, then exit')
    parser.add_argument('--no-import', action='store_true', help='Skip old source migration; completed fast checkpoints still reused')
    args = parser.parse_args(argv)
    if min(args.workers,args.bootstrap_workers,args.worker_threads)<1 or args.reserve_memory_gb<0:
        parser.error('Worker counts/threads must be positive, reserve memory nonnegative')
    import yaml
    cfg = yaml.safe_load(Path(args.config).read_text())
    if set(cfg['lineages']) != {'Fibro'}:
        raise ValueError('Only Fibro configuration is permitted')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key] = str(args.worker_threads)
    os.environ['CUDA_VISIBLE_DEVICES'] = ','.join(map(str,cfg['accelerated']['gpu_ids'])) if cfg['accelerated']['gpu_enabled'] else ''
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/matplotlib'))
    os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'runtime/numba'))
    os.environ['PYTHONPATH'] = os.pathsep.join([str(ROOT/'vendor'),str(ROOT),os.environ.get('PYTHONPATH','')])
    options = dict(workers=args.workers,bootstrap_workers=args.bootstrap_workers,
                   worker_threads=args.worker_threads,reserve_memory_gb=args.reserve_memory_gb)
    if args.dry_run:
        print(json.dumps(dict(config=str(Path(args.config).resolve()),output=cfg['output_root'],
              scientific_config_unchanged=True,runtime_options=options,
              phase_order=['combined_no_adapt','per_gse_no_adapt','combined_adaptation','per_gse_adaptation'],
              task_layout='tasks/{dataset}/K{K}/{model}/fold{fold}/{no_adapt|adapt_fX.X_rR}',
              resume='SHA256-validated existing model/feature/external checkpoints',
              Endo='Existing process and source untouched'),indent=2))
        return 0
    from tabulastromal import accelerator, pipeline
    from tabulastromal.storage import process_lock,write_json
    import fibro_external_consistent
    accelerator.install(cfg)
    fibro_external_consistent.install(options)
    root = Path(cfg['output_root'])
    with process_lock(root/'driver.lock'):
        ctx = pipeline.Context(cfg,'Fibro',True)
        pipeline.logger_for(ctx.root)
        if not args.no_import:
            accelerator.migrate(ctx)
        if args.prepare_only:
            return 0
        write_json(ctx.root/'reports/driver.json',dict(pid=os.getpid(),entry=str(Path(__file__).resolve()),
                   config=cfg,signature=ctx.cache.signature,runtime_options=options))
        write_json(ROOT/'reports/consistent_driver.json',dict(pid=os.getpid(),entry=str(Path(__file__).resolve()),
                   config_path=str(Path(args.config).resolve()),output=str(ctx.root),signature=ctx.cache.signature,
                   runtime_options=options,Endo_untouched=True))
        (ROOT/'logs').mkdir(exist_ok=True)
        (ROOT/'logs/Fibro_fast.pid').write_text(str(os.getpid())+'\n')
        try:
            pipeline.run_lineage(cfg,'Fibro',args.stage,True)
        except KeyboardInterrupt:
            return 130
        except Exception:
            import logging
            logging.exception('Fibro pipeline incomplete; successful task checkpoints retained')
            return 2
    return 0


if __name__=='__main__':
    raise SystemExit(main())
