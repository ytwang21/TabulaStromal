#!/usr/bin/env python
"""Complete isolated Fibro pipeline: parallel CPU, bounded GPU, validated old checkpoints."""
from pathlib import Path
import argparse,json,os,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
sys.path.insert(0,str(ROOT))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=str(ROOT/'config_fibro_fast.yaml'))
    parser.add_argument('--stage',choices=['audit','features','benchmark','external','stability','enrichment','figures','all'],default='all')
    parser.add_argument('--resume',action='store_true',help='All saved complete tasks are checked and reused')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--prepare-only',action='store_true',help='Validate/migrate input caches and old completed models, without training')
    parser.add_argument('--no-import',action='store_true',help='Fresh synthetic/smoke run only; do not import old artifacts')
    args=parser.parse_args(argv)
    import yaml
    cfg=yaml.safe_load(Path(args.config).read_text());a=cfg['accelerated']
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]=str(a['worker_threads'])
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,a['gpu_ids'])) if a['gpu_enabled'] else ''
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/matplotlib'));os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'runtime/numba'))
    os.environ['PYTHONPATH']=os.pathsep.join([str(ROOT/'vendor'),str(ROOT),os.environ.get('PYTHONPATH','')])
    if set(cfg['lineages'])!={'Fibro'}:raise ValueError('This entry point only accepts Fibro; no Endo writes are allowed')
    if args.dry_run:
        print(json.dumps(dict(output=cfg['output_root'],models=cfg['models']['names'],K=cfg['features']['K'],folds=cfg['cv']['n_splits'],main_model_combinations=len(cfg['models']['names'])*len(cfg['features']['K'])*cfg['cv']['n_splits'],hardware=a,
            bootstrap_per_lineage=(cfg['cv']['n_splits']+1)*cfg['stability']['n_bootstrap'],adaptive_repeats=cfg['external']['repeats'],fractions=cfg['external']['fractions'],Endo='original process and source files untouched'),indent=2));return 0
    from tabulastromal import pipeline
    from tabulastromal.storage import process_lock,write_json
    from tabulastromal import accelerator
    accelerator.install(cfg);pipeline.run_stage=accelerator.stage
    root=Path(cfg['output_root']);root.mkdir(parents=True,exist_ok=True)
    with process_lock(root/'driver.lock'):
        ctx=pipeline.Context(cfg,'Fibro',True);pipeline.logger_for(ctx.root)
        if not args.no_import:accelerator.migrate(ctx)
        write_json(ctx.root/'reports/driver.json',dict(pid=os.getpid(),entry=str(Path(__file__).resolve()),config=cfg,signature=ctx.cache.signature))
        if args.prepare_only:return 0
        try:pipeline.run_lineage(cfg,'Fibro',args.stage,True)
        except KeyboardInterrupt:return 130
        except Exception:
            import logging;logging.exception('Accelerated Fibro incomplete; see preserved task checkpoints');return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
