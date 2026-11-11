#!/usr/bin/env python
"""CLI entry point. Both lineages run sequentially with conservative thread limits."""
import argparse,os,sys
from pathlib import Path
import yaml


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lineage',choices=['Fibro','Endo','both'],default='both')
    parser.add_argument('--stage',choices=['audit','features','benchmark','external','stability','enrichment','figures','all'],default='all')
    parser.add_argument('--config',default=str(Path(__file__).with_name('config.yaml')))
    parser.add_argument('--resume',action='store_true',help='Reuse only complete matching checkpoints')
    parser.add_argument('--dry-run',action='store_true',help='Inspect configuration, input existence and work counts; do not load expression or train')
    parser.add_argument('--seed',type=int)
    parser.add_argument('--stop-after',choices=['audit','features','benchmark','external','stability','enrichment','figures'])
    args=parser.parse_args(argv);cfg=yaml.safe_load(Path(args.config).read_text())
    if args.seed is not None:cfg['seed']=args.seed
    expected=['XGBoost','ExtraTrees','ElasticNet','MLP','Ridge','Lasso','RandomForest','GradientBoosting','Supervised_PCA_Ridge','SVR','DeepSets','Pseudobulk','Ensemble']
    if set(cfg['models']['names'])!=set(expected):raise ValueError('All 13 named models are required; use smoke configuration only for smaller K/folds/iterations')
    if cfg['cv']['n_splits']<2:raise ValueError('At least 2 folds required')
    if len(set(cfg['features']['K']))!=len(cfg['features']['K']):raise ValueError('Duplicate feature K values')
    if cfg['stability']['n_bootstrap']<1:raise ValueError('n_bootstrap must be positive')
    threads=str(cfg['resources']['threads'])
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]=threads
    os.environ.setdefault('MPLCONFIGDIR',str(Path(cfg['output_root'])/'runtime/matplotlib'))
    os.environ.setdefault('NUMBA_CACHE_DIR',str(Path(cfg['output_root'])/'runtime/numba'))
    lineages=['Fibro','Endo'] if args.lineage=='both' else [args.lineage]
    if args.dry_run:
        import json
        report=dict(config=str(Path(args.config).resolve()),output_root=cfg['output_root'],stage=args.stage,seed=cfg['seed'],sequential_lineages=True,threads=int(threads),lineages={})
        for lin in lineages:
            p=cfg['lineages'][lin];files=sorted(Path(p['external_gse_dir']).rglob('*.h5ad'))
            report['lineages'][lin]=dict(inputs={k:dict(path=v,exists=Path(v).exists()) for k,v in p.items()},per_gse_files=len(files),
                models=cfg['models']['names'],K=cfg['features']['K'],folds=cfg['cv']['n_splits'],benchmark_tasks=len(expected)*len(cfg['features']['K'])*cfg['cv']['n_splits'],
                adaptation_upper_bound=(1+len(files))*len(expected)*len(cfg['features']['K'])*cfg['cv']['n_splits']*len(cfg['external']['fractions'])*cfg['external']['repeats'],
                bootstrap_fits=(cfg['cv']['n_splits']+1)*cfg['stability']['n_bootstrap'])
        print(json.dumps(report,indent=2));return 0
    from tabulastromal.pipeline import run_lineage
    failures=[]
    for lin in lineages:
        try:run_lineage(cfg,lin,args.stage,args.resume,args.stop_after)
        except KeyboardInterrupt:return 130
        except Exception as e:
            failures.append((lin,str(e)));print(f'FAILED {lin}: {e}',file=sys.stderr)
    if failures:
        print('Incomplete requested analysis: '+repr(failures),file=sys.stderr);return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
