#!/usr/bin/env python
"""Raw h5ad -> 14-model grouped CV -> external -> stable clocks -> enrichment/figures."""
import argparse,os,sys,json,copy,subprocess
from pathlib import Path
import yaml


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lineage',choices=['Fibro','Endo','both','Fibroblast','Endothelial'],default='both')
    stages=['audit','features','benchmark','external','stability','enrichment','figures']
    parser.add_argument('--stage',choices=stages+['all'],default='all')
    parser.add_argument('--config',default=str(Path(__file__).with_name('config.yaml')))
    parser.add_argument('--resume',action='store_true',help='Only reuse current input/config/source signature and verified artifacts')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--seed',type=int)
    parser.add_argument('--threads',type=int)
    parser.add_argument('--workers',type=int,help='Independent K/fold tasks per lineage')
    parser.add_argument('--bootstrap-workers',type=int)
    parser.add_argument('--gpu',help='Physical GPU IDs, e.g. 0,3,4; cpu disables CUDA; independent per-lineage defaults')
    parser.add_argument('--stop-after',choices=stages)
    args=parser.parse_args(argv);cfg=yaml.safe_load(Path(args.config).read_text())
    args.lineage={'Fibroblast':'Fibro','Endothelial':'Endo'}.get(args.lineage,args.lineage)
    if args.seed is not None:cfg['seed']=args.seed
    for key,val in [('threads',args.threads),('benchmark_workers',args.workers),('bootstrap_workers',args.bootstrap_workers)]:
        if val is not None:
            if val<1:parser.error(key+' must be positive')
            cfg['resources'][key]=val
    expected=['XGBoost','ExtraTrees','ElasticNet','MLP','Ridge','Lasso','RandomForest','GradientBoosting','Supervised_PCA_Ridge','SVR','IdealAge(DeepSets)','IdealAge(Pseudobulk)','IdealAge(Ensemble)','JungClock']
    if cfg['models']['names'] not in [expected,expected[:-1]]:raise ValueError('Use the canonical thirteen self-contained models or append optional JungClock as the fourteenth')
    active=cfg['models']['names']
    if cfg['cv']['n_splits']<2:raise ValueError('At least 2 folds required')
    if len(set(cfg['features']['K']))!=len(cfg['features']['K']):raise ValueError('Duplicate K')
    if cfg['stability']['n_bootstrap']<1:raise ValueError('Invalid bootstrap count')
    lineages=['Fibro','Endo'] if args.lineage=='both' else [args.lineage]
    if args.dry_run:
        report=dict(config=str(Path(args.config).resolve()),stage=args.stage,revision=cfg.get('execution_revision'),resources=cfg['resources'],lineages={})
        for lin in lineages:
            p=cfg['lineages'][lin];files=list(Path(p['external_gse_dir']).rglob('*.h5ad'))
            report['lineages'][lin]=dict(output=p.get('output_dir',str(Path(cfg['output_root'])/lin)),
                inputs={k:dict(path=v,exists=Path(v).exists()) for k,v in p.items() if k!='output_dir'},
                models=active,K=cfg['features']['K'],folds=cfg['cv']['n_splits'],benchmark_tasks=len(active)*len(cfg['features']['K'])*cfg['cv']['n_splits'],
                per_gse_files=len(files),external_scope=cfg['resources'].get('external_scope','all'),per_gse_external_deferred=cfg['resources'].get('external_scope')=='combined',adaptation_upper_bound=(1 if cfg['resources'].get('external_scope')=='combined' else 1+len(files))*len(active)*len(cfg['features']['K'])*cfg['cv']['n_splits']*len(cfg['external']['fractions'])*cfg['external']['repeats'],
                bootstrap_fits=(cfg['cv']['n_splits']+1)*cfg['stability']['n_bootstrap'])
        print(json.dumps(report,indent=2));return 0
    if args.lineage=='both':
        original=list(argv if argv is not None else sys.argv[1:]);procs=[]
        # Separate driver PIDs; each uses its own output tree, lock, GPU and log.
        for lin in lineages:
            child=original.copy()
            if '--lineage' in child:child[child.index('--lineage')+1]=lin
            else:child+=['--lineage',lin]
            procs.append(subprocess.Popen([sys.executable,'-u',str(Path(__file__).resolve())]+child))
        codes=[p.wait() for p in procs];return 0 if all(c==0 for c in codes) else 2
    threads=str(cfg['resources']['threads'])
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]=threads
    gpu=args.gpu or cfg['resources'].get('lineage_gpu',{}).get(args.lineage,'cpu')
    if cfg['resources']['device']=='cpu':gpu='cpu'
    os.environ['CUDA_VISIBLE_DEVICES']='' if gpu=='cpu' else str(gpu)
    if cfg['resources']['device']=='auto' or cfg['resources'].get('xgb_device')=='auto':
        import torch
        available=(gpu!='cpu' and torch.cuda.is_available())
        cfg['resources']['device']='cuda' if available else 'cpu'
        cfg['resources']['xgb_device']='cuda' if available else 'cpu'
        print(f'{args.lineage}: physical GPU={gpu}; CUDA available={available}; resolved device={cfg["resources"]["device"]}',flush=True)
    cfg['resources']['physical_gpu']=gpu
    for k,v in [('MPLCONFIGDIR','matplotlib'),('NUMBA_CACHE_DIR','numba')]:os.environ.setdefault(k,str(Path(cfg['output_root'])/'runtime'/v))
    vendor=cfg.get('jungclock',{}).get('python_vendor')
    if vendor:sys.path.insert(0,vendor)
    from tabulastromal.pipeline import run_lineage
    try:run_lineage(cfg,args.lineage,args.stage,args.resume,args.stop_after)
    except KeyboardInterrupt:return 130
    except Exception as e:
        import traceback
        traceback.print_exc();print(f'FAILED {args.lineage}: {e}',file=sys.stderr);return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
