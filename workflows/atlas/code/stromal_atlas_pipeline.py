#!/usr/bin/env python
"""StromalAtlas: read-only atlas audit, rendering, streamed signatures and state exploration."""
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import platform
import signal
import sys
import time
import traceback

# Set before numerical imports. User-provided thread settings take precedence.
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ.setdefault(name,'4')
os.environ.setdefault('MPLBACKEND','Agg')

HERE = Path(__file__).resolve().parent
for name, subdir in [('NUMBA_CACHE_DIR','numba_cache'),('MPLCONFIGDIR','matplotlib_cache')]:
    cache=HERE.parent/'tmp'/subdir
    cache.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault(name,str(cache))
STAGES = ['audit','atlas','markers','trajectory','figure']


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=HERE/'config.json')
    p.add_argument('--run-dir',type=Path,help='New directory or existing run to resume; default creates a unique timestamped run.')
    p.add_argument('--stage',choices=['all']+STAGES,default='all')
    p.add_argument('--lineage',choices=['both','fibro','endo'],default='both')
    p.add_argument('--resume',action='store_true',help='Skip completed stages only when input/config fingerprint matches.')
    p.add_argument('--reuse-atlas-run',type=Path,help='Copy validated audit/atlas/marker artifacts from a previous run; trajectory always uses the current configuration.')
    p.add_argument('--dpi',type=int,help='Override rendering DPI, useful for smoke tests')
    p.add_argument('--trajectory-cells',type=int)
    return p.parse_args(argv)


def main(argv=None):
    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'Received signal {signum}; incomplete stage will be rerun on resume')
    signal.signal(signal.SIGTERM, interrupted)
    args=parse_args(argv)
    from utils_data import audit, marker_signature, dump_json
    cfg=json.loads(args.config.read_text())
    for k in ['dpi','trajectory_cells']:
        if getattr(args,k) is not None: cfg[k]=getattr(args,k)
    for k,v in cfg['inputs'].items():
        cfg['inputs'][k]=str((args.config.resolve().parent/v).resolve())
    for k,v in cfg.get('metacell_inputs',{}).items():
        cfg['metacell_inputs'][k]=str((args.config.resolve().parent/v).resolve())
    if cfg['chunk_rows'] < 1 or cfg['trajectory_cells'] < 50:
        raise ValueError('chunk_rows >=1 and trajectory_cells >=50 required')
    if min(cfg['trajectory_neighbors']) < 2 or max(cfg['trajectory_neighbors']) >= cfg['trajectory_cells']:
        raise ValueError('Invalid trajectory neighbor counts')
    timestamp=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
    root=(args.run_dir or HERE.parent/'results'/'runs'/timestamp).resolve()
    if root.exists() and not args.resume:
        raise FileExistsError(f'{root} exists. Use --resume for this run, or choose a new run directory.')
    root.mkdir(parents=True,exist_ok=True)
    # O_EXCL prevents simultaneous writers to the same run; stale locks are documented.
    lock=root/'.running.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout),logging.FileHandler(root/'pipeline.log')],force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    try:
        inputs={k:dict(path=v,size=Path(v).stat().st_size,mtime_ns=Path(v).stat().st_mtime_ns) for k,v in cfg['inputs'].items()}
        for k,v in cfg.get('metacell_inputs',{}).items():
            inputs[f'metacell_{k}']=dict(path=v,size=Path(v).stat().st_size,mtime_ns=Path(v).stat().st_mtime_ns)
        fingerprint=hashlib.sha256(json.dumps(dict(config=cfg,inputs=inputs),sort_keys=True).encode()).hexdigest()
        manifest_path=root/'run_manifest.json'
        if manifest_path.exists():
            manifest=json.loads(manifest_path.read_text())
            if manifest['fingerprint'] != fingerprint: raise ValueError('Resume config/input mismatch. Start a new run.')
        else:
            manifest=dict(created_at=timestamp,timezone='Asia/Shanghai',fingerprint=fingerprint,
                          config=cfg,inputs=inputs,stages={})
        if args.reuse_atlas_run and not manifest['stages']:
            import shutil
            base=args.reuse_atlas_run.resolve();old=json.loads((base/'run_manifest.json').read_text())
            for lin in ['fibro','endo']:
                if old['inputs'][lin]!=inputs[lin]:raise ValueError('Reusable atlas input fingerprint mismatch')
                target=root/lin;target.mkdir(exist_ok=True)
                for source in (base/lin).iterdir():
                    if source.is_file() and not source.name.startswith('trajectory_'):
                        shutil.copy2(source,target/source.name)
                for stage in ['audit','atlas','markers']:
                    previous=old['stages'].get(f'{stage}:{lin}',{})
                    if previous.get('status')!='complete':raise ValueError('Source atlas stage incomplete')
                    manifest['stages'][f'{stage}:{lin}']=dict(status='complete',mode='reused',source_run=str(base),seconds=0)
            (root/'shared').mkdir(exist_ok=True)
            for source in (base/'shared').glob('data_audit*'):shutil.copy2(source,root/'shared'/source.name)
            manifest['reused_atlas_run']=str(base)
        dump_json(manifest_path,manifest);dump_json(root/'parameters.json',cfg)
        (root/'shared').mkdir(exist_ok=True)
        versions={}
        for name in ['numpy','pandas','scipy','h5py','scanpy','anndata','matplotlib','scikit-learn','igraph','umap-learn','pynndescent','numba','nbformat']:
            try: versions[name]=importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError: versions[name]='NOT_INSTALLED'
        dump_json(root/'shared'/'environment.json',dict(python=sys.executable,version=sys.version,platform=platform.platform(),packages=versions,
                  threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']}))
        (root/'shared'/'requirements-lock.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items() if v!='NOT_INSTALLED')+'\n')
        import shutil
        source_dir=root/'shared'/'source_snapshots'/timestamp
        source_dir.mkdir(parents=True,exist_ok=True)
        for source in HERE.iterdir():
            if source.is_file() and source.suffix in ['.py','.R','.json','.md','.txt','.sh']:
                shutil.copy2(source,source_dir/source.name)
        logging.info('RUN_DIR=%s',root)
        lineages=['fibro','endo'] if args.lineage=='both' else [args.lineage]
        stages=STAGES if args.stage=='all' else [args.stage]
        for stage in stages:
            if stage=='figure':
                tasks=['shared']
            else: tasks=lineages
            for lineage in tasks:
                key=f'{stage}:{lineage}'
                if args.resume and manifest['stages'].get(key,{}).get('status')=='complete':
                    logging.info('Already complete: %s',key);continue
                # Check concrete prerequisites before starting an expensive phase.
                required={'audit':[], 'atlas':['audit'], 'markers':['audit'],
                          'trajectory':['markers'], 'figure':[]}[stage]
                for dep in required:
                    if manifest['stages'].get(f'{dep}:{lineage}',{}).get('status')!='complete':
                        raise ValueError(f'{key} requires completed {dep}:{lineage}')
                if stage=='figure':
                    for lin in ['fibro','endo']:
                        for dep in ['atlas','markers','trajectory']:
                            if manifest['stages'].get(f'{dep}:{lin}',{}).get('status')!='complete':
                                raise ValueError(f'Integrated figure requires {dep}:{lin}')
                started=time.time()
                manifest['stages'][key]=dict(status='running',started_at=dt.datetime.now(dt.timezone.utc).isoformat())
                dump_json(manifest_path,manifest);logging.info('START %s',key)
                try:
                    if stage=='audit': audit(cfg['inputs'][lineage],lineage,cfg,root)
                    elif stage=='atlas':
                        from utils_plot import atlas_stage
                        atlas_stage(lineage,cfg,root)
                    elif stage=='markers':
                        marker_signature(cfg['inputs'][lineage],lineage,cfg,root)
                        from utils_plot import marker_stage
                        marker_stage(lineage,cfg,root)
                    elif stage=='trajectory':
                        if cfg.get('trajectory_engine')=='metacell_three_method':
                            from metacell_trajectory import trajectory_stage
                        else:
                            from archive.v1_utils_trajectory import trajectory_stage
                        trajectory_stage(cfg['inputs'][lineage],lineage,cfg,root)
                    else:
                        if cfg.get('trajectory_engine')!='metacell_three_method':
                            from archive.v1_utils_trajectory import refresh_reports
                            refresh_reports(root,cfg)
                        from utils_plot import integrated_figure
                        integrated_figure(cfg,root)
                    manifest['stages'][key].update(status='complete',seconds=round(time.time()-started,2))
                except BaseException as e:
                    manifest['stages'][key].update(status='failed',error=str(e),traceback=traceback.format_exc())
                    raise
                finally:
                    dump_json(manifest_path,manifest)
                    with (root/'phase_progress.md').open('a') as fp:
                        status=manifest['stages'][key]
                        fp.write(f'\n- {key}: {status["status"]}; output: {root / lineage}; elapsed: {time.time()-started:.1f}s; issue: {status.get("error","none")}. Next: resume remaining phases.\n')
            if stage=='audit':
                lines=[]
                for lin in ['fibro','endo']:
                    p=root/'shared'/f'data_audit_{lin}.json'
                    if p.exists():
                        d=json.loads(p.read_text());lines.append(f'{lin}: {d["n_cells"]:,} cells; {d["n_genes"]:,} genes; {d["n_tissues"]} tissues; {d["n_samples"]} samples; donor identifiers: {d["n_donors"]}; raw: {d["has_raw"]}; layers: {d["layers"]}; fields: {d["fields"]}')
                (root/'shared'/'data_audit_summary.txt').write_text('\n'.join(lines)+'\n')
        # Full hashes of code, not hundreds of GB of biological inputs.
        dump_json(root/'shared'/'code_sha256.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.is_file() and p.suffix in ['.py','.R','.json','.sh']})
        from report_run import report
        report(root)
        files=[dict(path=str(p.relative_to(root)),bytes=p.stat().st_size) for p in root.rglob('*') if p.is_file() and p.name not in ['file_manifest.csv','.running.lock']]
        import pandas as pd
        pd.DataFrame(files).to_csv(root/'file_manifest.csv',index=False)
        if manifest['stages'].get('figure:shared',{}).get('status')=='complete':
            from validate_run import validate
            result=validate(root)
            logging.info('Output validation: %s (%s checks)',result['status'],result['n_checks'])
        logging.info('SUCCESS requested phases finished. Results: %s',root)
        return root
    finally:
        lock.unlink(missing_ok=True)


if __name__=='__main__':
    main()
