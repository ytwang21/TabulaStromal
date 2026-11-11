#!/usr/bin/env python3
"""Command-line entry point; notebook imports the same stage functions."""
import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MKL_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'logs/matplotlib_cache'))
import signature_pipeline as sp


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.yaml'))
    p.add_argument('--lineage',choices=['fibro','endo','both'],default='both')
    p.add_argument('--stage',choices=['audit','markers','score','validate','figure','all'],default='all')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--seed',type=int)
    args=p.parse_args(); cfg=sp.read_config(args.config)
    if args.seed is not None: cfg['seed']=args.seed
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(sp.root(cfg)/'logs/pipeline.log')])
    lineages=['fibro','endo'] if args.lineage=='both' else [args.lineage]
    if args.dry_run:
        for lin in lineages:
            logging.info('%s input: %s',lin,sp.fingerprint(cfg['inputs'][lin]))
        logging.info('DRY RUN: stage=%s output=%s; no expression processing',args.stage,cfg['output_dir']);return
    provenance={'config':cfg,'inputs':{lin:sp.fingerprint(path) for lin,path in cfg['inputs'].items()},
                'notebooks':{lin:sp.fingerprint(path) for lin,path in cfg['notebooks'].items()},
                'prompt_sha256':hashlib.sha256(Path(cfg['prompt_file']).read_bytes()).hexdigest()}
    digest=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    manifest=sp.root(cfg)/'logs/config_fingerprint.json'
    if args.resume and manifest.exists() and sp.read_json(manifest)['digest']!=digest:
        raise ValueError('Config/input/seed changed: refuse stale --resume; rerun all without --resume or choose another output_dir.')
    sp.dump(manifest,{'digest':digest,'provenance':provenance})
    stages=['audit','markers','score','validate','figure'] if args.stage=='all' else [args.stage]
    for stage in stages:
        marker=sp.root(cfg)/f'logs/completed_{stage}_{args.lineage}.json'
        if args.resume and marker.exists():
            logging.info('Resume: stage %s already complete',stage);continue
        if stage=='audit': sp.stage_audit(cfg,lineages)
        elif stage=='markers': sp.marker_qc(cfg,lineages)
        elif stage=='score': sp.stage_score(cfg,lineages,resume=args.resume)
        elif stage=='validate': sp.stage_validate(cfg,lineages)
        elif stage=='figure': sp.stage_figure(cfg,lineages)
        sp.dump(marker,{'stage':stage,'lineages':lineages,'version':sp.VERSION,'seed':cfg['seed']})
    logging.info('Requested stages completed.')


if __name__=='__main__': main()
