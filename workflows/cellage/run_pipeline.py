#!/usr/bin/env python3
"""Single reproducible entry point for TabulaStromal CellAge Exploration."""
from __future__ import annotations
import argparse
from datetime import datetime
import hashlib
import importlib
import json
import logging
import os
from pathlib import Path
import sys
import signal
from zoneinfo import ZoneInfo

PHASES = {'audit': 'step00_audit', 'merge': 'step01_merge_scores',
          'subtypes': 'step02_celltype_senescence', 'age': 'step03_age_associations',
          'concordance': 'step04_method_concordance', 'clock': 'step05_genemodel_integration',
          'genes': 'step06_gene_pathway_analysis', 'figure': 'step07_figure'}


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', default=str(Path(__file__).with_name('config.yaml')))
    p.add_argument('--lineage', choices=['fibro', 'endo', 'both'], default='both')
    p.add_argument('--phase', choices=[*PHASES, 'all', 'smoke-test'], default='all')
    p.add_argument('--resume', action='store_true', help='Continue compatible latest run; completed stages remain immutable')
    p.add_argument('--run-id', help='Explicit immutable run ID (required to select an older run)')
    p.add_argument('--source-run', help='For --phase figure only: reuse audited numerical outputs in a NEW run; preserve original figures')
    p.add_argument('--checkpoints', action='store_true', help='Read-only validation of stage receipts; do not run analyses')
    p.add_argument('--dry-run', action='store_true', help='Validate configuration and paths; read no numeric tables; write nothing')
    p.add_argument('--output-root')
    p.add_argument('--n-jobs', type=int, default=1, help='BLAS thread cap; lineages execute sequentially to bound RAM')
    p.add_argument('--max-memory-gb', type=float)
    return p.parse_args(argv)


def run_pipeline(argv=None):
    args = arguments(argv)
    if args.source_run and args.phase != 'figure':
        raise ValueError('--source-run is only supported for --phase figure')
    if args.n_jobs < 1:
        raise ValueError('--n-jobs must be positive')
    for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[k] = str(args.n_jobs)
    import yaml
    from scripts.common import environment_record, fingerprint, memory_guard, read_json, write_json
    from scripts.checkpointing import stage_valid, stage_work, seal_stage, migrate_run
    cfg = yaml.safe_load(Path(args.config).read_text())
    if args.output_root:
        cfg['output_root'] = str(Path(args.output_root).resolve())
    if args.max_memory_gb:
        cfg['max_memory_gb'] = args.max_memory_gb
    cfg['max_memory_gb'] = float(cfg.get('max_memory_gb', 16))
    lineages = ['fibro', 'endo'] if args.lineage == 'both' else [args.lineage]
    root = Path(cfg['output_root']).resolve()
    os.environ.setdefault('MPLCONFIGDIR', str(root/'.runtime/matplotlib'))
    if args.dry_run:
        plan = {'lineages': lineages, 'phase': args.phase, 'output_root': str(root), 'n_jobs': args.n_jobs,
                'donor_identity_evidence': cfg.get('donor_identity_evidence'), 'inputs': []}
        for lin in lineages:
            c = cfg['lineages'][lin]
            for role, path in [(k, c[k]) for k in ['atlas', 'age_atlas', 'raw_atlas']] + [(k, s['path']) for k, s in c['scores'].items()]:
                p = Path(path)
                plan['inputs'].append({'lineage': lin, 'role': role, 'path': path, 'exists': p.is_file(), 'bytes': p.stat().st_size if p.is_file() else None})
        print(json.dumps(plan, indent=2, ensure_ascii=False))
        return None
    if args.phase == 'smoke-test':
        from scripts.smoke_test import run_smoke
        return run_smoke(cfg, root)
    root.mkdir(parents=True, exist_ok=True)
    config_hash = fingerprint(cfg)
    inputs = []
    for conf in cfg['lineages'].values():
        paths = [conf[k] for k in ['atlas', 'age_atlas', 'raw_atlas']]
        paths += [s['path'] for s in conf['scores'].values()]
        paths += [conf['genemodel'][k] for k in ['predictions', 'stable_features', 'feature_universe', 'gene_set_config']]
        for name in paths:
            p = Path(name)
            inputs.append((str(p), p.stat().st_size if p.exists() else None, p.stat().st_mtime_ns if p.exists() else None))
    inputs_hash = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    # A code digest invalidates checkpoints when implementation changes.
    code_hash = hashlib.sha256(b''.join(p.read_bytes() for p in sorted(Path(__file__).parent.glob('scripts/*.py'))) + Path(__file__).read_bytes()).hexdigest()
    run_id = args.run_id
    migration_source = None
    latest = read_json(root/'latest_run.json') if (root/'latest_run.json').exists() else {}
    legacy_migration = False
    compatibility_path = Path(__file__).with_name('checkpoint_upgrade_legacy.json')
    compatibility = read_json(compatibility_path) if compatibility_path.exists() else {}
    if args.resume and not run_id:
        if latest.get('config_hash') == config_hash and latest.get('code_hash') == code_hash:
            if latest.get('inputs_hash') == inputs_hash:
                run_id = latest['run_id']
            else:
                migration_source = Path(latest['path'])
        known_upgrade = compatibility.get('target_code_hash') == code_hash
        if not run_id and migration_source is None and known_upgrade and latest.get('config_hash') == config_hash:
            if latest.get('code_hash') in compatibility.get('compatible_checkpoint_code_hashes',[]):
                migration_source = Path(latest['path'])
        if not run_id and migration_source is None and known_upgrade:
            source = Path(compatibility['source_run'])
            previous = read_json(source/'run_identity.json')
            if previous.get('config_hash') == config_hash and previous.get('code_hash') == compatibility['legacy_code_hash']:
                migration_source,legacy_migration = source,True
    if args.checkpoints:
        if not run_id:
            candidate = root/'runs'/args.run_id if args.run_id else (Path(latest['path']) if (root/'latest_run.json').exists() else None)
        else:
            candidate = root/'runs'/run_id
        if candidate is None or not candidate.is_dir():
            raise ValueError('No run to inspect; specify --run-id or use --resume')
        for lin in lineages:
            for phase in PHASES:
                who = 'both' if phase == 'figure' else lin
                scope = lineages if phase == 'figure' else None
                print(f'{lin}/{phase}: '+('validated' if stage_valid(candidate,phase,who,scope) else 'not_validated_or_legacy'))
        print('RUN_DIR='+str(candidate))
        return candidate
    run_id = run_id or datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y%m%d_%H%M%S_%f')
    if Path(run_id).name != run_id or run_id in {'.', '..'}:
        raise ValueError('run-id must be a single directory name')
    run = root/'runs'/run_id
    identity = {'run_id': run_id, 'config_hash': config_hash, 'code_hash': code_hash, 'inputs_hash': inputs_hash}
    if run.exists():
        old = read_json(run/'run_identity.json')
        if not args.resume or old != identity:
            raise ValueError('Existing run is immutable or incompatible; use a new run ID')
    else:
        for sub in ['audit', 'tables', 'cache', 'figures/panels', 'figures/source_data', 'results', 'logs', 'checkpoints']:
            (run/sub).mkdir(parents=True, exist_ok=True)
        write_json(run/'run_identity.json', identity)
        write_json(run/'audit/input_signature.json', inputs)
        (run/'config.resolved.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))
        write_json(run/'logs/environment.json', environment_record())
        (run/'logs/requirements-lock.txt').write_text('\n'.join(f'{k}=={v}' for k, v in environment_record()['packages'].items()) + '\n')
        # Preserve source code used by this run.
        import shutil
        snapshot = run/'logs/code_snapshot'
        shutil.copytree(Path(__file__).parent/'scripts', snapshot/'scripts', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(__file__, snapshot/'run_pipeline.py')
    # Acquire the run lock BEFORE migrating/copying/checking any outputs.
    import fcntl
    lock_handle = open(run/'.pipeline.lock', 'a')
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise RuntimeError(f'Another process is already using {run}') from exc
    if migration_source:
        if legacy_migration:
            for name,expected in compatibility['modules'].items():
                archived=migration_source/'logs/code_snapshot/scripts'/name
                if not archived.exists() or hashlib.sha256(archived.read_bytes()).hexdigest()!=expected:
                    raise ValueError('Legacy code snapshot does not match checkpoint-upgrade allowlist')
        migrate_run(cfg,migration_source,run,inputs,PHASES,legacy=legacy_migration)
    if args.source_run:
        source = Path(args.source_run).resolve()
        previous = read_json(source/'run_identity.json')
        if previous['config_hash'] != config_hash or previous['inputs_hash'] != inputs_hash:
            raise ValueError('Source run configuration or selected input signature differs; numerical reuse rejected')
        # Figure-only revisions may reuse numerical results only if every
        # numerical analysis module is byte-identical to the original snapshot.
        for p in (Path(__file__).parent/'scripts').glob('*.py'):
            if p.name == 'step07_figure.py':
                continue
            old = source/'logs/code_snapshot/scripts'/p.name
            if not old.exists() or old.read_bytes() != p.read_bytes():
                compatibility = read_json(Path(__file__).with_name('checkpoint_upgrade_legacy.json'))
                legacy_matches = all(
                    (source/'logs/code_snapshot/scripts'/name).exists() and
                    hashlib.sha256((source/'logs/code_snapshot/scripts'/name).read_bytes()).hexdigest() == expected
                    for name,expected in compatibility['modules'].items() if name != 'step07_figure.py')
                if not legacy_matches or compatibility.get('target_code_hash') != code_hash:
                    raise ValueError(f'Numerical module changed: {p.name}; run the analysis again')
                break
        import shutil
        for sub in ['audit', 'tables', 'cache']:
            for original in (source/sub).rglob('*'):
                if original.is_file():
                    target = run/sub/original.relative_to(source/sub)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if not target.exists():
                        # Independent copies: a future resumed summary writer
                        # must never be able to modify the original run inode.
                        shutil.copy2(original, target)
        for original in (source/'checkpoints').glob('*.json'):
            if original.name != 'both_figure.json':
                shutil.copy2(original, run/'checkpoints'/original.name)
        if (source/'checkpoints/receipts').exists():
            shutil.copytree(source/'checkpoints/receipts',run/'checkpoints/receipts',dirs_exist_ok=True)
        write_json(run/'reused_results_provenance.json', {'source_run':str(source),
                   'source_identity':previous, 'verified_numerical_modules_unchanged':True,
                   'policy':'Existing numerical artifacts are immutable; figure outputs and reports regenerated in new run'})
    write_json(root/'latest_run.json', {**identity, 'path': str(run)})
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.StreamHandler(), logging.FileHandler(run/'logs/pipeline.log')], force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    log = logging.getLogger('cellage')
    log.info('RUN_DIR=%s', run)
    def interrupted(signum, _frame):
        log.warning('Received signal %s; committed work units remain reusable',signum)
        raise SystemExit(128+signum)
    old_term = signal.signal(signal.SIGTERM,interrupted)
    phases = list(PHASES) if args.phase == 'all' else [args.phase]
    if 'audit' in phases:
        from scripts.step00_audit import inventory
        if not (run/'audit/input_manifest.csv').exists():
            inventory(cfg, run)
    # Sequential execution bounds peak memory; threads controlled above.
    for phase in phases:
        executed_stage = False
        memory_guard(cfg)
        if phase == 'figure':
            module = importlib.import_module('scripts.' + PHASES[phase])
            marker = run/'checkpoints/both_figure.json'
            if args.resume and stage_valid(run,phase,'both',lineages):
                log.info('Resume skips completed figure')
            else:
                with stage_work(run,phase,'both'):
                    module.run(cfg, run, lineages, source_run=args.source_run)
                    seal_stage(run,phase,'both',lineages)
            continue
        for lin in lineages:
            marker = run/'checkpoints'/f'{lin}_{phase}.json'
            if args.resume and stage_valid(run,phase,lin):
                log.info('Resume skips %s/%s (%s)', lin, phase, read_json(marker)['status'])
                continue
            if phase != 'audit' and not (run/'checkpoints'/f'{lin}_audit.json').exists():
                raise RuntimeError('Run --phase audit first in the same compatible run, or use --phase all')
            module = importlib.import_module('scripts.' + PHASES[phase])
            with stage_work(run,phase,lin):
                module.run(cfg, run, lin)
                seal_stage(run,phase,lin)
                executed_stage = True
        if phase == 'audit':
            from scripts.step00_audit import finalize
            if executed_stage or not (run/'audit/audit_summary.md').exists():
                finalize(cfg, run, lineages)
    log.info('Finished requested phases. Conditional/blocked analyses are recorded in checkpoints/.')
    print('RUN_DIR=' + str(run))
    signal.signal(signal.SIGTERM,old_term)
    return run


if __name__ == '__main__':
    run_pipeline()
