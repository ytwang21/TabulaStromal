"""External evaluation with Endo's task contract and shared scientific functions.

Only execution order, process scheduling, and progress exports differ. Kept outside
tabulastromal/ so existing, scientifically identical completed fits keep their
validated signature. This executor's hashes are recorded separately in reports.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import itertools
import json
import logging
import multiprocessing
import os
import time
import traceback
from functools import lru_cache
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tabulastromal import accelerator as fast, external, models, pipeline, training
from tabulastromal.metrics import evaluate, METRICS
from tabulastromal.storage import write_json, write_csv, file_hash

LOG = logging.getLogger(__name__)
OPTIONS = {}
BASE_RUN_JOBS = fast.run_jobs


@lru_cache(maxsize=3)
def arrays(path, with_cells=False):
    if with_cells:
        return fast.dataset(path)
    root = Path(path)
    # The serial external functions use only X/meta/schema for regression and
    # calibration. Loading training cell reservoirs here would waste RAM and I/O.
    return dict(path=root, X=sp.load_npz(root/'pseudobulk.npz'),
                meta=pd.read_csv(root/'samples.csv', keep_default_na=False),
                schema=json.loads((root/'schema.json').read_text()))


def worker_init(cfg, threads):
    fast.worker_init(cfg, threads)


def run_pool(ctx, phase, fn, items, workers, threads, reserve_gb):
    """Bounded spawn pool, memory admission control, preserved completed units."""
    import psutil
    tasks = iter(items)
    total = len(items)
    completed, failures, active = 0, [], {}
    start = time.time()
    if not total:
        return
    with ProcessPoolExecutor(max_workers=min(workers, total),
                             mp_context=multiprocessing.get_context('spawn'),
                             initializer=worker_init, initargs=(ctx.cfg, threads)) as pool:
        exhausted = False
        def submit():
            nonlocal exhausted
            if exhausted:
                return False
            if psutil.virtual_memory().available / 2**30 < reserve_gb:
                return False
            try:
                item = next(tasks)
            except StopIteration:
                exhausted = True
                return False
            active[pool.submit(fn, ctx.cfg, item)] = item
            return True
        while active or not exhausted:
            while len(active) < min(workers, total) and submit():
                pass
            done = wait(active, timeout=ctx.cfg['accelerated']['heartbeat_seconds'],
                        return_when=FIRST_COMPLETED)[0] if active else []
            for future in done:
                item = active.pop(future)
                try:
                    future.result()
                except Exception as exc:
                    failures.append(dict(item=item, error=str(exc), traceback=traceback.format_exc()))
                    LOG.exception('External task failed: %s', item)
                completed += 1
            state = dict(stage='external', phase=phase, driver_pid=os.getpid(),
                         updated=time.time(), complete=completed, total=total,
                         failed=len(failures), running_or_queued=len(active),
                         workers=min(workers,total), worker_threads=threads,
                         elapsed_seconds=time.time()-start,
                         unit='dataset/K/model/fold within this phase, including cache hits')
            write_json(ctx.root/'reports/external_progress.json', state)
            write_json(ctx.root/f'reports/external_{phase}_progress.json', state)
            LOG.info('PROGRESS external/%s %d/%d failed=%d active=%d',
                     phase, completed, total, len(failures), len(active))
            if not active and not exhausted:
                # No queued unit is lost under low memory; admission resumes later.
                LOG.warning('Insufficient available RAM; preserving %.1f GiB reserve', reserve_gb)
                time.sleep(min(10, ctx.cfg['accelerated']['heartbeat_seconds']))
    write_json(ctx.root/f'reports/external_{phase}_failures.json', failures)
    if failures:
        raise RuntimeError(f'{phase}: {len(failures)} task failures; successful units remain resumable')


def base_info(cfg, item, meta, inventory):
    dsid, K, name, fold, mode = item
    kind = 'combined' if dsid == 'combined' else 'per_gse'
    return dict(lineage='Fibro', K=int(K), model=name, fold=int(fold),
                GSE='__combined__' if kind=='combined' else str(meta.GSE.iloc[0]),
                dataset=dsid, method=inventory['method'], source_kind=kind,
                primary=inventory['method'] in ['stored_log1p', cfg['external']['primary_method']])


def external_unit(cfg, item):
    dsid, K, name, fold, mode = item
    K, fold = int(K), int(fold)
    ctx = fast.context(cfg)
    fast.job_log(ctx, 'external_consistent', f'{dsid}_K{K}_{name}_fold{fold}_{mode}')
    root = ctx.root/'external/tasks'/dsid/f'K{K}'/name/f'fold{fold}'
    inventory = next(r for r in json.loads((ctx.root/'audit/inventory.json').read_text()) if r['dataset']==dsid)
    # An already complete no-adaptation task must not load arrays or fit an adapter.
    nop = root/'no_adapt'
    if mode == 'no_adapt' and ctx.cache.valid(nop/'complete.json'):
        return
    meta = pd.read_csv(ctx.root/'features'/dsid/'samples.csv', keep_default_na=False)
    if not len(meta):
        raise ValueError('Empty dataset should have been excluded at dataset level: '+dsid)
    info_base = base_info(cfg, item, meta, inventory)
    outputs = [nop/'complete.json']
    pending = []
    for fraction, rep in itertools.product(cfg['external']['fractions'], range(cfg['external']['repeats'])):
        ap = root/f'adapt_f{fraction:.1f}_r{rep}'
        marker = ap/'complete.json'
        if mode == 'no_adapt':
            break
        if ctx.cache.valid(marker):
            outputs.append(marker)
            continue
        too_small = dsid!='combined' and meta.group_id.nunique() < cfg['external']['min_gse_samples']
        split, reason = external.adaptation_split(meta, fraction, cfg['seed']+rep,
                                                cfg['external']['min_adapt_groups'], cfg['external']['min_test_groups'])
        if too_small or split is None:
            # Same canonical per-fraction skip files and fields as Endo.
            write_json(ap/'skipped.json', dict(**info_base, adaptation_fraction=fraction, repeat=rep,
                       reason='GSE has fewer than five independent groups' if too_small else reason))
            ctx.cache.mark(marker, [ap/'skipped.json'], outcome='skipped')
            outputs.append(marker)
        else:
            pending.append((fraction, rep, split))
    if mode=='adapt' and not pending and ctx.cache.valid(nop/'complete.json'):
        ctx.cache.mark(root/'group_complete.json', outputs,
                       executor='Endo-compatible task layout; existing completed fits retained')
        return
    source = training.task_path(ctx, K, name, fold)
    if not ctx.cache.valid(source/'complete.json'):
        raise RuntimeError('Incomplete benchmark source: '+str(source))
    ext = arrays(str(ctx.root/'features'/dsid), name in ['DeepSets','Ensemble'])
    _, model, pre, info = models.load_bundle(source, cfg)
    y = ext['meta'].Age.to_numpy(float)
    if ctx.cache.valid(nop/'complete.json'):
        prediction = pd.read_csv(nop/'predictions.csv').prediction.to_numpy()
    else:
        LOG.info('EXTERNAL %s K=%d %s fold=%d n=%d', dsid, K, name, fold, len(meta))
        bags = fast.bags(ext, pre, np.arange(len(meta))) if name in ['DeepSets','Ensemble'] else None
        prediction = models.predict_model(model, name, pre.transform(ext['X']), bags)
        record = dict(**info_base, adaptation_fraction=0., repeat=0, strategy='no_adaptation',
                      n_adapt=0, n_test=len(meta), **evaluate(y, prediction))
        write_csv(nop/'metrics.csv', pd.DataFrame([record]))
        table = meta.copy()
        table['prediction'] = prediction
        for key, value in info_base.items():
            table[key] = value
        table['adaptation_fraction'], table['repeat'], table['strategy'] = 0., 0, 'no_adaptation'
        write_csv(nop/'predictions.csv', table)
        ctx.cache.mark(nop/'complete.json', [nop/'metrics.csv', nop/'predictions.csv'])
    if mode == 'no_adapt':
        return
    train = arrays(str(ctx.root/'features/train'))
    tr = np.flatnonzero(train['meta'].sample_id.isin(info['train_ids']))
    for fraction, rep, (ai, ti) in pending:
        ap = root/f'adapt_f{fraction:.1f}_r{rep}'
        assert set(meta.group_id.iloc[ai]).isdisjoint(meta.group_id.iloc[ti])
        LOG.info('ADAPT %s K%d %s fold%d f%.1f repeat%d n_adapt%d n_test%d',
                 dsid, K, name, fold, fraction, rep, len(ai), len(ti))
        write_json(ap/'split.json', dict(adapt_ids=meta.sample_id.iloc[ai].tolist(), test_ids=meta.sample_id.iloc[ti].tolist(),
                   adapt_groups=meta.group_id.iloc[ai].unique().tolist(), test_groups=meta.group_id.iloc[ti].unique().tolist(),
                   n_groups_adapt=meta.group_id.iloc[ai].nunique(), n_groups_test=meta.group_id.iloc[ti].nunique(),
                   actual_fraction=meta.group_id.iloc[ai].nunique()/meta.group_id.nunique(), seed=cfg['seed']+rep,
                   adapt_age=meta.Age.iloc[ai].describe().to_dict(), test_age=meta.Age.iloc[ti].describe().to_dict(),
                   adapt_GSE=meta.GSE.iloc[ai].value_counts().to_dict(), test_GSE=meta.GSE.iloc[ti].value_counts().to_dict()))
        # These are the same unmodified scientific routines used by Endothelial.
        # The accelerated linear solver hooks retain the pre-existing convex objective.
        adapter, strategy, detail = external.adapt_model(name, model, pre, train, ext, tr, ai, ctx, info, ap)
        after = external.predict_adapted(adapter, strategy, model, name, pre, ext, ti, ctx)
        before_metrics, after_metrics = evaluate(y[ti], prediction[ti]), evaluate(y[ti], after)
        record = dict(**info_base, adaptation_fraction=fraction, repeat=rep, strategy=strategy,
                      strategy_detail=detail, n_adapt=len(ai), n_test=len(ti),
                      small_sample_note='held-out n<20: unstable; repeat SD is descriptive' if len(ti)<20 else '',
                      **after_metrics)
        record.update({'baseline_'+m: before_metrics[m] for m in METRICS})
        record.update({'delta_'+m: after_metrics[m]-before_metrics[m] for m in METRICS})
        write_csv(ap/'metrics.csv', pd.DataFrame([record]))
        table = meta.iloc[ti].copy()
        table['prediction'], table['baseline_prediction'] = after, prediction[ti]
        for key, value in info_base.items():
            table[key] = value
        table['adaptation_fraction'], table['repeat'], table['strategy'] = fraction, rep, strategy
        write_csv(ap/'predictions.csv', table)
        marker = ap/'complete.json'
        ctx.cache.mark(marker, [ap/'metrics.csv', ap/'predictions.csv', ap/'split.json']+list(ap.glob('adapter*')))
        outputs.append(marker)
        del adapter
    ctx.cache.mark(root/'group_complete.json', outputs,
                   executor='Endo-compatible task layout; shared adaptation/prediction functions')


def collect(ctx):
    # The original collector produces the same file names/fields as Endo.
    external.collect_external(ctx)


def external_consistent(ctx):
    opts = OPTIONS
    inventory = json.loads((ctx.root/'audit/inventory.json').read_text())[1:]
    datasets, excluded = [], []
    for it in inventory:
        meta = pd.read_csv(ctx.root/'features'/it['dataset']/'samples.csv', keep_default_na=False)
        if not len(meta):
            path = ctx.root/'external/tasks'/it['dataset']/'excluded.json'
            write_json(path, dict(reason='no nonoverlapping eligible samples'))
            excluded.append(dict(dataset=it['dataset'], reason='no nonoverlapping eligible samples'))
        else:
            datasets.append(it['dataset'])
    def tasks(ids, mode):
        return [(dsid, K, name, fold, mode) for dsid, K, name, fold in
                itertools.product(ids, ctx.cfg['features']['K'], ctx.cfg['models']['names'],
                                  range(1, ctx.cfg['cv']['n_splits']+1))]
    combined = ['combined'] if 'combined' in datasets else []
    others = [dsid for dsid in datasets if dsid!='combined']
    phases = [('combined_no_adapt', tasks(combined, 'no_adapt')),
              ('per_gse_no_adapt', tasks(others, 'no_adapt')),
              ('combined_adaptation', tasks(combined, 'adapt')),
              ('per_gse_adaptation', tasks(others, 'adapt'))]
    contract = dict(signature=ctx.cache.signature, executor_sha256=file_hash(__file__),
                    entry_sha256=file_hash(Path(__file__).parent/'run_fibro_consistent.py'),
                    task_pattern='external/tasks/{dataset}/K{K}/{model}/fold{fold}/{no_adapt|adapt_fX.X_rR}',
                    dataset_order=datasets, excluded=excluded,
                    phase_order=[name for name,_ in phases], runtime_options=opts,
                    scientific_functions='tabulastromal.external.adaptation_split/adapt_model/predict_adapted/collect_external',
                    reuse_policy='Existing current-signature outputs are SHA256 validated; source training models are not refitted.',
                    solver_note='Existing Celer acceleration retained for EN/Lasso and bulk components. Same convex objective/tolerance checks; not bitwise identical to sklearn Endo.',
                    skip_policy='Canonical skipped.json + complete.json per skipped fraction/repeat, matching Endo')
    write_json(ctx.root/'reports/external_execution_contract.json', contract)
    try:
        for phase, items in phases:
            write_json(ctx.root/'reports/external_phase_status.json', dict(phase=phase,status='running',updated=time.time(),groups=len(items)))
            run_pool(ctx, phase, external_unit, items, opts['workers'], opts['worker_threads'], opts['reserve_memory_gb'])
            collect(ctx)
            if phase.endswith('no_adapt'):
                kind = 'combined' if phase.startswith('combined') else 'per_gse'
                frame = pd.read_csv(ctx.root/f'external/external_{kind}_no_adapt.csv')
                summary = frame.groupby(['lineage','K','model'])[METRICS].agg(['mean','std','count']) if len(frame) else pd.DataFrame()
                if len(summary):
                    summary.columns = ['_'.join(col) for col in summary.columns]
                    summary = summary.reset_index()
                write_csv(ctx.root/f'external/external_{kind}_no_adapt_summary.csv', summary)
            write_json(ctx.root/'reports/external_phase_status.json', dict(phase=phase,status='complete',updated=time.time(),groups=len(items)))
    except BaseException:
        # Canonical checkpoints always survive. Partial aggregate snapshots are not stage completion.
        try:
            collect(ctx)
        except Exception:
            LOG.exception('Partial external collection failed; individual checkpoints remain intact')
        raise
    write_csv(ctx.root/'external/failed_combinations.csv', pd.DataFrame(columns=[
        'lineage','K','model','fold','GSE','dataset','method','source_kind','primary','error','traceback']))
    # pipeline.run_lineage adds every canonical child complete.json itself, exactly
    # as in Endo. Do not duplicate the same checks through group_complete.json.
    return [ctx.root/'external'/f for f in ['external_combined_no_adapt.csv','external_per_gse_no_adapt.csv',
            'external_combined_adaptation_all.csv','external_per_gse_adaptation.csv','skipped_gse.csv']]


def install(options):
    global OPTIONS
    OPTIONS = dict(options)
    def stages(ctx, name):
        return external_consistent(ctx) if name=='external' else fast.stage(ctx, name)
    pipeline.run_stage = stages
    def runtime_jobs(ctx, kind, fn, items, workers):
        if kind.startswith('bootstrap'):
            workers = OPTIONS['bootstrap_workers']
        return BASE_RUN_JOBS(ctx, kind, fn, items, workers)
    fast.run_jobs = runtime_jobs
