"""Crash-safe output commits, validated stage receipts and reusable work units.

A marker alone is never sufficient to skip work. SIGKILL/power loss may leave
*.partial files, which are never treated as committed outputs.
"""
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import tempfile
import uuid

LOG = logging.getLogger('cellage')
ACTIVE = ContextVar('checkpoint_store', default=None)
TERMINAL = {'complete', 'conditional', 'blocked'}


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


@contextmanager
def atomic_path(path):
    """Commit on the same filesystem only after the writer closes successfully."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.'+path.stem+'.', suffix=path.suffix+'.partial', dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        yield tmp
        with open(tmp, 'rb') as handle:
            os.fsync(handle.fileno())
        os.replace(tmp, path)
        fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    except BaseException:
        # An incomplete temporary file is harmless; keep it for diagnosis.
        raise


def _default(value):
    if hasattr(value, 'item'):
        return value.item()
    if value.__class__.__name__ in {'NAType', 'NaTType'}:
        return None
    return str(value)


def atomic_json(path, obj):
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=_default))


def artifacts(paths, root):
    root = Path(root)
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            files.extend(x for x in p.rglob('*') if x.is_file() and not x.name.endswith('.partial'))
        elif p.is_file():
            files.append(p)
    return [{'path': str(p.relative_to(root)), 'size': p.stat().st_size, 'sha256': digest(p)}
            for p in sorted(set(files))]


def valid_artifacts(entries, root):
    root = Path(root).resolve()
    try:
        for item in entries:
            p = (root/item['path']).resolve()
            if not p.is_relative_to(root) or not p.is_file() or p.stat().st_size != item['size'] or digest(p) != item['sha256']:
                return False
    except (OSError, KeyError, ValueError):
        return False
    return True


def stage_files(run, stage, lineage):
    """Explicit ownership prevents unrelated output damage invalidating audit."""
    r = Path(run)
    if stage == 'audit':
        paths = [r/'cache'/f'{lineage}_audited_cells.parquet']
        paths += [p for p in (r/'audit').glob(f'{lineage}_*')
                  if not any(k in p.name for k in ['prediction_snapshot', 'clock_', 'genemodel_join_qc', 'gene_set_provenance', 'unparsed_stable'])]
    elif stage == 'merge':
        paths = [r/'cache'/f'{lineage}_cells.parquet']
        paths += [r/'tables'/f'{lineage}_{s}' for s in ['cell_scores','schema.json','key_mapping.csv','score_qc.csv','cohort.csv','donor_celltype_scores.parquet','donor_celltype_sub_scores.parquet']]
    elif stage == 'subtypes':
        paths = [r/'tables'/f'{lineage}_{s}.csv' for s in ['subtype_ranking','subtype_contrasts','confounder_support']]
    elif stage == 'age':
        paths = list((r/'tables').glob(f'{lineage}_age_*')) + list((r/'tables').glob(f'{lineage}_composition*'))
    elif stage == 'concordance':
        paths = [r/'tables'/f'{lineage}_{s}.csv' for s in ['SID_recommendation_landscape','hUSI_recSID_associations','true_score_concordance']]
    elif stage == 'clock':
        paths = list((r/'audit').glob(f'{lineage}_clock_*'))
        paths += [r/'audit'/f'{lineage}_{s}.csv' for s in ['prediction_snapshot','genemodel_join_qc']]
        paths += list((r/'tables').glob(f'{lineage}_clock_*'))
    elif stage == 'genes':
        names = ['stable_features_parsed','tested_feature_universe','feature_subtype_integration','score_feature_gene_overlap','feature_evidence','pathway_enrichment','pseudobulk_score_gene_associations','pseudobulk_coverage']
        paths = [r/'tables'/f'{lineage}_{s}.csv' for s in names]
        paths += [r/'audit'/f'{lineage}_{s}.csv' for s in ['gene_set_provenance','unparsed_stable_features']]
        paths += [r/'results'/f'{lineage}_subtype_rank_comparison.json']
    elif stage == 'figure':
        paths = list((r/'figures').glob('*')) + [r/'results'/s for s in ['biological_findings.md','statistical_methods.md','limitations.md']]
    else:
        paths = []
    return paths


def required_outputs(stage, lineage):
    return {
        'audit': [f'cache/{lineage}_audited_cells.parquet',f'audit/{lineage}_mapping_qc.csv',f'audit/{lineage}_h5ad_schema.json'],
        'merge': [f'cache/{lineage}_cells.parquet',f'tables/{lineage}_cohort.csv'],
        'subtypes': [f'tables/{lineage}_subtype_ranking.csv'],
        'age': [f'tables/{lineage}_age_associations.csv',f'tables/{lineage}_age_donor_scores.parquet'],
        'clock': [f'audit/{lineage}_clock_provenance.json'],
        'genes': [f'tables/{lineage}_stable_features_parsed.csv'],
        'figure': ['figures/Fig_CellAge_Stromal_Integration.png','figures/Fig_CellAge_Stromal_Integration.pdf','figures/source_data/panel_manifest.csv'],
    }.get(stage, [])


def seal_stage(run, stage, lineage, scope=None):
    run = Path(run)
    marker = run/'checkpoints'/f'{lineage}_{stage}.json'
    state = json.loads(marker.read_text())
    if state.get('status') not in TERMINAL:
        raise RuntimeError(f'Stage did not commit a terminal diagnosis: {marker}')
    if state['status'] == 'complete':
        for name in required_outputs(stage, lineage):
            if not (run/name).is_file():
                raise RuntimeError(f'Completed stage missing required output: {name}')
    entries = artifacts(stage_files(run, stage, lineage)+[marker], run)
    atomic_json(run/'checkpoints/receipts'/f'{lineage}_{stage}.json',
                {'version':1,'stage':stage,'lineage':lineage,'scope':scope,'artifacts':entries})


def stage_valid(run, stage, lineage, scope=None):
    r = Path(run)
    try:
        marker = json.loads((r/'checkpoints'/f'{lineage}_{stage}.json').read_text())
        receipt = json.loads((r/'checkpoints/receipts'/f'{lineage}_{stage}.json').read_text())
        return (marker['status'] in TERMINAL and receipt.get('scope') == scope
                and bool(receipt['artifacts']) and valid_artifacts(receipt['artifacts'],r))
    except (OSError, ValueError, KeyError):
        return False


class WorkStore:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.hits = 0
        self.computed = 0

    def _paths(self, key):
        token = hashlib.sha256(json.dumps(key,sort_keys=True,default=_default).encode()).hexdigest()
        return self.root/(token+'.json'), self.root/(token+'.receipt.json')

    def value(self, key, compute, rng=None):
        data, marker = self._paths(key)
        try:
            receipt = json.loads(marker.read_text())
            if valid_artifacts(receipt['artifacts'],self.root):
                payload = json.loads(data.read_text())
                if rng is not None:
                    rng.bit_generator.state = payload['rng_after']
                self.hits += 1
                LOG.info('Resume work unit %s', str(key)[:150])
                return payload['value']
        except (OSError, ValueError, KeyError):
            pass
        value = compute()
        atomic_json(data, {'key':key,'value':value,'rng_after':rng.bit_generator.state if rng is not None else None})
        atomic_json(marker, {'artifacts':artifacts([data],self.root)})
        self.computed += 1
        return value

    def frame(self, key, compute):
        import pandas as pd
        data, marker = self._paths(['dataframe',key])
        data = data.with_suffix('.parquet')
        try:
            receipt = json.loads(marker.read_text())
            if valid_artifacts(receipt['artifacts'],self.root):
                self.hits += 1
                LOG.info('Resume table unit %s',str(key)[:150])
                return pd.read_parquet(data)
        except (OSError, ValueError, KeyError):
            pass
        value = compute()
        from .common import prepare_parquet
        with atomic_path(data) as tmp:
            prepare_parquet(value).to_parquet(tmp,index=False)
        atomic_json(marker, {'artifacts':artifacts([data],self.root)})
        self.computed += 1
        return value

    def action(self, key, paths, compute):
        _, marker = self._paths(['action',key])
        try:
            receipt = json.loads(marker.read_text())
            if receipt.get('artifacts') and valid_artifacts(receipt['artifacts'],self.root.parent.parent.parent):
                self.hits += 1
                LOG.info('Resume output unit %s',str(key)[:150])
                return
        except (OSError, ValueError, KeyError):
            pass
        compute()
        # action paths belong to the run, rather than units/<stage>/.
        run = self.root.parent.parent.parent
        entries = artifacts(paths,run)
        if not entries:
            raise RuntimeError(f'Output unit produced no files: {key}')
        atomic_json(marker, {'artifacts':entries})
        self.computed += 1


def cached_value(key, compute, rng=None):
    store = ACTIVE.get()
    return store.value(key,compute,rng) if store is not None else compute()


def cached_frame(key, compute):
    store = ACTIVE.get()
    return store.frame(key,compute) if store is not None else compute()


def cached_action(key, paths, compute):
    store = ACTIVE.get()
    return store.action(key,paths,compute) if store is not None else compute()


@contextmanager
def stage_work(run, stage, lineage):
    run = Path(run)
    path = run/'checkpoints/execution'/f'{lineage}_{stage}.json'
    store = WorkStore(run/'checkpoints/units'/f'{lineage}_{stage}')
    token = ACTIVE.set(store)
    atomic_json(path,{'state':'running','pid':os.getpid(),'stage':stage,'lineage':lineage,'started_utc':datetime.now(timezone.utc).isoformat()})
    try:
        yield store
    except BaseException as e:
        atomic_json(path,{'state':'interrupted' if isinstance(e,(KeyboardInterrupt,SystemExit)) else 'failed',
                          'stage':stage,'lineage':lineage,'exception':repr(e),'work_hits':store.hits,'work_computed':store.computed})
        raise
    else:
        atomic_json(path,{'state':'complete','stage':stage,'lineage':lineage,'work_hits':store.hits,'work_computed':store.computed})
    finally:
        ACTIVE.reset(token)


def atomic_partitioned_frame(df,path,partition_cols):
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_name('.'+path.name+'.'+uuid.uuid4().hex+'.partial')
    from .common import prepare_parquet
    prepare_parquet(df).to_parquet(tmp,partition_cols=partition_cols,index=False)
    if path.exists():
        os.replace(path,path.with_name('.'+path.name+'.previous.'+uuid.uuid4().hex))
    os.replace(tmp,path)


def save_figure(fig,path,**kwargs):
    with atomic_path(path) as tmp:
        fig.savefig(tmp,format=Path(path).suffix.lstrip('.'),**kwargs)


def invalidated_stages(cfg, source, current_inputs):
    """Propagate only changed input dependencies, never discard valid Atlas work
    just because a benchmark prediction CSV gained another completed fold.
    """
    source = Path(source)
    old = {p:(size,mtime) for p,size,mtime in json.loads((source/'audit/input_signature.json').read_text())}
    new = {p:(size,mtime) for p,size,mtime in current_inputs}
    changed = {p for p in old.keys() | new.keys() if old.get(p) != new.get(p)}
    invalid = set()
    all_stages = ['audit','merge','subtypes','age','concordance','clock','genes']
    for lin,conf in cfg['lineages'].items():
        atlas_inputs = {conf[k] for k in ['atlas','age_atlas','raw_atlas']} | {s['path'] for s in conf['scores'].values()}
        if atlas_inputs & changed:
            invalid.update((lin,s) for s in all_stages)
        g = conf['genemodel']
        if g['predictions'] in changed:
            invalid.add((lin,'clock'))
        if {g[k] for k in ['stable_features','feature_universe','gene_set_config']} & changed:
            invalid.add((lin,'genes'))
    if changed:
        invalid.add(('both','figure'))
    return invalid, sorted(changed)


def migrate_run(cfg, source, target, current_inputs, phases, legacy=False):
    """Independent copies of reusable artifacts, with invalid stage outputs omitted."""
    source,target = Path(source),Path(target)
    invalid,changed = invalidated_stages(cfg,source,current_inputs)
    excluded = set()
    for who,stage in invalid:
        for p in stage_files(source,stage,who):
            if p.is_dir():excluded.update(x for x in p.rglob('*') if x.is_file())
            else:excluded.add(p)
    global_audit = {'input_signature.json','input_manifest.csv','audit_summary.md','data_dictionary.md',
                    'sample_qc.csv','mapping_qc.csv','age_checks.csv','score_schema.csv','available_analyses.json'}
    for sub in ['audit','tables','cache','figures','results']:
        for p in (source/sub).rglob('*'):
            retained_embedding = any(p.name.endswith(f'_{lin}_embedding.csv') and (lin,'audit') not in invalid
                                     for lin in cfg['lineages'])
            if not p.is_file() or (p in excluded and not retained_embedding) or p.name.endswith('.partial'):
                continue
            if sub == 'audit' and p.name in global_audit:
                continue
            dest = target/sub/p.relative_to(source/sub)
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(p,dest)
    # Refresh lightweight prediction schema if a CSV grew; retain Atlas audit.
    for lin,conf in cfg['lineages'].items():
        g = conf['genemodel']
        if g['predictions'] in changed and (lin,'audit') not in invalid:
            from .common import read_table
            d=read_table(g['predictions'])
            atomic_json(target/'audit'/f'{lin}_genemodel_schema.json',{
                'path':g['predictions'],'shape':d.shape,'columns':list(d),
                'splits':d[g['split_column']].value_counts().to_dict(),
                'split_evidence':g['split_evidence'],'donor_disjoint_verified':False,
                'refresh_note':'Prediction input changed; Atlas audit reused, clock stage invalidated'})
    for who in ['fibro','endo','both']:
        for stage in phases:
            if (who=='both') != (stage=='figure') or (who,stage) in invalid:
                continue
            marker=source/'checkpoints'/f'{who}_{stage}.json'
            unit_dir=source/'checkpoints/units'/f'{who}_{stage}'
            if unit_dir.exists():
                shutil.copytree(unit_dir,target/'checkpoints/units'/unit_dir.name,dirs_exist_ok=True)
            if not marker.exists():continue
            scope=None
            if stage=='figure':
                import pandas as pd
                scope=pd.read_csv(source/'figures/source_data/panel_manifest.csv').lineage.drop_duplicates().tolist()
            if not legacy and not stage_valid(source,stage,who,scope):continue
            shutil.copy2(marker,target/'checkpoints'/marker.name)
            try:
                # Validate required legacy Parquet footers before establishing
                # their first checksum baseline (historical hashes do not exist).
                if legacy:
                    import pyarrow.parquet as pq
                    for name in required_outputs(stage,who):
                        if name.endswith('.parquet'):pq.read_metadata(target/name)
                seal_stage(target,stage,who,scope)
            except (OSError,RuntimeError,ValueError):
                pass
    atomic_json(target/'checkpoint_upgrade_provenance.json',{
        'source_run':str(source),'source_identity':json.loads((source/'run_identity.json').read_text()),
        'changed_inputs':changed,'invalidated_stages':[list(x) for x in sorted(invalid)],
        'reason':'Checkpoint/atomic-I/O upgrade or changed input dependencies; independent copies; numerical algorithms unchanged',
        'legacy_checksum_note':'Legacy artifacts acquire their first integrity baseline on migration; no claim of historical tamper detection' if legacy else None})
