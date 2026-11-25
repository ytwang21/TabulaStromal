"""Orchestration with explicit phase status, prerequisite checks and resumable units."""
from pathlib import Path
from dataclasses import dataclass
import importlib,json,logging,os,sys,time,traceback,shutil
from .storage import Cache,digest,file_info,file_hash,write_json,write_text,process_lock,annotation_identity
STAGES=['audit','features','benchmark','external','stability','enrichment','figures']


def scientific_config(cfg):
    """Runtime allocation is excluded only for a verified execution-only resume."""
    return {k:v for k,v in cfg.items() if k!='resources'}


def execution_identity(cfg,lineage,sources,inputs):
    fresh=digest(dict(config=cfg,lineage=lineage,code=sources,inputs=inputs))
    policy=cfg.get('resources',{}).get('execution_only_resume')
    if not policy:return fresh,dict(mode='normal')
    if sources!=policy['approved_sources']:raise RuntimeError('Execution-only resume refused: code differs from reviewed acceleration sources')
    manifest=json.loads(Path(policy['manifest']).read_text())['acceleration_resume']
    old=manifest['original_identities'][lineage]
    if scientific_config(cfg)!=scientific_config(old['config']):raise RuntimeError('Execution-only resume refused: scientific configuration changed')
    if inputs!=old['inputs']:raise RuntimeError('Execution-only resume refused: input/annotation/resource identity changed')
    if sys.executable!=old['python']:raise RuntimeError('Execution-only resume refused: Python environment changed')
    return old['signature'],dict(mode='verified_execution_only_resume',original_signature=old['signature'],
        execution_signature=fresh,scientific_configuration_unchanged=True,input_identities_unchanged=True,
        note='Reviewed scheduling/BLAS allocation changes only; existing artifact SHA checks remain enabled')


def task_context(ctx,index):
    """Assign an available logical GPU per task without changing model parameters."""
    import copy
    child=copy.copy(ctx);child.cfg={**ctx.cfg,'resources':dict(ctx.cfg['resources']),'jungclock':dict(ctx.cfg.get('jungclock',{}))}
    resources=child.cfg['resources'];gpus=[v for v in str(resources.get('physical_gpu','cpu')).split(',') if v!='cpu']
    if str(resources.get('device','cpu')).startswith('cuda') and gpus:
        device=f'cuda:{index%len(gpus)}';resources['device']=device;resources['xgb_device']=device
    return child


@dataclass
class Context:
    cfg:dict
    lineage:str
    resume:bool=True
    def __post_init__(self):
        self.root=Path(self.cfg['lineages'][self.lineage].get('output_dir',str(Path(self.cfg['output_root'])/self.lineage)));self.paths=self.cfg['lineages'][self.lineage]
        for d in ['audit','features','models','benchmark','external','stability','enrichment','figures','logs','reports','checkpoints']: (self.root/d).mkdir(parents=True,exist_ok=True)
        code_root=Path(__file__).resolve().parents[1]
        sources={str(p.relative_to(code_root)):file_hash(p) for p in sorted((code_root/'tabulastromal').glob('*.py'))}
        sources['run_tabulastromal.py']=file_hash(code_root/'run_tabulastromal.py')
        inputs=[]
        for p in [Path(self.paths['train']),Path(self.paths['external_combined'])]+sorted(Path(self.paths['external_gse_dir']).rglob('*.h5ad')):
            inputs.append(dict(**file_info(p),annotation_sha256=annotation_identity(p,self.cfg)) if p.exists() else dict(path=str(p),missing=True))
        for path in [self.cfg.get('jungclock',{}).get('ontology'),self.cfg['enrichment']['provided_sets'],self.cfg['enrichment'].get('hgnc_tsv'),self.cfg['enrichment'].get('stromal_sasp_evidence_csv')]:
            if path:inputs.append(file_info(path) if Path(path).exists() else dict(path=path,missing=True))
        for p in sorted(Path(self.cfg['enrichment']['gmt_dir']).glob('*')):
            if p.is_file():inputs.append(file_info(p))
        signature,execution=execution_identity(self.cfg,self.lineage,sources,inputs)
        self.cache=Cache(signature,self.resume)
        write_json(self.root/'reports/run_identity.json',dict(signature=self.cache.signature,code=sources,inputs=inputs,config=self.cfg,python=sys.executable,execution=execution))


def logger_for(root):
    log=logging.getLogger();log.setLevel(logging.INFO)
    for handler in log.handlers[:]:handler.close();log.removeHandler(handler)
    fmt=logging.Formatter('%(asctime)s | %(levelname)s | %(name)s | %(message)s')
    for handler in [logging.StreamHandler(sys.stdout),logging.FileHandler(root/'logs/pipeline.log')]:handler.setFormatter(fmt);log.addHandler(handler)
    logging.captureWarnings(True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def run_stage(ctx,stage):
    funcs={'audit':('data','audit'),'features':('data','features'),'benchmark':('training','benchmark'),
           'external':('external','external'),'stability':('stability','stability'),
           'enrichment':('biology','enrichment'),'figures':('figures','figures')}
    module,fn=funcs[stage];return getattr(importlib.import_module('tabulastromal.'+module),fn)(ctx)


def update_summary(ctx,status):
    import pandas as pd
    lines=[f'# Run summary: {ctx.lineage}',f'PID: {os.getpid()}',f'Python: {sys.executable}',f'Signature: {ctx.cache.signature}',
           '\nOnly completed phase markers attest execution. Files existing without a current matching checkpoint are not reported as completed.']
    for stage,info in status.items():lines.append(f'- {stage}: {info["status"]}'+(' — '+info['error'] if 'error' in info else ''))
    comp=ctx.root/'benchmark/completeness.json'
    if comp.exists():lines.append('\nBenchmark counts: '+comp.read_text())
    for f in ['stability/selection_decision.json','enrichment/limitations.json','figures/figure_manifest.json']:
        p=ctx.root/f
        if p.exists():lines.append('\n'+f+':\n'+p.read_text())
    lines+=['\n## Outputs and recovery',f'Lineage outputs: {ctx.root}',f'Input audit: {ctx.root}/audit/audit_report.md',
            f'Persistent stage/traceback log: {ctx.root}/logs/pipeline.log',f'Progress JSON: {ctx.root}/reports/run_status.json',
            'Restart the same launcher with --resume. Completed matching model/fold/adapter/bootstrap checkpoints are retained; an interrupted individual sklearn fit restarts from its last completed task. DeepSets resumes from its last completed epoch.',
            'No complete scientific findings are claimed until all requested stages complete. Missing mapping resources, empty stable sets, untestable GSEs and model failures remain explicit in CSV/reports.']
    write_text(ctx.root/'reports/run_summary.md','\n\n'.join(lines)+'\n')
    write_json(ctx.root/'reports/run_status.json',dict(lineage=ctx.lineage,pid=os.getpid(),signature=ctx.cache.signature,updated=time.time(),stages=status))


def run_lineage(cfg,lineage,stage='all',resume=True,stop_after=None):
    # Acquire the lock before writing any shared run identity or stage files.
    root=Path(cfg['lineages'][lineage].get('output_dir',str(Path(cfg['output_root'])/lineage)))
    with process_lock(root/'checkpoints/pipeline.lock'):
        ctx=Context(cfg,lineage,resume);logger_for(ctx.root);log=logging.getLogger(__name__)
        requested=STAGES if stage=='all' else [stage]
        status={s:dict(status='pending') for s in STAGES}
        for s in STAGES:
            if ctx.cache.valid(ctx.root/f'checkpoints/stage_{s}.json'):status[s]=dict(status='complete_cached')
        update_summary(ctx,status)
        for s in requested:
            marker=ctx.root/f'checkpoints/stage_{s}.json'
            if ctx.cache.valid(marker):log.info('RESUME stage %s already complete',s);continue
            # Stages can be rerun independently after their prerequisites complete.
            prerequisites={'audit':[],'features':['audit'],'benchmark':['audit','features'],
                           'external':['audit','features','benchmark'],'stability':['audit','features','benchmark','external'],
                           'enrichment':['stability'],'figures':['benchmark','external','stability','enrichment']}
            verification=Cache(ctx.cache.signature,True)
            missing=[p for p in prerequisites[s] if not verification.valid(ctx.root/f'checkpoints/stage_{p}.json')]
            if missing:raise RuntimeError(f'{s} requires completed current-signature stages: {missing}; run --stage all --resume')
            status[s]=dict(status='running',started=time.time());update_summary(ctx,status);log.info('STAGE %s START',s)
            try:
                outputs=run_stage(ctx,s)
                # Include child task markers to detect deleted/corrupt model/cache files during resume.
                patterns={'benchmark':'models/K*/*/fold*/complete.json','external':'external/tasks/*/*/*/*/*/complete.json','stability':'stability/*/bootstrap/*.complete.json'}
                if s in patterns:
                    children=[p for p in ctx.root.glob(patterns[s]) if verification.valid(p)]
                    if s=='external' and ctx.cfg['resources'].get('external_scope')=='combined':
                        children=[p for p in children if p.relative_to(ctx.root/'external/tasks').parts[0]=='combined']
                    outputs=list(outputs)+children
                ctx.cache.mark(marker,outputs);status[s]=dict(status='complete',finished=time.time());update_summary(ctx,status);log.info('STAGE %s COMPLETE',s)
            except BaseException as e:
                status[s]=dict(status='interrupted' if isinstance(e,(KeyboardInterrupt,SystemExit)) else 'failed',error=str(e),traceback=traceback.format_exc())
                update_summary(ctx,status);log.exception('STAGE %s FAILED; checkpoints retained',s);raise
            if stop_after==s:log.info('Intentional stop after %s; remaining stages pending',s);break
        update_summary(ctx,status)
        if status['figures']['status'].startswith('complete'):
            from .delivery import combined_figures
            combined_figures(ctx)
        return status
