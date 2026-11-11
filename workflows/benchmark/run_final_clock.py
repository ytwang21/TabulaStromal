#!/usr/bin/env python
"""Complete Phase3/4 from audited checkpoints without rerunning 13-model training."""
import argparse,importlib,json,logging,os,subprocess,sys
from pathlib import Path
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/endo_matplotlib'))

def main():
    from final_clock_common import Context,STAGES,own_valid,mark,LOG,file_hash
    from tabulastromal.storage import process_lock,write_json
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--lineage',choices=['Fibro','Endo','both'],default='both')
    p.add_argument('--stage',choices=STAGES+['all','validate'],default='all');p.add_argument('--config');p.add_argument('--resume',action='store_true')
    p.add_argument('--workers',type=int);p.add_argument('--smoke-test',action='store_true')
    a=p.parse_args();stages=STAGES if a.stage=='all' else [a.stage]
    if a.stage=='validate':
        from validate_final_clock import validate
        validate(a.lineage);return 0
    if a.lineage=='both':
        if a.config:p.error('--config applies to a single lineage')
        for stage in stages:
            jobs=[]
            for lineage in ['Fibro','Endo']:
                cmd=[sys.executable,'-u',str(ROOT/'run_final_clock.py'),'--lineage',lineage,'--stage',stage]
                if a.resume:cmd+=['--resume']
                if a.workers:cmd+=['--workers',str(a.workers)]
                if a.smoke_test:cmd+=['--smoke-test']
                from final_clock_common import ROOTS
                log=(ROOTS[lineage]/'logs/final_clock_pipeline.log').open('ab')
                jobs.append((lineage,subprocess.Popen(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT),log))
            codes={lineage:job.wait() for lineage,job,_ in jobs}
            for _,_,log in jobs:log.close()
            print(json.dumps(dict(stage=stage,exit_codes=codes)),flush=True)
            if any(codes.values()):return 2
        if a.stage in ['all','figures'] and not a.smoke_test:
            from final_clock_figures import combine,write_final_report
            combine();write_final_report()
            from validate_final_clock import validate
            validate('both')
        return 0
    ctx=Context(a.lineage,a.config,a.resume)
    if a.workers:ctx.options['workers']=a.workers
    logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s')
    modules=dict(audit='final_clock_audit',bootstrap='final_clock_stability',final_model='final_clock_training',
                 external='final_clock_external',enrichment='final_clock_enrichment',figures='final_clock_figures')
    with process_lock(ctx.path('checkpoints/final_clock.lock')):
        for stage in stages:
            module=importlib.import_module(modules[stage]);code=Path(module.__file__)
            marker=ctx.path(f'checkpoints/final_clock_{stage}.complete.json')
            if own_valid(ctx,marker,code) and ctx.json(str(marker.relative_to(ctx.root))).get('common_sha256')==file_hash(ROOT/'final_clock_common.py'):
                LOG.info('%s %s valid checkpoint reused',ctx.lineage,stage);continue
            LOG.info('%s STAGE %s START',ctx.lineage,stage)
            try:
                outputs=module.run(ctx,smoke=a.smoke_test)
                if not a.smoke_test:mark(ctx,marker,outputs,code)
                failure=ctx.path('reports/final_clock_failure.json')
                if failure.exists():
                    prior=json.loads(failure.read_text())
                    if prior.get('stage')==stage and prior.get('status')!='resolved':
                        from final_clock_common import now
                        write_json(failure,dict(prior,status='resolved',resolved_at=now()))
                LOG.info('%s STAGE %s COMPLETE artifacts=%d',ctx.lineage,stage,len(outputs))
            except Exception as e:
                LOG.exception('STAGE %s FAILED; previous artifacts retained',stage)
                write_json(ctx.path('reports/final_clock_failure.json'),dict(stage=stage,error=str(e)));return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
