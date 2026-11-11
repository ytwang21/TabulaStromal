#!/usr/bin/env python
"""Resume saved lineage pipelines from external with one fixed FT30 evaluation."""
from pathlib import Path
import argparse,json,os,sys,subprocess,time
ROOT=Path(__file__).resolve().parent

def package_path(lineage):
    package=ROOT/'Fibro_Accelerated_20261102' if lineage=='Fibro' else ROOT
    sys.path[:0]=[str(package/'vendor'),str(package),str(ROOT)]

# Spawn workers must import the same package/classes as the corresponding driver.
if os.environ.get('TABULASTROMAL_FT30_LINEAGE'):
    package_path(os.environ['TABULASTROMAL_FT30_LINEAGE'])

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lineage',choices=['both','Fibro','Endo'],default='both')
    p.add_argument('--ft30-config',default=str(ROOT/'ft30_config.yaml'))
    p.add_argument('--config',help='Optional scientific config for a single lineage/test')
    p.add_argument('--stage',choices=['all','external','stability','enrichment','figures'],default='all')
    p.add_argument('--workers',type=int)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--smoke-test',action='store_true')
    p.add_argument('--external-eval-mode',choices=['adapt30_all100','no_adapt_all100'])
    args=p.parse_args(argv)
    import yaml
    options=yaml.safe_load(Path(args.ft30_config).read_text())
    if args.workers:options['workers_per_lineage']=args.workers
    if args.external_eval_mode:options['external_eval_mode']=args.external_eval_mode
    if options['adaptation_fraction']!=.30 or options['repeats']!=1 or options['adaptation_seed']!=42:
        p.error('This protocol requires fraction=.30, one fixed split, seed=42')
    if args.lineage=='both':
        if args.config:p.error('--config requires one lineage')
        if args.dry_run:
            print(json.dumps(options,indent=2));return 0
        children=[]
        for lineage in ['Fibro','Endo']:
            command=[sys.executable,'-u',str(Path(__file__).resolve()),'--lineage',lineage,
                '--ft30-config',str(Path(args.ft30_config).resolve()),'--stage',args.stage,'--resume']
            if args.workers:command+=['--workers',str(args.workers)]
            if args.smoke_test:command+=['--smoke-test']
            if args.external_eval_mode:command+=['--external-eval-mode',args.external_eval_mode]
            log=(ROOT/'logs'/f'{lineage}_ft30.nohup.log').open('ab')
            child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
            log.close();children.append((lineage,child))
            print(f'{lineage} FT30 driver PID={child.pid}',flush=True)
        codes={lin:child.wait() for lin,child in children}
        print(json.dumps(codes),flush=True)
        return 0 if all(v==0 for v in codes.values()) else 2
    lineage=args.lineage
    os.environ['TABULASTROMAL_FT30_LINEAGE']=lineage
    package_path(lineage)
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    os.environ['CUDA_VISIBLE_DEVICES']=','.join(map(str,options['gpu_ids'])) if options['gpu_enabled'] else ''
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runtime/endo_matplotlib'))
    cfgpath=Path(args.config) if args.config else ROOT/options['lineages'][lineage]['config']
    cfg=yaml.safe_load(cfgpath.read_text())
    out=Path(cfg['output_root'])/lineage
    if args.dry_run:
        print(json.dumps(dict(lineage=lineage,scientific_config=str(cfgpath),output=str(out),
            options=options,stages='external FT30 -> cached/resumed stability -> enrichment -> FT30 figures',
            retrain_development_models=False,new_versioned_result_directory=False),indent=2));return 0
    from tabulastromal import pipeline
    from tabulastromal.storage import process_lock,write_json,file_hash,digest
    import ft30,ft30_figures
    ft30.configure(cfg,lineage,options)
    with process_lock(out/'checkpoints/pipeline.lock'):
        oldidentity=json.loads((out/'reports/run_identity.json').read_text())
        ctx=pipeline.Context(cfg,lineage,True)
        if ctx.cache.signature!=oldidentity['signature']:
            write_json(out/'reports/run_identity.json',oldidentity)
            raise RuntimeError('Scientific signature changed; refusing to invalidate saved development models')
        pipeline.logger_for(out);log=__import__('logging').getLogger(__name__)
        for stage in ['audit','features','benchmark']:
            if not ctx.cache.valid(out/f'checkpoints/stage_{stage}.json'):
                raise RuntimeError(f'Requires completed {stage}; never retrains development models')
        protocol_file=out/'external/ft30_protocol.json'
        formal_driver=out/'reports/ft30_driver.json'
        reusable=False
        if args.smoke_test and protocol_file.exists() and formal_driver.exists():
            prior=json.loads(formal_driver.read_text());record=json.loads(protocol_file.read_text())
            reusable=(prior.get('source_signature')==ctx.cache.signature and
                prior.get('source',{}).get('ft30.py')==file_hash(ROOT/'ft30.py') and
                record['options']['protocol']==options['protocol'])
        protocol=json.loads(protocol_file.read_text())['ft30_protocol'] if reusable else ft30.prepare(ctx)
        write_json(out/'reports'/('ft30_smoke_driver.json' if args.smoke_test else 'ft30_driver.json'),dict(pid=os.getpid(),lineage=lineage,ft30_protocol=protocol,
            source_signature=ctx.cache.signature,options=options,config=str(cfgpath),
            source={n:file_hash(ROOT/n) for n in ['ft30.py','ft30_figures.py','ft30_figure_panels.py','run_tabulastromal_ft30.py']},
            policy='Same original directories; no development CV/model or scientific config changes'))
        if args.smoke_test:
            # A smoke test must never replace the full 260-task aggregate CSVs.
            item=next(i for i in ft30.selected_folds(ctx) if i[0]==cfg['features']['K'][0] and i[1]=='ElasticNet')
            ft30.task(cfg,item,protocol)
            _,target,prefix=ft30.source_and_target(ctx,item)
            from tabulastromal.storage import write_csv
            from tabulastromal.metrics import evaluate,METRICS
            import pandas as pd,numpy as np
            predictions=pd.read_csv(target/'predictions_ft30.csv')
            split=json.loads((target/'split_ft30.json').read_text())
            assert set(predictions.sample_id)==set(split['evaluation_ids'])
            assert not set(split['fit_ids'])&set(split['holdout_ids'])
            assert set(split['external_fit_ids'])==set(split['adaptation_ids'])
            assert predictions.used_for_adaptation.sum()==len(split['adaptation_ids'])
            for source,dest,column in [('metrics_ft30.csv','external_combined_ft30_smoke_per_fold.csv','prediction'),
                                       ('baseline_metrics_ft30.csv','external_combined_ft30_smoke_baseline.csv','baseline_prediction')]:
                row=pd.read_csv(target/source);computed=evaluate(predictions.Age,predictions[column])
                for metric in METRICS:np.testing.assert_allclose(computed[metric],row.iloc[0][metric],rtol=1e-8,atol=1e-8,equal_nan=True)
                write_csv(out/'external'/dest,row)
            write_csv(out/'external/external_combined_ft30_smoke_predictions.csv',predictions)
            write_csv(out/'external/external_combined_ft30_smoke_heldout.csv',pd.read_csv(target/'heldout_metrics_ft30.csv'))
            ft30_figures.benchmark_plots(ctx,partial=True)
            write_json(out/'reports/ft30_smoke_test.json',dict(status='passed',lineage=lineage,
                ft30_protocol=protocol,model='ElasticNet',K=item[0],fold=item[2],evaluation_scope='all100',
                predictions_validated=True,figure_data_source_verified=True,
                adaptation_strategy=pd.read_csv(target/'metrics_ft30.csv').iloc[0].adaptation_strategy,
                adaptation_status=pd.read_csv(target/'metrics_ft30.csv').iloc[0].adaptation_status,
                n_all=len(predictions),n_heldout=len(split['holdout_ids']),EN_reference_adapter_verified=True))
            return 0
        figure_code=digest({n:file_hash(ROOT/n) for n in ['ft30_figures.py','ft30_figure_panels.py']})
        def stage_valid(s):
            marker=out/f'checkpoints/stage_{s}.json'
            # Do not scan 271,445 obsolete grid outputs before replacing External.
            if s in ['external','figures']:
                if not marker.exists() or json.loads(marker.read_text()).get('ft30_protocol')!=protocol:return False
            if s=='figures' and json.loads(marker.read_text()).get('figure_code')!=figure_code:return False
            return ctx.cache.valid(marker)
        status={s:dict(status='complete_cached' if stage_valid(s) else 'pending') for s in pipeline.STAGES}
        requested=['external','stability','enrichment','figures'] if args.stage=='all' else [args.stage]
        def update():pipeline.update_summary(ctx,status)
        update()
        for stage in requested:
            marker=out/f'checkpoints/stage_{stage}.json'
            needs=(not stage_valid(stage))
            if stage in ['external','figures']:
                needs=needs or json.loads(marker.read_text()).get('ft30_protocol')!=protocol if marker.exists() else True
                if stage=='figures' and marker.exists():
                    needs=needs or json.loads(marker.read_text()).get('external_eval_mode')!=options['external_eval_mode']
            if not needs:
                log.info('FT30 RESUME %s complete',stage);continue
            status[stage]=dict(status='running',started=time.time());update()
            try:
                if stage=='external':
                    outputs=ft30.run_external(ctx,protocol)
                    # Publish benchmark/K figures immediately; biology need not finish first.
                    ft30_figures.backup_benchmarks(ctx)
                    ft30_figures.benchmark_plots(ctx)
                elif stage=='stability':
                    if status['external']['status'] not in ['complete','complete_cached']:raise RuntimeError('external FT30 required')
                    from tabulastromal.stability import stability
                    outputs=stability(ctx)+list(out.glob('stability/*/bootstrap/*.complete.json'))
                elif stage=='enrichment':
                    if not ctx.cache.valid(out/'checkpoints/stage_stability.json'):raise RuntimeError('stability required')
                    from tabulastromal.biology import enrichment
                    outputs=enrichment(ctx)
                else:
                    if not ctx.cache.valid(out/'checkpoints/stage_stability.json') or not ctx.cache.valid(out/'checkpoints/stage_enrichment.json'):
                        raise RuntimeError('stability and enrichment required')
                    outputs=ft30.final_external(ctx,protocol)+ft30_figures.render(ctx)
                if marker.exists() and not json.loads(marker.read_text()).get('ft30_protocol'):
                    backup=marker.with_name(marker.stem+'.pre_ft30.json')
                    if not backup.exists():__import__('shutil').copy2(marker,backup)
                ctx.cache.mark(marker,outputs,ft30_protocol=protocol,external_eval_mode=options['external_eval_mode'],
                               **({'figure_code':figure_code} if stage=='figures' else {}))
                status[stage]=dict(status='complete',finished=time.time());update()
                log.info('FT30 STAGE %s COMPLETE',stage)
            except BaseException as e:
                status[stage]=dict(status='failed',error=str(e),traceback=__import__('traceback').format_exc());update();log.exception('FT30 stage failed');raise
        if args.stage=='all':
            write_json(out/'reports/ft30_complete.json',dict(status='complete',ft30_protocol=protocol,
                signature=ctx.cache.signature,finished=time.time(),pid=os.getpid(),stages=status))
        return 0

if __name__=='__main__':raise SystemExit(main())
