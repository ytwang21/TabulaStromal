"""Exercise actual process death and verify committed work is not repeated."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pytest

from scripts.checkpointing import (ACTIVE, WorkStore, atomic_path, cached_value,
    seal_stage, stage_valid, stage_work, invalidated_stages)
from scripts.common import save_table, stage_status, write_json


def test_atomic_output_keeps_previous_file_if_writer_dies(tmp_path):
    p=tmp_path/'result.csv';p.write_text('previous complete output')
    with pytest.raises(RuntimeError):
        with atomic_path(p) as partial:
            partial.write_text('incomplete')
            raise RuntimeError('interruption')
    assert p.read_text()=='previous complete output'
    with atomic_path(p) as partial:
        partial.write_text('new complete output')
    assert p.read_text()=='new complete output'


def test_corrupt_work_cache_recomputes_only_that_unit(tmp_path):
    store=WorkStore(tmp_path/'units');calls=[]
    def work():
        calls.append(1);return {'effect':.123,'p':float('nan')}
    first=store.value('one',work)
    assert store.value('one',work)['effect']==first['effect'] and len(calls)==1
    data,_=store._paths('one');data.write_text('broken')
    assert store.value('one',work)['effect']==.123 and len(calls)==2


def test_rng_position_is_restored_when_bootstrap_unit_is_skipped(tmp_path):
    s=WorkStore(tmp_path);rng=np.random.default_rng(42)
    values=s.value('bootstrap',lambda:rng.normal(size=20).tolist(),rng=rng)
    next_value=rng.random()
    rng=np.random.default_rng(42)
    restored=s.value('bootstrap',lambda:pytest.fail('should reuse completed unit'),rng=rng)
    assert restored==values and rng.random()==next_value


def test_missing_or_modified_output_invalidates_completed_stage(tmp_path):
    save_table(pd.DataFrame({'value':[1]}),tmp_path/'tables/fibro_subtype_ranking.csv')
    stage_status(tmp_path,'subtypes','fibro','complete')
    seal_stage(tmp_path,'subtypes','fibro')
    assert stage_valid(tmp_path,'subtypes','fibro')
    p=tmp_path/'tables/fibro_subtype_ranking.csv';p.write_text('value\n2\n')
    assert not stage_valid(tmp_path,'subtypes','fibro')
    p.unlink()
    assert not stage_valid(tmp_path,'subtypes','fibro')


def test_frame_cache_survives_interrupted_stage(tmp_path):
    with pytest.raises(KeyboardInterrupt):
        with stage_work(tmp_path,'merge','fibro') as store:
            store.frame('donor_aggregate',lambda:pd.DataFrame({'n':[20],'mean':[.5]}))
            raise KeyboardInterrupt()
    assert json.loads((tmp_path/'checkpoints/execution/fibro_merge.json').read_text())['state']=='interrupted'
    with stage_work(tmp_path,'merge','fibro') as store:
        d=store.frame('donor_aggregate',lambda:pytest.fail('expensive aggregate should not run again'))
    assert d['n'].iloc[0]==20 and store.hits==1


def test_sigkill_and_restart_preserves_first_completed_unit(tmp_path):
    root=Path(__file__).resolve().parents[1]
    code=r'''
from pathlib import Path
import sys,time
from scripts.checkpointing import stage_work,cached_value
r=Path(sys.argv[1]);mode=sys.argv[2]
def first():
    with (r/'first.log').open('a') as f:f.write('computed\n')
    return {'n':42}
def second():
    with (r/'second.log').open('a') as f:f.write('attempt\n')
    (r/'ready').write_text('ready')
    if mode=='interrupt':time.sleep(60)
    return {'done':True}
with stage_work(r,'age','fibro'):
    assert cached_value('first',first)['n']==42
    assert cached_value('second',second)['done']
'''
    env={**os.environ,'PYTHONPATH':str(root)}
    proc=subprocess.Popen([sys.executable,'-c',code,str(tmp_path),'interrupt'],env=env)
    try:
        limit=time.monotonic()+15
        while not (tmp_path/'ready').exists():
            if proc.poll() is not None or time.monotonic()>limit:
                pytest.fail('worker did not reach interruption point')
            time.sleep(.02)
        proc.kill();proc.wait(timeout=5)
    finally:
        if proc.poll() is None:proc.kill();proc.wait(timeout=5)
    subprocess.run([sys.executable,'-c',code,str(tmp_path),'resume'],env=env,check=True,timeout=15)
    assert (tmp_path/'first.log').read_text().splitlines()==['computed']
    assert len((tmp_path/'second.log').read_text().splitlines())==2


def test_cli_resume_skips_finished_lineage_and_reuses_interrupted_unit(tmp_path,monkeypatch):
    import yaml
    import run_pipeline as entry
    from scripts import step00_audit as module
    f=tmp_path/'input.txt';f.write_text('fixed input')
    conf={'output_root':str(tmp_path/'output'),'project_root':str(tmp_path),'lineages':{}}
    for lin in ['fibro','endo']:
        conf['lineages'][lin]={'atlas':str(f),'age_atlas':str(f),'raw_atlas':str(f),
            'scores':{},'genemodel':{k:str(f) for k in ['predictions','stable_features','feature_universe','gene_set_config']}}
    config=tmp_path/'config.yaml';config.write_text(yaml.safe_dump(conf))
    attempts={'fibro':0,'endo':0,'unit':0}
    def audit(cfg,run,lin):
        attempts[lin]+=1
        if lin=='endo':
            def unit():attempts['unit']+=1;return 123
            assert cached_value('costly audit subunit',unit)==123
            if attempts[lin]==1:raise KeyboardInterrupt()
        save_table(pd.DataFrame({'cell_key':['C1']}),run/'cache'/f'{lin}_audited_cells.parquet')
        save_table(pd.DataFrame({'passed':[True]}),run/'audit'/f'{lin}_mapping_qc.csv')
        write_json(run/'audit'/f'{lin}_h5ad_schema.json',{'n_obs':1})
        stage_status(run,'audit',lin,'complete')
    monkeypatch.setattr(module,'run',audit)
    monkeypatch.setattr(module,'inventory',lambda *a:None)
    monkeypatch.setattr(module,'finalize',lambda *a:None)
    argv=['--config',str(config),'--phase','audit','--lineage','both']
    previous=signal.getsignal(signal.SIGTERM)
    try:
        with pytest.raises(KeyboardInterrupt):entry.run_pipeline(argv)
        first=json.loads((tmp_path/'output/latest_run.json').read_text())['path']
        run=entry.run_pipeline(argv+['--resume'])
        assert str(run)==first and attempts=={'fibro':1,'endo':2,'unit':1}
        assert stage_valid(run,'audit','fibro') and stage_valid(run,'audit','endo')
    finally:
        signal.signal(signal.SIGTERM,previous)


def test_new_prediction_fold_does_not_invalidate_atlas_or_other_lineage(tmp_path):
    cfg={'lineages':{}}
    rows=[]
    for lin in ['fibro','endo']:
        c={'atlas':lin+'_atlas','age_atlas':lin+'_age','raw_atlas':lin+'_raw','scores':{},
           'genemodel':{k:lin+'_'+k for k in ['predictions','stable_features','feature_universe','gene_set_config']}}
        cfg['lineages'][lin]=c
        rows.extend((p,10,100) for p in [c[k] for k in ['atlas','age_atlas','raw_atlas']]+list(c['genemodel'].values()))
    write_json(tmp_path/'audit/input_signature.json',rows)
    current=[(p,20 if p=='endo_predictions' else size,101 if p=='endo_predictions' else mtime) for p,size,mtime in rows]
    invalid,changed=invalidated_stages(cfg,tmp_path,current)
    assert changed==['endo_predictions']
    assert invalid=={('endo','clock'),('both','figure')}
