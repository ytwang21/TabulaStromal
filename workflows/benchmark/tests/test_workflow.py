"""Meaningful synthetic leakage, artifact and complete-stage integration checks."""
import copy,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import yaml
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tabulastromal.metrics import evaluate
from tabulastromal.models import Preprocessor,NAMES,fit_model,predict_model,save_bundle,load_bundle
from tabulastromal.training import grouped_splits
from tabulastromal.external import adaptation_split
from tabulastromal.storage import Cache,write_text


def toy_config(base):
    cfg=yaml.safe_load((ROOT/'config.yaml').read_text());cfg['output_root']=str(base/'output')
    cfg['resources']['threads']=1;cfg['models'].update(xgb_estimators=5,tree_estimators=5,boost_estimators=5,mlp_max_iter=10,
        deep_epochs=2,deep_hidden=4,deep_batch_samples=4,max_iter=5000,linear_alphas=[.1],l1_ratios=[.5])
    cfg['features'].update(K=[8],selection='variance',min_cells_per_sample=2,chunk_cells=17,cells_per_subtype=2,variance_threshold=1e-6)
    cfg['cv'].update(n_splits=2,inner_splits=2);cfg['stability'].update(n_bootstrap=2,min_fold_improvements=2)
    cfg['external'].update(repeats=1,fractions=[.1,.3,.7]);cfg['figures']['dpi']=300
    cfg['enrichment']['gmt_dir']=str(base/'gmts');Path(cfg['enrichment']['gmt_dir']).mkdir()
    (Path(cfg['enrichment']['gmt_dir'])/'toy_GO.gmt').write_text('test term\ttoy_generated_not_biological_evidence\tIL6\tTP53\tCDKN1A\tMMP1\tVEGFA\n')
    cfg['enrichment']['min_set_size']=1
    rng=np.random.default_rng(16);genes=['IL6','TP53','CDKN1A','MMP1','VEGFA','TAGLN','COL1A1','DCN','VIM','GAPDH','ACTB','CLDN5']
    def make(n,prefix,gse,repeat_samples=False):
        obs=[];expr=[]
        for donor in range(n):
            age=10+donor*5
            for sample in range(2 if repeat_samples else 1):
                for ct in ['Subtype_one_pan','Subtype_two_skin']:
                    for cell in range(3):
                        obs.append(dict(Age=float(age),sample=f'{prefix}{donor}_{sample}',GSE=gse if isinstance(gse,str) else gse[donor],GSM=f'{prefix}{donor}',Celltype=ct,tissue='skin',sex='female' if donor%2 else 'male'))
                        expr.append(np.maximum(0,rng.normal(1,.2,len(genes))+age/40*np.linspace(.1,1,len(genes))))
        return ad.AnnData(sp.csr_matrix(np.asarray(expr,np.float32)),obs=pd.DataFrame(obs,index=[f'{prefix}_c{i}' for i in range(len(obs))]),var=pd.DataFrame(index=genes))
    train=make(12,'TR','GSE_train',True);ext=make(12,'EXT',['GSE_small']*4+['GSE_large']*8)
    trainp=base/'train.h5ad';extp=base/'combined.h5ad';train.write_h5ad(trainp);ext.write_h5ad(extp)
    gsedir=base/'per_gse'
    for gse in ['GSE_small','GSE_large']:
        path=gsedir/gse;path.mkdir(parents=True)
        sub=ext[ext.obs.GSE==gse].copy();sub.obs['sample']=sub.obs['sample'].astype(str)+'_query';sub.raw=sub.copy();sub.X=sp.csr_matrix(sub.X.toarray()-2)
        sub.write_h5ad(path/'adata_ext_new_bbknn_labeled.h5ad')
        alt=sub.copy();alt.X=alt.raw.X.copy();alt.write_h5ad(path/'adata_ext_new_scanorama_labeled.h5ad')
    paths=dict(train=str(trainp),external_combined=str(extp),external_gse_dir=str(gsedir));cfg['lineages']={'Fibro':paths,'Endo':paths}
    return cfg


class Invariants(unittest.TestCase):
    def test_grouping_preprocessing_adaptation(self):
        rng=np.random.default_rng(4);X=sp.csr_matrix(rng.normal(size=(24,12)).astype(np.float32));y=np.repeat(np.arange(12)*5,2)
        meta=pd.DataFrame(dict(Age=y,group_id=np.repeat([f'g{i}' for i in range(12)],2),sample_id=[f's{i}' for i in range(24)],GSE=['one']*24))
        cfg=yaml.safe_load((ROOT/'config.yaml').read_text());cfg['features']['selection']='variance'
        for tr,va in grouped_splits(meta,5):
            self.assertFalse(set(meta.group_id.iloc[tr])&set(meta.group_id.iloc[va]))
            pre=Preprocessor.fit(X[tr],y[tr],[f'f{i}' for i in range(12)],meta.sample_id.iloc[tr],8,cfg)
            self.assertFalse(set(pre.fit_ids)&set(meta.sample_id.iloc[va]))
            np.testing.assert_allclose(pre.scaler.mean_,X[tr][:,pre.indices].toarray().mean(0),rtol=1e-6,atol=1e-7)
            mean=pre.scaler.mean_.copy();pre.transform(X[va]*100000);np.testing.assert_array_equal(pre.scaler.mean_,mean)
        for fraction in [.1,.3,.7]:
            split,reason=adaptation_split(meta,fraction,42)
            if split:
                a,t=split;self.assertFalse(set(meta.group_id.iloc[a])&set(meta.group_id.iloc[t]));self.assertFalse(set(a)&set(t))
        self.assertIsNone(adaptation_split(meta.iloc[:4],.1,42)[0])

    def test_metrics_and_checkpoint_integrity(self):
        m=evaluate([10,20,30],[-1,22,32]);self.assertEqual(m['MAD'],20/3);self.assertEqual(m['MedAE'],2)
        self.assertTrue(np.isnan(evaluate([2],[3])['R2']));self.assertTrue(np.isnan(evaluate([2,2],[3,4])['PCC']))
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'data';marker=Path(td)/'complete.json';write_text(p,'abc');cache=Cache('test');cache.mark(marker,[p]);self.assertTrue(cache.valid(marker))
            write_text(p,'xyz');self.assertFalse(cache.valid(marker))
            fresh=Cache('test',False);self.assertFalse(fresh.valid(marker));fresh.mark(marker,[p]);self.assertTrue(fresh.valid(marker))

    def test_all_thirteen_predict_and_reload(self):
        rng=np.random.default_rng(2);X=rng.normal(size=(16,8)).astype(np.float32);y=np.arange(16.)
        cfg=yaml.safe_load((ROOT/'config.yaml').read_text());cfg['resources']['threads']=1
        cfg['models'].update(xgb_estimators=3,tree_estimators=3,boost_estimators=3,mlp_max_iter=4,deep_epochs=1,deep_hidden=4,max_iter=5000)
        names=[f'ct_G{i}' for i in range(8)];pre=Preprocessor.fit(sp.csr_matrix(X),y,names,[str(i) for i in range(16)],8,{**cfg,'features':{**cfg['features'],'selection':'variance'}})
        Z=pre.transform(sp.csr_matrix(X));bags=[(sp.csr_matrix(np.vstack([z,z*.9])),np.array([0,0])) for z in Z]
        with tempfile.TemporaryDirectory() as td:
            for name in NAMES:
                model=fit_model(name,Z,y,cfg,42,bags=bags)
                pred=predict_model(model,name,Z,bags);self.assertEqual(pred.shape,(16,));self.assertTrue(np.isfinite(pred).all())
                save_bundle(Path(td)/name,name,model,pre,dict(seed=42));nn,loaded,pp,info=load_bundle(Path(td)/name,cfg)
                np.testing.assert_allclose(predict_model(loaded,nn,Z,bags),pred,rtol=1e-6)
                if name in ['DeepSets','Ensemble']:self.assertTrue((Path(td)/name/'state_dict.pt').exists())

    def test_end_to_end_and_resume(self):
        with tempfile.TemporaryDirectory(prefix='tabulastromal_smoke_') as td:
            base=Path(td);cfg=toy_config(base);config=base/'toy.yaml';config.write_text(yaml.safe_dump(cfg))
            command=[sys.executable,str(ROOT/'run_tabulastromal.py'),'--config',str(config),'--lineage','Fibro','--stage','all']
            run=subprocess.run(command,text=True,capture_output=True)
            log=ROOT/'tests/last_smoke.log';log.write_text(run.stdout+'\n'+run.stderr)
            self.assertEqual(run.returncode,0,'See '+str(log)+'\n'+run.stdout[-3000:]+'\n'+run.stderr[-2000:])
            out=base/'output/Fibro';cv=pd.read_csv(out/'benchmark/cv_fold_metrics.csv');self.assertEqual(set(cv.model),set(NAMES));self.assertEqual(len(cv),52)
            self.assertTrue(set(['R2','MAE','MAD','Age_MAD','MedAE','RMSE','PCC']).issubset(cv))
            merged=pd.read_csv(out/'external/external_combined_adaptation_all.csv');self.assertIn(0.,set(merged.adaptation_fraction))
            paired=merged[merged.adaptation_fraction>0];self.assertGreater(len(paired),0);self.assertIn('baseline_MAE',paired);self.assertIn('delta_MAE',paired)
            for splitfile in (out/'external/tasks').rglob('split.json'):
                split=json.loads(splitfile.read_text());self.assertFalse(set(split['adapt_groups'])&set(split['test_groups']))
            skips=pd.read_csv(out/'external/skipped_gse.csv');self.assertIn('GSE_small',set(skips.GSE))
            final=out/'models/final_model/model.json';before=final.stat().st_mtime_ns
            again=subprocess.run(command+['--resume'],text=True,capture_output=True);self.assertEqual(again.returncode,0,again.stderr);self.assertEqual(before,final.stat().st_mtime_ns)
            self.assertTrue((out/'figures/Fig3_TabulaStromal_Fibro.pdf').stat().st_size>1000)
            # Retain compact evidence without claiming toy values are real study results.
            evidence=dict(label='SYNTHETIC_SMOKE_ONLY',status='passed',benchmark_tasks=26,models=NAMES,cv_metric_rows=len(cv),paired_adaptation_rows=len(paired),stage_status=json.loads((out/'reports/run_status.json').read_text())['stages'],resume_preserved_model_mtime=True)
            (ROOT/'tests/smoke_evidence.json').write_text(json.dumps(evidence,indent=2))


if __name__=='__main__':unittest.main(verbosity=2)
