"""End-to-end synthetic smoke test; never mixed with biological findings."""
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import yaml


def run_smoke(base_cfg, output_root):
    root = Path(output_root)/'smoke_tests'
    root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='synthetic_', dir=root))
    rng = np.random.default_rng(42)
    cfg = dict(base_cfg)
    cfg.update(project_root=str(work), output_root=str(work/'output'), bootstrap_n=100, min_donors=12,
               min_gse=2, embedding_max_cells=1000, figure_dpi=90, max_memory_gb=8,
               donor_column='sample', donor_crosswalk=None, donor_identity_evidence='Synthetic donor labels; smoke test only',
               expression={'enabled':False}, lineages={})
    for lineage in ['fibro', 'endo']:
        directory = work/lineage
        directory.mkdir()
        genes = ['CDKN1A', 'CDKN2A', 'SERPINE1', 'GLB1', 'LMNB1'] + [f'GENE_{i}' for i in range(15)]
        types = ['Type_with_under_score', 'Type_B']
        records = []
        for donor in range(48):
            age = float(rng.uniform(5, 90))
            for subtype in types:
                for cell in range(25):
                    records.append({'sample':f'D{donor:03d}', 'GSE':f'GSE{donor%3}', 'GSM':f'GSM{donor}', 'Age':age, 'age':age,
                        'Celltype':subtype, 'Celltype_sub':subtype, 'Celltype_general':subtype,
                        'tissue':['lung','skin'][donor%2], 'sex':['female','male'][donor%2],
                        'Disease_Status':'healthy' if donor%4 else 'disease', 'treat':'untreated', 'platform':'10X',
                        'total_counts':1000.+donor, 'pct_counts_mt':2.+donor/100})
        obs = pd.DataFrame(records, index=[f'{lineage}_cell_{i}' for i in range(len(records))])
        a = ad.AnnData(sparse.csr_matrix(rng.poisson(2, size=(len(obs), len(genes))).astype(np.float32)), obs=obs, var=pd.DataFrame(index=genes))
        a.obsm['X_umap'] = rng.normal(size=(len(obs),2)).astype(np.float32)
        atlas = directory/'atlas.h5ad'
        a.write_h5ad(atlas)
        h = obs[['sample','Celltype','Celltype_sub','Celltype_general','Age']].copy()
        h.insert(0, 'cellid', h.index)
        h['hUSI'] = np.clip(.2 + .003*h.Age + .1*h.Celltype.eq(types[0]) + rng.normal(0,.05,len(h)),0,1)
        h['husi_SSE'] = (h.hUSI>.45).astype(int)
        h['husi_GMM'] = 1 + (h.hUSI>.35).astype(int) + (h.hUSI>.5).astype(int)
        hp = directory/'husi.csv'; h.to_csv(hp,index=False)
        s = h[['cellid','sample','Celltype','Age']].copy()
        for i in range(1,7): s[f'SID{i}'] = rng.normal(i/3, .4, len(s))
        sp=directory/'sencid.tsv'; s.to_csv(sp,sep='\t',index=False)
        predictions=[]
        dm=obs[['sample','Age']].drop_duplicates()
        folds={d:i%3 for i,d in enumerate(dm['sample'])}
        for fold in range(3):
            for _,r in dm.iterrows():
                predictions.append({'sample':r['sample'],'Age':r.Age,'prediction':r.Age+rng.normal(0,8),
                    'model':'ElasticNet','K':8000,'fold':fold,'split':'Valid' if folds[r['sample']]==fold else 'Train'})
        pp=directory/'predictions.csv';pd.DataFrame(predictions).to_csv(pp,index=False)
        feats=[]
        for t in types:
            for j,gene in enumerate(genes):
                feats.append({'Feature':t+'_'+gene,'Celltype':t,'Gene':gene,'Mean_coefficient':.05 if j%2 else -.05,
                              'Selection_frequency':.9 if j<8 else .1,'In_candidate_panel':True})
        f=pd.DataFrame(feats); fp=directory/'stable.csv'; f.loc[f.Selection_frequency>.7].to_csv(fp,index=False)
        up=directory/'universe.csv'; f.to_csv(up,index=False)
        gmt=directory/'synthetic.gmt';gmt.write_text('Synthetic_markers\tunit-test\t'+'\t'.join(genes[:8])+'\nSynthetic_background\tunit-test\t'+'\t'.join(genes[8:])+'\n')
        gp=directory/'genes.yaml';gp.write_text(yaml.safe_dump({'gene_sets':{'Synthetic':{'path':str(gmt),'version':'fixture','source':'synthetic smoke test'}}}))
        cfg['lineages'][lineage]={'atlas':str(atlas),'age_atlas':str(atlas),'raw_atlas':str(atlas),
          'scores':{'hUSI':{'path':str(hp),'cell_column':'cellid','values':{'hUSI':'hUSI','husi_SSE':'husi_SSE','husi_GMM':'husi_GMM'},'semantics':'senescence_score','source_script':str(Path(__file__))},
                    'SenCID':{'path':str(sp),'cell_column':'cellid','values':{f'SID{i}':f'SenCID_rec_SID{i}' for i in range(1,7)},'semantics':'SID_recommendation_not_senescence_prediction','source_script':str(Path(__file__))}},
          'genemodel':{'predictions':str(pp),'model_column':'model','model_value':'ElasticNet','feature_count_column':'K','feature_count':8000,
                       'split_column':'split','allowed_splits':['Valid'],'sample_column':'sample','age_column':'Age','prediction_column':'prediction','fold_column':'fold',
                       'split_evidence':'synthetic split','donor_disjoint_verified':False,'preprocessing_train_only_verified':True,
                       'stable_features':str(fp),'feature_universe':str(up),'gene_set_config':str(gp),'score_gene_lists':{}}}
    config = work/'config.synthetic.yaml'
    config.write_text(yaml.safe_dump(cfg,sort_keys=False))
    command = [sys.executable,str(Path(__file__).parents[1]/'run_pipeline.py'),'--config',str(config),'--phase','all','--n-jobs','1']
    subprocess.run(command,check=True)
    import json
    run=Path(json.loads((work/'output/latest_run.json').read_text())['path'])
    for lin in ['fibro','endo']:
        qc=pd.read_csv(run/'audit'/f'{lin}_mapping_qc.csv')
        assert qc.passed.all() and qc.match_rate.eq(1).all()
        agg=pd.read_parquet(run/'tables'/f'{lin}_donor_celltype_scores.parquet')
        assert agg.donor_key.nunique()==48 and len(agg)==96
        slopes=pd.read_csv(run/'tables'/f'{lin}_age_associations.csv')
        primary=slopes.loc[slopes.sensitivity.eq('min_cells_20') & slopes.adjusted.eq(True) & slopes.status.eq('ok')]
        assert len(primary)==2 and primary.effect.gt(0).all()
        assert all(str(c).startswith('SenCID_rec_') for c in pd.read_csv(run/'tables'/f'{lin}_SID_recommendation_landscape.csv').recommendation.unique())
        clock=json.loads((run/'checkpoints'/f'{lin}_clock.json').read_text())
        assert clock['status']=='complete'
    assert (run/'figures/Fig_CellAge_Stromal_Integration.pdf').stat().st_size>1000
    print('SYNTHETIC_SMOKE_PASSED=' + str(run))
    return run
