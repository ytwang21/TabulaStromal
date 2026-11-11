#!/usr/bin/env python
"""Acceptance tests of actual final-clock artifacts, not simulated results."""
from final_clock_common import *
from scipy.stats import hypergeom
from tabulastromal.biology import bh

def check_metrics(frame,y,p):
    s=score(y,p)
    assert int(frame['n'])==len(y)
    for k in METRICS:np.testing.assert_allclose(float(frame[k]),s[k],rtol=1e-8,atol=1e-8,equal_nan=True)

def validate(lineage='both'):
    from PIL import Image
    from pypdf import PdfReader
    Image.MAX_IMAGE_PIXELS=None
    results={}
    modules=dict(audit='final_clock_audit',bootstrap='final_clock_stability',final_model='final_clock_training',external='final_clock_external',enrichment='final_clock_enrichment',figures='final_clock_figures')
    for l in ['Fibro','Endo'] if lineage=='both' else [lineage]:
        c=Context(l);LOG.info('Validating actual final clock %s',l)
        for stage in STAGES:
            marker=c.path(f'checkpoints/final_clock_{stage}.complete.json')
            assert own_valid(c,marker,ROOT/(modules[stage]+'.py')),str(marker)
            assert json.loads(marker.read_text())['common_sha256']==file_hash(ROOT/'final_clock_common.py')
        pins=c.json('audit/final_clock_source_pins.json')
        assert all(file_hash(p['path'])==p['sha256'] for p in pins),'Historical source changed'
        b=load_final_clock(c.path('models/final_clock.joblib'));v=versions(c)
        train=c.csv('features/train/samples.csv');ext=c.csv('features/combined/samples.csv')
        assert b['training_scope']=='full_development' and set(b['training_ids'])==set(train.sample_id)
        assert b['model_id']==v['versions'][v['selected_version']]['model_id']
        assert not set(b['training_ids'])&set(ext.sample_id)
        stable=c.csv('stability/stable_celltype_gene_features.csv');cells=c.csv('stability/stable_celltypes_ranked.csv');genes=c.csv('stability/stable_genes_unique.csv')
        assert len(stable)==stable.feature.nunique() and (stable.selection_frequency>=.7).all()
        assert set(stable.Celltype)==set(cells.Celltype) and set(stable.Gene)==set(genes.Gene)
        count=c.csv('stability/stable_feature_counts_summary.csv');count=count[(count.scope=='full_train')&(count.threshold==.7)].iloc[0]
        assert [len(stable),len(cells),len(genes)]==[int(count.n_stable_Celltype_Gene),int(count.n_stable_Celltypes),int(count.n_unique_stable_Genes)]
        manifest=c.csv('stability/bootstrap_runs_manifest.csv')
        assert len(manifest)==1800 and manifest.groupby('scope').size().eq(300).all() and manifest.status.eq('success').all()
        assert all(file_hash(r.path)==r.sha256 for r in manifest.itertuples())
        cv=pd.read_csv(c.path('stability/final_clock_internal_cv_metrics.csv'));oof=c.csv('stability/final_clock_internal_oof_predictions.csv')
        for (model,fold),part in oof.groupby(['model_version','fold']):
            name='original_EN' if model=='original_EN' else 'stable_EN_0.70'
            row=cv[(cv.model_version==name)&(cv.fold==fold)].iloc[0]
            check_metrics(row,part.Age,part.prediction)
            info=c.json(f'models/K{c.K}/ElasticNet/fold{fold}/model.json')
            assert set(part.sample_id)==set(info['valid_ids'])
            assert not set(part.sample_id)&set(info['train_ids'])
        orig=cv[cv.model_version=='original_EN'];st=cv[cv.model_version=='stable_EN_0.70']
        assert v['selected_version']=='original_EN' if st.MAE.mean()>orig.MAE.mean() else True
        p=c.csv('external/final_clock_external_combined_predictions.csv');m=c.csv('external/final_clock_external_combined_metrics.csv').iloc[0]
        assert len(p)==len(ext) and set(p.sample_id)==set(ext.sample_id) and not p.sample_id.duplicated().any()
        actual=predict_age(c.path('features/combined'),b).set_index('sample_id')
        np.testing.assert_allclose(p.age_pred,actual.loc[p.sample_id,'age_pred'],rtol=1e-10,atol=1e-10)
        check_metrics(m,p.age_true,p.age_pred)
        ft=c.csv('external/external_ft30_all_100_predictions.csv');fm=c.csv('external/external_ft30_all_100_metrics.csv').iloc[0]
        assert set(ft.sample_id)==set(p.sample_id) and len(ft)==len(p) and not ft.sample_id.duplicated().any()
        adapter=load_final_clock(c.path('models/final_clock_FT30.joblib'))
        assert adapter['source_full_development_model_id']==b['model_id']
        assert adapter['feature_schema']==b['feature_schema']
        for field in ['mean_','scale_','var_']:np.testing.assert_array_equal(getattr(adapter['preprocessor'],field),getattr(b['preprocessor'],field))
        np.testing.assert_allclose(predict_age(c.path('features/combined'),adapter).set_index('sample_id').loc[ft.sample_id,'age_pred'],ft.age_pred,rtol=1e-10,atol=1e-10)
        split=c.csv('external/external_ft30_split.csv');adapt=set(split.loc[split.set=='adaptation','sample_id']);held=set(split.loc[split.set=='holdout','sample_id'])
        assert set(adapter['training_ids'])==adapt and not adapt&held
        assert set(ft.loc[ft.used_for_adaptation,'sample_id'])==adapt
        check_metrics(fm,ft.age_true,ft.age_pred)
        hm=c.csv('external/external_ft30_heldout70_metrics.csv');h=ft[ft.sample_id.isin(held)]
        for row in hm.to_dict('records'):
            check_metrics(row,h.age_true,h.age_pred if row['model_version']=='final_clock_FT30_refit' else h.age_pred_noadapt)
        # The new best-fold reference must reproduce its historical no-adapt
        # predictions, showing final no-adapt weakness is not a column mix-up.
        reference=c.csv('external/final_clock_all_version_external_predictions.csv')
        reference=reference[(reference.dataset=='combined')&(reference.model_version=='EN_best_fold_reference')]
        old=c.csv('external/external_ft30_all_folds_predictions.csv')
        old=old[(old.model=='ElasticNet')&(old.K==c.K)&(old.fold==v['best_reference_fold'])].set_index('sample_id')
        np.testing.assert_allclose(reference.age_pred,old.loc[reference.sample_id,'baseline_prediction'],rtol=1e-7,atol=1e-6)
        for path in ['enrichment/final_clock_pathway_enrichment.csv','enrichment/aging_gene_set_overlap.csv']:
            f=c.csv(path)
            np.testing.assert_allclose(f.p_value,hypergeom.sf(f.overlap_count.to_numpy()-1,f.universe_size.to_numpy(),f.set_size.to_numpy(),f.query_size.to_numpy()),rtol=1e-8,atol=1e-12)
            for _,part in f.groupby(['scope','library']):np.testing.assert_allclose(part.FDR,bh(part.p_value),rtol=1e-8,atol=1e-12)
        ag=c.csv('enrichment/aging_gene_set_overlap.csv');ag=ag[ag.scope=='global']
        assert set(ag.term)=={'SASP','SenMayo','CellAge','SigRS','SenUp','AgingAtlas','ASIG','GenAge'}
        assert ag.source_label.eq('provided_unverified').all()
        figure=c.json('figures/final_clock_figure_manifest.json')
        assert figure['final_model_id']==b['model_id'] and figure['external_primary']=='no adaptation'
        for path in map(Path,figure['figures']):
            if path.suffix=='.pdf':assert all(len(page.images)==0 for page in PdfReader(path).pages),str(path)
            else:
                with Image.open(path) as image:assert min(image.info['dpi'])>=300
        assert c.json('reports/final_clock_external_smoke.json')['status']=='passed'
        assert c.json('reports/final_clock_training_smoke.json')['status']=='passed'
        results[l]=dict(status='passed',n_bootstrap_per_scope=300,n_actual_bootstrap_records=1800,
            stable_Celltype_Gene=len(stable),stable_Celltypes=len(cells),stable_Genes=len(genes),
            final_model_id=b['model_id'],n_development=len(train),n_external=len(p),n_adapt=len(adapt),n_heldout=len(held),
            final_bundle_reload_prediction_verified=True,original_source_files_unchanged=True,OOF_provenance_verified=True,
            best_fold_reference_reproduces_historical_noadapt=True,all100_same_pairs_for_all_metrics=True,
            final_origin_FT30_and_frozen_preprocessing_verified=True,hypergeometric_and_BH_verified=True,vector_PDFs_verified=True)
        write_json(c.path('reports/final_clock_validation.json'),results[l])
    if lineage=='both':
        out=ROOT/'Stromal_Combined_Figure/outputs';done=json.loads((out/'FinalClock_FINAL_v2.completion.json').read_text())
        assert done['status']=='complete'
        assert all(file_hash(x['path'])==x['sha256'] for x in done['artifacts'])
        assert all(file_hash(x['path'])==x['sha256'] for x in done['sources'])
        for x in done['artifacts']:
            if x['path'].endswith('.pdf'):assert all(len(page.images)==0 for page in PdfReader(x['path']).pages)
        assert len(done['panel_sources'])==18
    write_json(ROOT/'final_clock_validation.json',dict(status='passed',validated_at=now(),lineages=results))
    print(json.dumps(results,ensure_ascii=False,indent=2),flush=True)
    return results

if __name__=='__main__':validate()
