"""Final-clock no-adapt predictions and separately identified final-origin FT30."""
from concurrent.futures import ThreadPoolExecutor
import re
from final_clock_common import *

def metadata_with_counts(ctx,name,meta):
    m=meta.copy();counts=ctx.csv('audit/'+name+'_sample_subtype_counts.csv')
    known=set(ctx.json('models/final_feature_schema.json')['celltypes'])
    n=counts[counts.Celltype.isin(known)].groupby('sample_id').n_cells.sum()
    m['n_cells']=m.sample_id.map(n).fillna(0).astype(int)
    m['age_true']=pd.to_numeric(m.Age,errors='raise');m['donor_id']=''
    m['donor_identity_status']='unavailable; group_id is a GSE/GSM proxy'
    m['lineage']=ctx.lineage
    return m

def task(ctx,item,clocks):
    name=item['dataset'];stem='final_clock_task_'+re.sub(r'[^a-zA-Z0-9_-]','_',name)
    marker=ctx.path('external/'+stem+'.complete.json')
    pp=ctx.path('external/'+stem+'_predictions.csv');mp=ctx.path('external/'+stem+'_metrics.csv')
    if own_valid(ctx,marker,__file__):return pd.read_csv(pp,keep_default_na=False),pd.read_csv(mp),[]
    raw=sp.load_npz(ctx.path('features/'+name+'/pseudobulk.npz'))
    meta=metadata_with_counts(ctx,name,ctx.csv('features/'+name+'/samples.csv'))
    if len(meta)==0:return pd.DataFrame(),pd.DataFrame(),[dict(dataset=name,path=item['path'],reason='no eligible labeled samples after original QC/overlap exclusions')]
    schema=ctx.json('features/'+name+'/schema.json');lookup=pd.Index(schema['features'])
    rows=[];metrics=[]
    combined=ctx.csv('features/combined/samples.csv');combined_ids=set(combined.sample_id);combined_gsm=set(combined.GSM)
    for version,b in clocks.items():
        ix=lookup.get_indexer(b['selected_features']);assert (ix>=0).all()
        pred=bundle_predict(b,raw[:,ix].toarray())
        q=meta.copy();q['age_pred']=pred;q['prediction']=pred;q['age_acceleration_raw']=pred-q.age_true
        q['model_version']=version;q['final_model_id']=b['model_id'];q['final_model_type']='ElasticNet'
        q['K']=ctx.K;q['dataset']=name;q['method']=item['method'];q['n_features']=len(ix)
        q['evaluation_source']='final_clock_external_noadapt' if version=='final_clock' else 'EN_version_comparison_noadapt'
        q['overlaps_combined']=q.sample_id.isin(combined_ids)|q.GSM.isin(combined_gsm)
        q['interpretation']='chronological age prediction; age_acceleration_raw is an unadjusted residual'
        rows.append(q)
        metrics.append(dict(lineage=ctx.lineage,dataset=name,method=item['method'],model_version=version,
            final_model_id=b['model_id'],K=ctx.K,n_features=len(ix),evaluation_source=q.evaluation_source.iloc[0],
            n_overlapping_combined=int(q.overlaps_combined.sum()),**score(meta.age_true,pred)))
    predictions=pd.concat(rows,ignore_index=True);metric=pd.DataFrame(metrics)
    write_csv(pp,predictions);write_csv(mp,metric)
    pins=[ctx.path('features/'+name+'/'+f) for f in ['pseudobulk.npz','samples.csv','schema.json']]
    pins += [ctx.path('models/'+f) for f in ['final_clock.joblib','EN_full_train_original_features.joblib','EN_full_train_stable_features.joblib','EN_best_fold_reference.joblib']]
    mark(ctx,marker,[pp,mp]+pins,__file__)
    LOG.info('%s no-adapt %s: %d samples; all model versions evaluated',ctx.lineage,name,len(meta))
    return predictions,metric,[]

def run(ctx,smoke=False):
    decision=versions(ctx)
    clocks={name:load_final_clock(info['path']) for name,info in decision['versions'].items()}
    clocks['EN_best_fold_reference']=load_final_clock(ctx.path('models/EN_best_fold_reference.joblib'))
    clocks['final_clock']=load_final_clock(ctx.path('models/final_clock.joblib'))
    items=ctx.json('audit/inventory.json')[1:];outputs=[]
    if smoke:
        available=[i for i in items if i['dataset']!='combined' and len(ctx.csv('features/'+i['dataset']+'/samples.csv'))>0]
        item=min(available,key=lambda i:len(ctx.csv('features/'+i['dataset']+'/samples.csv')))
        name=item['dataset'];b=clocks['final_clock']
        cached=predict_age(ctx.path('features/'+name),b)
        actual=predict_age(item['path'],b,expression=item['expression']).set_index('sample_id')
        p=actual.loc[cached.sample_id,'age_pred'].to_numpy();delta=float(np.max(abs(cached.age_pred.to_numpy()-p)))
        np.testing.assert_allclose(cached.age_pred,p,atol=.005,rtol=1e-4)
        target=ctx.path('reports/final_clock_external_smoke.json')
        write_json(target,dict(status='passed',smoke_only=True,production_complete=False,dataset=name,n=len(cached),
            actual_h5ad_path=item['path'],cache_vs_h5ad_max_difference_years=delta,expression=item['expression']))
        return [target]
    with ThreadPoolExecutor(max_workers=min(ctx.options['io_workers'],len(items))) as pool:
        results=list(pool.map(lambda it:task(ctx,it,clocks),items))
    preds=pd.concat([r[0] for r in results if len(r[0])],ignore_index=True)
    metrics=pd.concat([r[1] for r in results if len(r[1])],ignore_index=True)
    skipped=[q for r in results for q in r[2]]
    def save(rel,frame):
        p=ctx.path(rel);write_csv(p,frame);outputs.append(p)
    final=preds[preds.model_version=='final_clock'];fm=metrics[metrics.model_version=='final_clock']
    combined=final[final.dataset=='combined'].copy()
    assert not combined.sample_id.duplicated().any()
    save('external/final_clock_external_combined_predictions.csv',combined)
    save('external/final_clock_external_combined_metrics.csv',fm[fm.dataset=='combined'])
    save('external/final_clock_per_gse_predictions.csv',final[final.dataset!='combined'])
    save('external/final_clock_per_gse_metrics.csv',fm[fm.dataset!='combined'])
    save('external/final_clock_skipped_gse.csv',pd.DataFrame(skipped,columns=['dataset','path','reason']))
    save('external/final_clock_model_comparison_external.csv',metrics[metrics.dataset=='combined'])
    save('external/final_clock_all_version_external_predictions.csv',preds)
    split=ctx.csv('external/external_ft30_split.csv');allids=set(split.sample_id)
    assert set(combined.sample_id)==allids and not set(split.loc[split.set=='adaptation','group_id'])&set(split.loc[split.set=='holdout','group_id'])
    assert not set(combined.sample_id)&set(clocks['final_clock']['training_ids'])
    b=clocks['final_clock'];ds=dataset(ctx,'combined');ix=pd.Index(ds['schema']['features']).get_indexer(b['selected_features']);assert (ix>=0).all()
    raw=ds['X'][:,ix].toarray().astype(np.float32);m=metadata_with_counts(ctx,'combined',ds['meta'])
    adapt=m.sample_id.isin(split.loc[split.set=='adaptation','sample_id']).to_numpy()
    assert int(adapt.sum())==int((split.set=='adaptation').sum())
    np.testing.assert_allclose(m.age_true,split.set_index('sample_id').loc[m.sample_id,'Age'],rtol=0,atol=1e-10)
    before=bundle_predict(b,raw);baseline=combined.set_index('sample_id').loc[m.sample_id,'age_pred'].to_numpy()
    np.testing.assert_allclose(before,baseline,atol=1e-10,rtol=1e-10)
    updated=copy.deepcopy(b);updated['model']=fit_checked(b['model'],b['preprocessor'].transform(raw[adapt]),m.age_true.to_numpy()[adapt])
    updated['model_id']=b['model_id']+'_FT30_regularized_refit'
    updated['training_scope']='external_adaptation_only_regularized_refit'
    updated['training_ids']=m.sample_id[adapt].tolist();updated['source_full_development_model_id']=b['model_id']
    updated['adaptation_strategy']='regularized_refit_adapt_only_EN_reference_not_incremental'
    adapter=ctx.path('models/final_clock_FT30.joblib');dump(adapter,updated);outputs.append(adapter)
    after=bundle_predict(updated,raw);np.testing.assert_allclose(after,bundle_predict(load_final_clock(adapter),raw),atol=1e-10,rtol=1e-10)
    q=m.copy();q['age_pred']=after;q['prediction']=after;q['age_pred_noadapt']=before;q['age_acceleration_raw']=after-q.age_true
    q['adaptation_membership']=np.where(adapt,'adapt30','heldout70');q['used_for_adaptation']=adapt
    q['evaluation_source']='external_ft30_all_100';q['final_model_id']=updated['model_id'];q['source_full_development_model_id']=b['model_id']
    q['final_model_type']='ElasticNet';q['K']=ctx.K;q['adaptation_strategy']=updated['adaptation_strategy']
    assert len(q)==len(allids) and set(q.sample_id)==allids and not q.sample_id.duplicated().any()
    save('external/external_ft30_all_100_predictions.csv',q)
    common=dict(lineage=ctx.lineage,K=ctx.K,n_adapt=int(adapt.sum()),n_heldout=int((~adapt).sum()),
                final_model_id=updated['model_id'],source_full_development_model_id=b['model_id'],adaptation_strategy=updated['adaptation_strategy'])
    all_metric=dict(common,evaluation_source='external_ft30_all_100',performance_note='includes fitted30%; mixed apparent/unseen, NOT independent external validation',**score(m.age_true,after))
    save('external/external_ft30_all_100_metrics.csv',pd.DataFrame([all_metric]))
    held=[]
    for state,pred in [('final_clock_noadapt',before),('final_clock_FT30_refit',after)]:
        held.append(dict(common,model_version=state,evaluation_source='external_ft30_heldout70',**score(m.age_true.to_numpy()[~adapt],pred[~adapt])))
    save('external/external_ft30_heldout70_metrics.csv',pd.DataFrame(held))
    save('external/external_ft30_heldout70_predictions.csv',q[~adapt])
    protocols=[dict(common,model_version='final_clock_noadapt',evaluation_source='final_clock_noadapt_all100',**score(m.age_true,before)),*held,all_metric]
    save('external/final_clock_external_protocol_comparison.csv',pd.DataFrame(protocols))
    sensitivity=[]
    for c in ['GSE','tissue','sex','disease_global','Disease_Status']:
        if c not in combined:continue
        for name,part in combined.groupby(c,dropna=False):
            sensitivity.append(dict(covariate=c,level=name,evaluation_source='primary_noadapt_stratum_descriptive',
                independence_note='strata overlap; GSE-level integration alternatives are not independent cohorts',**score(part.age_true,part.age_pred)))
    save('external/final_clock_prediction_covariate_sensitivity.csv',pd.DataFrame(sensitivity))
    p=ctx.path('reports/final_clock_external_validation.json')
    write_json(p,dict(status='passed',n_combined=len(combined),n_adapt=int(adapt.sum()),n_heldout=int((~adapt).sum()),
        all100_same_pairs_for_all_metrics=True,n_per_gse_datasets=int((fm.dataset!='combined').sum()),
        n_per_gse_sample_representations=int((final.dataset!='combined').sum()),n_skipped=len(skipped),
        adapter_source='this final full-development clock; never old best fold',
        adaptation_definition=updated['adaptation_strategy'],primary_external_labels_used_for_training=False,
        secondary_all100_includes_fitted30=True,negative_predictions_retained=True,donor_independence_verified=False))
    outputs.append(p)
    return outputs+[ctx.path('models/model_versions.json'),ctx.path('models/final_clock.joblib'),ctx.path('external/external_ft30_split.csv')]
