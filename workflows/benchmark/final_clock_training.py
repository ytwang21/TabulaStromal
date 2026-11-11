"""Reuse honest outer TRAIN stability fits; refit both full-development EN versions."""
from final_clock_common import *

def full_preprocessor(ctx,ds):
    ids=ds['meta'].sample_id.tolist()
    p=ctx.path(f'checkpoints/preprocessing/full_train/K{ctx.K}.joblib')
    if p.is_file() and source_valid(ctx,p.with_suffix('.complete.json')):
        old=joblib.load(p)
        if old.fit_ids==ids:return old.indices,old.scaler,str(p)
    for tag in ['full_train','full_train_final']:
        p=ctx.path(f'checkpoints/preprocessing/{tag}/ranking.joblib')
        if p.is_file() and source_valid(ctx,p.with_suffix('.complete.json')):
            # Rank is model-independent MI computed on ALL development only.
            proof=list(ctx.path(f'checkpoints/preprocessing/{tag}').glob('K*.joblib'))
            if not any(joblib.load(q).fit_ids==ids for q in proof):continue
            ix=np.sort(np.asarray(joblib.load(p),int)[:ctx.K])
            scaler=StandardScaler().fit(ds['X'][:,ix].toarray().astype(np.float32))
            return ix,scaler,str(p)
    raise RuntimeError('No validated full-development feature rank; refuse external-derived features')

def run(ctx,smoke=False):
    ds=dataset(ctx);meta=ds['meta'];y=meta.Age.to_numpy(float);outputs=[];metric_rows=[];oof=[];pins=[]
    def save(rel,frame):
        p=ctx.path(rel);write_csv(p,frame);outputs.append(p)
    cv=ctx.csv('benchmark/cv_fold_metrics.csv');v=cv[(cv.model=='ElasticNet')&(cv.K==ctx.K)&(cv.split=='Valid')]
    best=v.sort_values(['R2','MAE','fold'],ascending=[False,True,True]).iloc[0]
    bestpath=ctx.path(f'models/K{ctx.K}/ElasticNet/fold{int(best.fold)}')
    bestmodel=joblib.load(bestpath/'estimator.joblib');bestpre=joblib.load(bestpath/'preprocessor.joblib')
    params=ctx.json(f'checkpoints/tuning/full_train/K{ctx.K}_ElasticNet.json')
    assert set(params['tuning_ids'])==set(meta.sample_id)
    indices,scaler,rank_source=full_preprocessor(ctx,ds)
    if smoke:
        n=min(80,len(y));raw=ds['X'][:n,indices[:min(200,len(indices))]].toarray()
        scale=StandardScaler().fit(raw);m=copy.deepcopy(bestmodel);m.set_params(**params['params'])
        trained=fit_checked(m,scale.transform(raw),y[:n]);p=trained.predict(scale.transform(raw))
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(dir=ctx.path('models'),suffix='.joblib') as f:
            joblib.dump(trained,f.name);np.testing.assert_allclose(joblib.load(f.name).predict(scale.transform(raw)),p,atol=1e-10)
        target=ctx.path('reports/final_clock_training_smoke.json')
        write_json(target,dict(status='passed',smoke_only=True,production_complete=False,n=n,n_features=raw.shape[1],
            solver=trained.__class__.__module__,dual_gap=float(trained.dual_gap_),read_reload_verified=True))
        return [target]
    # Re-evaluate all saved source and fold-stable estimators on the same unseen IDs.
    for fold in range(1,6):
        source=ctx.path(f'models/K{ctx.K}/ElasticNet/fold{fold}');info=json.loads((source/'model.json').read_text())
        assert source_valid(ctx,source/'complete.json')
        valid=meta.sample_id.isin(info['valid_ids']).to_numpy();train=meta.sample_id.isin(info['train_ids']).to_numpy()
        assert not set(meta.loc[valid,'group_id'])&set(meta.loc[train,'group_id'])
        source_model=joblib.load(source/'estimator.joblib');source_pre=joblib.load(source/'preprocessor.joblib')
        assert set(source_pre.fit_ids)==set(info['train_ids'])
        original_pred=source_model.predict(source_pre.transform(ds['X'][valid]))
        metric_rows.append(dict(model_version='original_EN',threshold=np.nan,fold=fold,K=ctx.K,
            n_features=len(source_pre.indices),evaluation_source='outer_valid_conditional_selected_K',**score(y[valid],original_pred)))
        for th in [ctx.options['threshold'],ctx.options['sensitivity_threshold']]:
            path=ctx.path(f'stability/outer{fold}/threshold_{th:.2f}')
            if not (path/'estimator.joblib').exists():
                metric_rows.append(dict(model_version=f'stable_EN_{th:.2f}',threshold=th,fold=fold,K=ctx.K,
                    n_features=0,evaluation_source='outer_valid_conditional_selected_K',status='no_stable_features',**{k:np.nan for k in METRICS}))
                continue
            m=joblib.load(path/'estimator.joblib');pre=joblib.load(path/'preprocessor.joblib')
            assert set(pre.fit_ids)==set(info['train_ids'])
            table=pd.read_csv(ctx.path(f'stability/outer{fold}/bootstrap_all_features.csv'),usecols=['feature','selection_frequency'])
            assert set(pre.names)==set(table.loc[table.selection_frequency>=th,'feature'])
            pred=m.predict(pre.transform(ds['X'][valid]));s=score(y[valid],pred)
            metric_rows.append(dict(model_version=f'stable_EN_{th:.2f}',threshold=th,fold=fold,K=ctx.K,
                n_features=len(pre.indices),evaluation_source='outer_valid_conditional_selected_K',**s))
            pins += [path/'estimator.joblib',path/'preprocessor.joblib']
            if th==ctx.options['threshold']:
                q=meta.loc[valid].copy();q['prediction']=pred;q['model_version']='stable_EN';q['fold']=fold;oof.append(q)
        q=meta.loc[valid].copy();q['prediction']=original_pred;q['model_version']='original_EN';q['fold']=fold;oof.append(q)
        pins += [source/'estimator.joblib',source/'preprocessor.joblib']
    metrics=pd.DataFrame(metric_rows);predictions=pd.concat(oof,ignore_index=True)
    for name,q in predictions.groupby('model_version'):
        assert len(q)==len(meta) and not q.sample_id.duplicated().any() and set(q.sample_id)==set(meta.sample_id)
    save('stability/final_clock_internal_cv_metrics.csv',metrics)
    save('stability/final_clock_internal_oof_predictions.csv',predictions)
    summary=metrics.groupby('model_version')[METRICS].agg(['mean','std','count']);summary.columns=['_'.join(k) for k in summary.columns]
    save('stability/final_clock_internal_cv_summary.csv',summary.reset_index())
    primary=metrics[metrics.model_version==f'stable_EN_{ctx.options["threshold"]:.2f}'].set_index('fold')
    original=metrics[metrics.model_version=='original_EN'].set_index('fold')
    paired=original[['MAE','R2']].join(primary[['MAE','R2']],lsuffix='_original',rsuffix='_stable')
    improved=(paired.MAE_stable<paired.MAE_original)&(paired.R2_stable>=paired.R2_original-ctx.options['stable_promotion_r2_tolerance'])
    promote=(int(improved.sum())>=ctx.options['stable_promotion_min_fold_improvements'] and primary.MAE.notna().all()
        and primary.MAE.mean()<original.MAE.mean() and primary.R2.mean()>=original.R2.mean()-ctx.options['stable_promotion_r2_tolerance'])
    selected='stable_EN' if promote else 'original_EN'
    save('stability/final_clock_paired_cv_comparison.csv',paired.reset_index().assign(improved=improved.to_numpy()))
    stable=ctx.csv('stability/stable_celltype_gene_features.csv');lookup=pd.Index(ds['schema']['features'])
    stable_indices=lookup.get_indexer(stable.feature);assert (stable_indices>=0).all()
    assert set(meta.sample_id).isdisjoint(set(ctx.csv('features/combined/samples.csv').sample_id))
    model_specs={};development=[]
    for name,ix,scale in [('original_EN',indices,scaler),('stable_EN',stable_indices,None)]:
        raw=ds['X'][:,ix].toarray().astype(np.float32)
        if scale is None:scale=StandardScaler().fit(raw)
        m=copy.deepcopy(bestmodel);m.set_params(**params['params'])
        if 'random_state' in m.get_params():m.set_params(random_state=ctx.cfg['seed'])
        LOG.info('%s fitting full-development %s: %d samples × %d features',ctx.lineage,name,len(meta),len(ix))
        fitted=fit_checked(m,scale.transform(raw),y)
        model_id=f'{ctx.lineage}_{name}_K{ctx.K}_{ctx.signature[:12]}'
        b=make_bundle(ctx,fitted,scale,ix,ds['schema'],model_id)
        b['training_ids']=meta.sample_id.tolist();b['training_groups']=sorted(set(meta.group_id));b['params']=params['params']
        path=ctx.path('models/'+('EN_full_train_original_features.joblib' if name=='original_EN' else 'EN_full_train_stable_features.joblib'))
        dump(path,b);outputs.append(path);loaded=load_final_clock(path)
        before=bundle_predict(b,raw);np.testing.assert_allclose(before,bundle_predict(loaded,raw),atol=1e-10,rtol=1e-10)
        s=score(y,before);metric_rows.append(dict(model_version=name,evaluation_source='full_development_in_sample',fold='all',K=ctx.K,n_features=len(ix),**s))
        q=meta.copy();q['prediction']=before;q['model_version']=name;q['evaluation_source']='full_development_in_sample';development.append(q)
        model_specs[name]=dict(path=str(path.resolve()),model_id=model_id,n_features=len(ix),solver=fitted.__class__.__module__,
            alpha=float(fitted.alpha),l1_ratio=float(fitted.l1_ratio),dual_gap=float(fitted.dual_gap_),training_samples=len(meta),training_groups=meta.group_id.nunique())
    reference=make_bundle(ctx,bestmodel,bestpre.scaler,bestpre.indices,ds['schema'],f'{ctx.lineage}_EN_best_fold_reference_K{ctx.K}_fold{int(best.fold)}')
    reference['training_scope']='best_fold_reference';reference['source_path']=str(bestpath.resolve())
    ref=ctx.path('models/EN_best_fold_reference.joblib');dump(ref,reference);outputs.append(ref)
    final=load_final_clock(model_specs[selected]['path']);path=ctx.path('models/final_clock.joblib');dump(path,final);outputs.append(path)
    path=ctx.path('models/final_preprocessor.joblib');dump(path,final['preprocessor']);outputs.append(path)
    path=ctx.path('models/final_feature_schema.json');write_json(path,final['feature_schema']);outputs.append(path)
    decision=dict(selected_version=selected,selection_basis='ONLY paired outer-valid original/stable EN at locked internal-CV K',
        original_mean_R2=float(original.R2.mean()),stable_mean_R2=float(primary.R2.mean()),original_mean_MAE=float(original.MAE.mean()),
        stable_mean_MAE=float(primary.MAE.mean()),improved_folds=int(improved.sum()),K=ctx.K,primary_threshold=ctx.options['threshold'],
        best_reference_fold=int(best.fold),full_train_hyperparameters=params['params'],feature_rank_source=rank_source,
        stable_features_role='predictive and explanatory' if promote else 'explanatory; reduced model predicts materially worse in internal CV',
        versions=model_specs,final_clock_path=str(ctx.path('models/final_clock.joblib').resolve()),locked_before_external=now(),
        n_bootstrap=ctx.options['n_bootstrap'],n_stable_features=len(stable),selection_affected_K_CV=True)
    p=ctx.path('models/model_versions.json');write_json(p,decision);outputs.append(p)
    p=ctx.path('models/final_training_manifest.json');write_json(p,dict(**decision,training_ids=meta.sample_id.tolist(),training_group_ids=sorted(set(meta.group_id)),
        independent_donor_verified=False,source_signature=ctx.source_signature,age_unit='years',target='chronological age, not independently validated biological age',
        external_labels_used_for_primary_model=False,external_used_for_model_selection=False,source_tuning_path=str(ctx.path(f'checkpoints/tuning/full_train/K{ctx.K}_ElasticNet.json'))));outputs.append(p)
    save('stability/final_clock_development_predictions.csv',pd.concat(development,ignore_index=True))
    save('stability/final_clock_model_comparison_metrics.csv',pd.DataFrame(metric_rows))
    p=ctx.path('models/model_decision_report.md')
    write_text(p,f'''# Final ElasticNet decision: {ctx.lineage}

Final: **{selected}**, full-development K={ctx.K}; alpha={params['params']['alpha']}, l1_ratio={params['params']['l1_ratio']}.
Internal fivefold mean R2: original={original.R2.mean():.6f}, stable={primary.R2.mean():.6f}; MAE: original={original.MAE.mean():.6f}, stable={primary.MAE.mean():.6f}. Improved folds={int(improved.sum())}/5; rule requires at least4 and mean MAE improvement with R2 loss<=0.02. External is not consulted.

Stable selection was independently refitted in each outer TRAIN, as proved by bootstrap draws and preprocessing IDs. K was selected from the existing internal benchmark, so CV is conditional/selection-affected, not a fresh unbiased nested-K experiment. Full-development train scores are explicitly in-sample, not evidence of improvement. Threshold0.75 is sensitivity-only and never replaces0.70 using external performance.

Both full-development EN bundles and a clearly marked best-fold reference are saved. Final clock is a chronological-age predictor. Stable features remain interpretable predictive associations; positive coefficients are not causal pro-aging drivers.
''');outputs.append(p)
    return outputs+pins+[Path(rank_source),ctx.path(f'checkpoints/tuning/full_train/K{ctx.K}_ElasticNet.json')]
