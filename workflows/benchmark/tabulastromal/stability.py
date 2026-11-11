"""Group bootstraps, fold-internal stable EN comparison, final training refit."""
from pathlib import Path
import gc,itertools,json,logging
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.preprocessing import StandardScaler
from .storage import dump,write_csv,write_json
from .models import Preprocessor,rank_features,fit_model,predict_model,bags_for,load_bundle,save_bundle
from .training import grouped_splits,cached_pre,tune_linear,task_path
from .data import load_dataset
from .metrics import evaluate,METRICS
LOG=logging.getLogger(__name__)


def bootstrap_en(ctx,ds,rows,K,params,tag):
    root=ctx.root/'stability'/tag;root.mkdir(parents=True,exist_ok=True)
    meta=ds['meta'].iloc[rows];groups=meta.group_id.unique();B=ctx.cfg['stability']['n_bootstrap'];p=ds['X'].shape[1]
    coefficients=[];exposures=[]
    for b in range(B):
        path=root/'bootstrap'/f'boot_{b:04d}.joblib';marker=path.with_suffix('.complete.json')
        if ctx.cache.valid(marker):record=joblib.load(path)
        else:
            rng=np.random.default_rng(ctx.cfg['seed']+b);draw=rng.choice(groups,len(groups),replace=True)
            # Duplicate all samples of each selected group, preserving donor/sample grouping.
            ix=np.concatenate([rows[np.flatnonzero(meta.group_id.to_numpy()==g)] for g in draw])
            LOG.info('BOOTSTRAP %s %d/%d, %d sample rows',tag,b+1,B,len(ix))
            rank=rank_features(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],ctx.cfg)
            pre=Preprocessor.fit(ds['X'][ix],ds['meta'].Age.to_numpy(float)[ix],ds['schema']['features'],ds['meta'].sample_id.iloc[ix],K,ctx.cfg,rank)
            model=fit_model('ElasticNet',pre.transform(ds['X'][ix]),ds['meta'].Age.to_numpy(float)[ix],ctx.cfg,ctx.cfg['seed']+b,params)
            record=dict(indices=pre.indices,coefficient=np.asarray(model.coef_,np.float32),
                        input_coefficient=np.asarray(model.coef_/pre.scaler.scale_,np.float32),
                        bootstrap=b,draw_groups=draw,fit_ids=pre.fit_ids,params=params)
            dump(path,record);ctx.cache.mark(marker,[path]);del pre,model
        beta=np.where(np.abs(record['coefficient'])>ctx.cfg['stability']['coefficient_epsilon'],record['coefficient'],0)
        coefficients.append(sp.csr_matrix((beta,(np.zeros(len(beta),int),record['indices'])),shape=(1,p)))
        exposures.append(sp.csr_matrix((np.ones(len(beta)),(np.zeros(len(beta),int),record['indices'])),shape=(1,p)))
    mat=sp.vstack(coefficients,format='csc');exp=sp.vstack(exposures,format='csc')
    nz=np.asarray((mat!=0).sum(0)).ravel();positive=np.asarray((mat>0).sum(0)).ravel();negative=np.asarray((mat<0).sum(0)).ravel()
    exposure=np.asarray(exp.sum(0)).ravel();freq=nz/B
    lo=np.zeros(p);hi=np.zeros(p);median=np.zeros(p);mean=np.asarray(mat.mean(0)).ravel()
    variable=np.flatnonzero(nz)
    for off in range(0,len(variable),512):
        idx=variable[off:off+512];arr=mat[:,idx].toarray();lo[idx],median[idx],hi[idx]=np.quantile(arr,[.025,.5,.975],axis=0)
    ng=len(ds['schema']['genes']);ct=np.repeat(ds['schema']['celltypes'],ng);genes=np.tile(ds['schema']['genes'],len(ds['schema']['celltypes']))
    # Opportunity universe is train expression with nonzero sample-to-sample variance.
    X=ds['X'][rows];means=np.asarray(X.mean(0)).ravel();variance=np.asarray(X.power(2).mean(0)).ravel()-means**2
    table=pd.DataFrame(dict(feature=ds['schema']['features'],Celltype=ct,Gene=genes,selection_frequency=freq,
              selected_count=nz,exposure_count=exposure,conditional_selection_frequency=np.divide(nz,exposure,out=np.zeros(p),where=exposure>0),
              positive_count=positive,negative_count=negative,sign_consistency=np.divide(np.maximum(positive,negative),nz,out=np.zeros(p),where=nz>0),
              coefficient_mean=mean,coefficient_median=median,coefficient_q025=lo,coefficient_q975=hi,candidate=variance>ctx.cfg['features']['variance_threshold']))
    write_csv(root/'bootstrap_all_features.csv',table)
    write_json(root/'bootstrap_definition.json',dict(n_bootstrap=B,fit_rows=ds['meta'].sample_id.iloc[rows].tolist(),
                  definition='nonzero coefficient count / total bootstraps; filtering, MI selection and scaling refit in each resample',
                  interval='2.5/97.5 coefficient percentiles including zero for unselected features; descriptive post-selection interval'))
    return table


def fixed_subset_pre(ds,rows,indices,K):
    scaler=StandardScaler().fit(ds['X'][rows][:,indices].toarray())
    return Preprocessor(indices,np.asarray(ds['schema']['features'])[indices],scaler,ds['meta'].sample_id.iloc[rows],K)


def stability(ctx):
    ds=load_dataset(ctx.root/'features/train');meta=ds['meta'];y=meta.Age.to_numpy(float)
    ranking=pd.read_csv(ctx.root/'benchmark/model_comparison.csv')
    if not len(ranking):raise RuntimeError('No benchmark ranking')
    winner=ranking.iloc[0];en=ranking[ranking.model=='ElasticNet']
    if en.empty:raise RuntimeError('ElasticNet missing from completed benchmark')
    en=en.iloc[0];K=int(en.K);threshold=ctx.cfg['stability']['threshold'];comparisons=[]
    splits=list(grouped_splits(meta,ctx.cfg['cv']['n_splits'],ctx.cfg['seed'],ctx.cfg['cv']['age_bins']))
    for fold,(tr,va) in enumerate(splits,1):
        path=ctx.root/'stability'/f'outer{fold}';marker=path/'comparison.complete.json'
        if ctx.cache.valid(marker):comparisons.extend(pd.read_csv(path/'comparison.csv').to_dict('records'));continue
        _,base,pre,info=load_bundle(task_path(ctx,K,'ElasticNet',fold),ctx.cfg)
        table=bootstrap_en(ctx,ds,tr,K,info['params'],f'outer{fold}')
        original=evaluate(y[va],base.predict(pre.transform(ds['X'][va])));foldrows=[]
        for th in [threshold]+ctx.cfg['stability']['sensitivity_thresholds']:
            indices=np.flatnonzero(table.selection_frequency.to_numpy()>=th)
            record=dict(fold=fold,K=K,threshold=th,n_stable=len(indices),n_subtypes=table.iloc[indices].Celltype.nunique(),
                        comparison='new_benchmark_full_feature_EN versus EN_stable_refit',baseline_label='best_original_EN_within_this_rebenchmark')
            record.update({'original_'+m:original[m] for m in METRICS})
            if len(indices):
                subpre=fixed_subset_pre(ds,tr,indices,K);model=fit_model('ElasticNet',subpre.transform(ds['X'][tr]),y[tr],ctx.cfg,ctx.cfg['seed']+fold,info['params'])
                stable=evaluate(y[va],model.predict(subpre.transform(ds['X'][va])));record.update({'stable_'+m:stable[m] for m in METRICS})
                record.update({'delta_'+m:stable[m]-original[m] for m in METRICS})
                record['improved']=stable['MAE']<original['MAE'] and stable['R2']>=original['R2']-ctx.cfg['stability']['r2_tolerance']
                save_bundle(path/f'threshold_{th:.2f}','ElasticNet',model,subpre,dict(lineage=ctx.lineage,K=K,fold=fold,params=info['params'],train_ids=meta.sample_id.iloc[tr].tolist(),valid_ids=meta.sample_id.iloc[va].tolist()))
            else:record.update(improved=False,reason='no feature reaches prespecified stability threshold')
            foldrows.append(record)
        write_csv(path/'comparison.csv',pd.DataFrame(foldrows));ctx.cache.mark(marker,[path/'comparison.csv']);comparisons.extend(foldrows)
    compare=pd.DataFrame(comparisons);write_csv(ctx.root/'stability/stable_vs_original_EN_cv.csv',compare)
    primary=compare[np.isclose(compare.threshold,threshold)]
    promote=(int(primary.improved.sum())>=ctx.cfg['stability']['min_fold_improvements'] and
             'stable_MAE' in primary and primary.stable_MAE.notna().all() and primary.stable_MAE.mean()<primary.original_MAE.mean() and
             primary.stable_R2.mean()>=primary.original_R2.mean()-ctx.cfg['stability']['r2_tolerance'])
    allrows=np.arange(len(meta));params=tune_linear(ctx,ds,allrows,K,'ElasticNet','full_train')
    table=bootstrap_en(ctx,ds,allrows,K,params,'full_train');stable=table[table.selection_frequency>=threshold].copy()
    write_csv(ctx.root/'stability/stable_Celltype_Gene.csv',stable)
    for th in ctx.cfg['stability']['sensitivity_thresholds']:write_csv(ctx.root/'stability'/f'stable_features_threshold_{th:.2f}.csv',table[table.selection_frequency>=th])
    cell=table.groupby('Celltype').agg(candidate_features=('candidate','sum'),raw_stable_features=('selection_frequency',lambda s:int((s>=threshold).sum())),
                    total_selection_frequency=('selection_frequency','sum'),mean_selection_frequency=('selection_frequency','mean')).reset_index()
    cell['stable_fraction_of_candidates']=cell.raw_stable_features/np.maximum(cell.candidate_features,1)
    coverage=pd.read_csv(ds['path']/'celltype_coverage.csv');cell=cell.merge(coverage,on='Celltype')
    cell['stable_features_per_1000_cells']=1000*cell.raw_stable_features/np.maximum(cell.n_cells,1)
    cell['stable_features_per_observed_sample']=cell.raw_stable_features/np.maximum(cell.n_samples,1)
    write_csv(ctx.root/'stability/stable_Celltype.csv',cell)
    gene=table.groupby('Gene').agg(raw_stable_pairs=('selection_frequency',lambda s:int((s>=threshold).sum())),max_frequency=('selection_frequency','max'),
                 mean_frequency=('selection_frequency','mean'),candidate_celltypes=('candidate','sum')).reset_index()
    gene['stable_fraction_of_candidate_celltypes']=gene.raw_stable_pairs/np.maximum(gene.candidate_celltypes,1)
    write_csv(ctx.root/'stability/stable_Gene.csv',gene);write_csv(ctx.root/'stability/stable_celltype_specific_genes.csv',stable[['Celltype','Gene','selection_frequency','sign_consistency']])
    name=str(winner.model);finalK=int(winner.K);scheme='full_feature_refit'
    if name=='ElasticNet' and promote and finalK==K and len(stable):scheme='stable_EN_refit';indices=np.flatnonzero(table.selection_frequency>=threshold);finalpre=fixed_subset_pre(ds,allrows,indices,K)
    else:finalpre=cached_pre(ctx,ds,allrows,finalK,'full_train_final')
    finalparams=tune_linear(ctx,ds,allrows,finalK,name,'full_train_final')
    finalpath=ctx.root/'models/final_model';marker=finalpath/'complete.json'
    if not ctx.cache.valid(marker):
        bags=bags_for(ds,finalpre,allrows) if name in ['DeepSets','Ensemble'] else None
        epoch=finalpath/'epochs'/ctx.cache.signature/'last.pt';epoch.parent.mkdir(parents=True,exist_ok=True)
        model=fit_model(name,finalpre.transform(ds['X']),y,ctx.cfg,ctx.cfg['seed'],finalparams,bags,epoch,ctx.resume)
        outputs=save_bundle(finalpath,name,model,finalpre,dict(lineage=ctx.lineage,K=finalK,scheme=scheme,params=finalparams,seed=ctx.cfg['seed'],train_ids=meta.sample_id.tolist(),valid_ids=[],selection_basis='internal CV only'))
        finalschema=dict(ds['schema'],selected_features=finalpre.names,scheme=scheme,model=name,K=finalK,signature=ctx.cache.signature)
        write_json(finalpath/'feature_schema.json',finalschema)
        write_json(finalpath/'training_config.json',ctx.cfg)
        from .storage import write_text
        write_text(finalpath/'environment.txt',(ctx.root/'audit/environment.txt').read_text())
        ctx.cache.mark(marker,outputs+[finalpath/'feature_schema.json',finalpath/'training_config.json',finalpath/'environment.txt'])
    _,model,finalpre,info=load_bundle(finalpath,ctx.cfg)
    final_external=[];final_predictions=[]
    for it in json.loads((ctx.root/'audit/inventory.json').read_text())[1:]:
        ext=load_dataset(ctx.root/'features'/it['dataset'])
        if not len(ext['meta']):continue
        ix=np.arange(len(ext['meta']));bags=bags_for(ext,finalpre,ix) if name in ['DeepSets','Ensemble'] else None
        p=predict_model(model,name,finalpre.transform(ext['X']),bags)
        record=dict(lineage=ctx.lineage,model=name,K=finalK,scheme=scheme,dataset=it['dataset'],method=it['method'],role=it['role'],selection='fixed internal-CV-selected full-training refit',**evaluate(ext['meta'].Age,p))
        final_external.append(record);t=ext['meta'].copy();t['prediction']=p;t['dataset']=it['dataset'];t['method']=it['method'];t['model']=name;t['K']=finalK;final_predictions.append(t)
        del ext;gc.collect()
    write_csv(ctx.root/'external/final_refit_external_metrics.csv',pd.DataFrame(final_external));write_csv(ctx.root/'external/final_refit_external_predictions.csv',pd.concat(final_predictions,ignore_index=True) if final_predictions else pd.DataFrame(columns=['sample_id','Age','prediction','dataset']))
    write_json(ctx.root/'stability/selection_decision.json',dict(internal_winner=name,K=finalK,EN_internal_rank=int(en.internal_rank),
             EN_K=K,EN_stable_eligible=bool(promote),chosen_scheme=scheme,n_stable_pairs=len(stable),n_stable_genes=stable.Gene.nunique(),
             n_stable_subtypes=stable.Celltype.nunique(),interpretation_label='selected clock features' if name=='ElasticNet' else 'exploratory EN stability; selected clock is '+name,
             archived_baseline_note='best_original_EN labels the un-reduced EN of this fresh benchmark, not a reconstructed historical pickle; archived results remain separate and untouched'))
    write_text_report(ctx,table,stable,scheme)
    return [ctx.root/'stability'/f for f in ['stable_Celltype_Gene.csv','stable_Celltype.csv','stable_Gene.csv','stable_vs_original_EN_cv.csv','selection_decision.json']]+[marker]


def write_text_report(ctx,table,stable,scheme):
    from .storage import write_text
    write_text(ctx.root/'stability/stability_methods.md',f'''# Stability definitions

Primary frequency threshold: {ctx.cfg['stability']['threshold']}; bootstrap resamples: {ctx.cfg['stability']['n_bootstrap']} per outer train and full training set.
Stable pairs: {len(stable)}, genes: {stable.Gene.nunique()}, subtypes: {stable.Celltype.nunique()}.
Final scheme: {scheme}. No external performance enters this decision.
Each bootstrap repeats variance filtering, MI ranking, topK selection, scaling and EN fitting. Full-training source hyperparameters are fixed within bootstraps.
Sign consistency is max(n_positive,n_negative)/n_nonzero; intervals include coefficient zero in unselected bootstraps. Full coefficient distributions and sample draws are in bootstrap/*.joblib.
Celltype stable fraction divides stable pairs by variance-eligible candidates. Coverage, cells and sample-normalized rates are also shown as descriptive bias diagnostics, never a quota to force subtype inclusion.
Gene stability reports max and mean pair frequency, raw stable-pair count and candidate-celltype-normalized coverage; none is a new bootstrap gene-level selection probability.
The selected-K/selected-model CV scores are selection-affected estimates. External fixed refit evaluation is the independent confirmation; A/B and alternate integration methods share samples.
''')
