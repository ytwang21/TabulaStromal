"""Export current-data stable EN artifacts; assemble both-lineage measured figures."""
from pathlib import Path
import json,logging
import numpy as np
import pandas as pd
from .storage import write_csv,write_json,dump,atomic_path,process_lock,annotation_identity


def stable_exports(ctx,ds,table,stable,params,K):
    """Current stable set is determined by training bootstrap, never external age."""
    from .stability import fixed_subset_pre
    from .models import fit_model,predict_model
    from .data import load_dataset
    from .metrics import evaluate
    from .external import adaptation_split
    out=ctx.root/'stable_clock_optimization';out.mkdir(exist_ok=True)
    prefix='performance_first';f=stable.copy()
    write_csv(out/f'{prefix}_features.csv',f)
    for field,label in [('Celltype','celltypes'),('Gene','genes')]:
        t=f.groupby(field).agg(n_features=('feature','size'),mean_frequency=('selection_frequency','mean')).reset_index()
        write_csv(out/f'{prefix}_{label}.csv',t)
    write_json(out/'current_protocol.json',dict(signature=ctx.cache.signature,revision=ctx.cfg.get('execution_revision'),
        stable_rule=f'bootstrap_frequency >= {ctx.cfg["stability"]["threshold"]}',n_bootstrap=ctx.cfg['stability']['n_bootstrap'],
        selection_basis='internal CV + training bootstraps only; no external outcome optimization',n_features=len(f),
        n_celltypes=f.Celltype.nunique(),n_genes=f.Gene.nunique(),membership='unique projection of stable pairs; normalized coverage in stability/stable_Celltype.csv',
        previous_external_selected_results='obsolete after annotation change; do not load files outside current-signature stage manifests'))
    if not len(stable):
        write_json(out/'stable_model_unavailable.json',dict(reason='no pair meets prespecified stability threshold',signature=ctx.cache.signature));return []
    rows=np.arange(len(ds['meta']));indices=np.flatnonzero(table.selection_frequency.to_numpy()>=ctx.cfg['stability']['threshold'])
    pre=fixed_subset_pre(ds,rows,indices,K)
    model=fit_model('ElasticNet',pre.transform(ds['X']),ds['meta'].Age.to_numpy(float),ctx.cfg,ctx.cfg['seed'],params)
    bundle=dict(model=model,preprocessor=pre,selected_features=pre.names,stable_feature_table=f,
                schema=ds['schema'],signature=ctx.cache.signature,selection_basis='training bootstrap; external-independent feature selection')
    modelpath=ctx.root/'models/final_stable_EN_bundle.joblib';dump(modelpath,bundle);dump(out/f'{prefix}_source_bundle.joblib',bundle)
    f['source_model_coefficient']=model.coef_;outputs=[modelpath];metrics=[];predictions=[]
    from .external import selected_external_items
    for it in selected_external_items(ctx):
        ext=load_dataset(ctx.root/'features'/it['dataset']);meta=ext['meta'];y=meta.Age.to_numpy(float)
        if not len(meta):continue
        X=pre.transform(ext['X']);pred=model.predict(X)
        metrics.append(dict(dataset=it['dataset'],scope='no_adaptation_ALL100',strategy='no_adaptation',**evaluate(y,pred)))
        t=meta.copy();t['prediction']=pred;t['dataset']=it['dataset'];t['scope']='no_adaptation_ALL100';predictions.append(t)
        if it['dataset']=='combined':
            split,reason=adaptation_split(meta,.3,ctx.cfg['seed'])
            if split is None:
                write_json(out/'FT30_skipped.json',dict(reason=reason));continue
            ai,ti=split
            # Preserve the reference EN .fit() adaptation and source scaler/params.
            adapted=fit_model('ElasticNet',X[ai],y[ai],ctx.cfg,ctx.cfg['seed'],params)
            ab={**bundle,'model':adapted,'adaptation_ids':meta.sample_id.iloc[ai].tolist(),
                'strategy':'refit','adaptation_detail':'fixed source scaler/EN parameters; fit on external adaptation30 only',
                'independent_evaluation_ids':meta.sample_id.iloc[ti].tolist()}
            adaptedpath=ctx.root/'models/final_stable_EN_FT30_bundle.joblib';dump(adaptedpath,ab);dump(out/f'{prefix}_FT30_bundle.joblib',ab);outputs.append(adaptedpath)
            after=adapted.predict(X);f['adapted_coefficient']=adapted.coef_
            q=meta.copy();q['prediction']=after;q['baseline_prediction']=pred;q['used_for_adaptation']=False;q.loc[ai,'used_for_adaptation']=True
            q['scope']='ALL100_post_adaptation_EXPLORATORY';q['dataset']='combined'
            write_csv(out/f'{prefix}_all100_predictions.csv',q);predictions.append(q)
            write_csv(out/f'{prefix}_heldout70_predictions.csv',q.iloc[ti])
            for scope,ix in [('ALL100_post_adaptation_EXPLORATORY',np.arange(len(meta))),('heldout70_adapted',ti),('heldout70_baseline',ti)]:
                pp=pred if scope.endswith('baseline') else after
                metrics.append(dict(dataset='combined',scope=scope,strategy='refit' if pp is after else 'no_adaptation',**evaluate(y[ix],pp[ix])))
    write_csv(out/f'{prefix}_features.csv',f)
    for field,label in [('Celltype','celltypes'),('Gene','genes')]:
        if 'adapted_coefficient' in f:
            t=f.assign(nonzero=f.adapted_coefficient.ne(0)).groupby(field).agg(n_features=('feature','size'),n_adapted_nonzero=('nonzero','sum'),mean_frequency=('selection_frequency','mean')).reset_index()
            write_csv(out/f'{prefix}_{label}.csv',t)
    write_csv(out/f'{prefix}_metrics.csv',pd.DataFrame(metrics));write_csv(ctx.root/'external/final_stable_EN_predictions.csv',pd.concat(predictions,ignore_index=True))
    return outputs+[out/f'{prefix}_metrics.csv',out/f'{prefix}_features.csv',out/f'{prefix}_celltypes.csv',out/f'{prefix}_genes.csv',out/'current_protocol.json']


def combined_figures(ctx):
    """Last finishing lineage renders shared outputs under one OS lock."""
    from .figures import style,plot_benchmark,scatter,dotplot,read
    from .storage import file_info
    import matplotlib.pyplot as plt
    base=Path(ctx.cfg['output_root']);roots={lin:Path(v.get('output_dir',str(base/lin))) for lin,v in ctx.cfg['lineages'].items()}
    other='Endo' if ctx.lineage=='Fibro' else 'Fibro';identity=roots[other]/'reports/run_identity.json';status=roots[other]/'reports/run_status.json'
    if not identity.exists() or not status.exists():return []
    ident=json.loads(identity.read_text());state=json.loads(status.read_text())
    if ident['config'].get('execution_revision')!=ctx.cfg.get('execution_revision') or not state['stages'].get('figures',{}).get('status','').startswith('complete'):return []
    code=json.loads((ctx.root/'reports/run_identity.json').read_text())['code']
    if ident['code']!=code:return []
    for item in ident['inputs']:
        if 'annotation_sha256' in item:
            p=item['path']
            if file_info(p)['mtime_ns']!=item['mtime_ns'] or annotation_identity(p,ctx.cfg)!=item['annotation_sha256']:return []
    dest=base/'Stromal_Combined_Figure/outputs';dest.mkdir(parents=True,exist_ok=True);style();outputs=[]
    with process_lock(dest/'combined_render.lock',blocking=True):
        signatures={lin:json.loads((root/'reports/run_identity.json').read_text())['signature'] for lin,root in roots.items()}
        manifest=dest/'current_14models_figure_manifest.json'
        if manifest.exists():
            previous=json.loads(manifest.read_text())
            if previous.get('lineage_signatures')==signatures and all(Path(p).exists() for p in previous['paths']):return []
        for fold in [None]+list(range(1,ctx.cfg['cv']['n_splits']+1)):
            fig,axes=plt.subplots(6,2,figsize=(24,30));summaries=[]
            for col,lin in enumerate(['Fibro','Endo']):
                root=roots[lin];rank=read(root/'benchmark/model_comparison.csv');cv=read(root/'benchmark/cv_fold_metrics.csv');ex=read(root/'external/external_combined_no_adapt.csv')
                selected=rank.sort_values('internal_rank').groupby('model',sort=False).head(1)[['model','K']]
                bm=cv.merge(selected,on=['model','K']);ep=ex.merge(selected,on=['model','K']);ep['split']='External-valid';bm=pd.concat([bm,ep],ignore_index=True)
                if fold is not None:bm=bm[bm.fold==fold]
                summaries.append(bm.assign(lineage=lin))
                for r,metric in enumerate(['R2','MAE']):
                    plot_benchmark(axes[r,col],bm,metric,ctx.cfg['models']['names'],fold is None)
                    axes[r,col].set_title(f'{chr(65+r*2+col)}. {lin}: {metric}',loc='left')
                    for container in axes[r,col].containers:
                        if hasattr(container,'datavalues'):axes[r,col].bar_label(container,fmt='%.2f',rotation=90,fontsize=4,padding=2)
                en=cv[(cv.model=='ElasticNet')&(cv.split=='Valid')]
                if fold is not None:en=en[en.fold==fold]
                k=int(en.groupby('K').R2.mean().idxmax());en=en[en.K==k]
                chosen=int(en.sort_values('R2',ascending=False).iloc[0].fold) if fold is None else fold
                pred=read(root/'benchmark/cv_predictions.csv');pred=pred[(pred.model=='ElasticNet')&(pred.K==k)&(pred.fold==chosen)]
                ext=read(root/'external/external_combined_no_adapt_predictions.csv');ext=ext[(ext.model=='ElasticNet')&(ext.K==k)&(ext.fold==chosen)].copy();ext['split']='External-valid'
                scatter(axes[2,col],pd.concat([pred,ext],ignore_index=True),f'{chr(69+col)}. {lin} EN K={k}, fold={chosen}; no adaptation')
                ct=read(root/'stability/stable_Celltype.csv');ct=ct[ct.raw_stable_features>0].sort_values('stable_fraction_of_candidates').tail(50)
                axes[3,col].barh(ct.Celltype,ct.stable_fraction_of_candidates,color='#8ebcd9');axes[3,col].tick_params(axis='y',labelsize=6);axes[3,col].set_title(f'{chr(71+col)}. Stable subtype coverage',loc='left')
                genes=read(root/'stability/stable_Gene.csv');genes=genes[genes.raw_stable_pairs>0].sort_values('max_frequency').tail(50)
                axes[4,col].barh(genes.Gene,genes.max_frequency,color='#b9d9ee');axes[4,col].tick_params(axis='y',labelsize=5);axes[4,col].set_title(f'{chr(73+col)}. Stable genes (up to 50)',loc='left')
                enrich=read(root/'enrichment/pathway_enrichment.csv')
                if len(enrich):dotplot(axes[5,col],enrich[enrich.scope=='global'],f'{chr(75+col)}. Measured enrichment')
                else:axes[5,col].text(.5,.5,'No completed pathway enrichment; see limitations.json',ha='center',transform=axes[5,col].transAxes)
            axes[0,0].legend(fontsize=7,frameon=False);fig.suptitle('TabulaStromal: current Celltype annotations; 14 models\nExternal-valid = combined, no adaptation; internal CV selects K/fold',fontsize=15)
            fig.tight_layout(rect=[0,0,1,.97]);name='Stromal_individual_aging_clock_Fibroblast_Endothelial_'+('FINAL_v2' if fold is None else str(fold))
            for ext in ['png','pdf']:
                path=dest/f'{name}.{ext}'
                with atomic_path(path) as tmp:fig.savefig(tmp,format=ext,dpi=ctx.cfg['figures']['dpi'],bbox_inches='tight')
                outputs.append(path)
            plt.close(fig);write_csv(dest/('benchmark_fivefold_summary.csv' if fold is None else f'benchmark_fold_{fold}.csv'),pd.concat(summaries,ignore_index=True))
        write_csv(dest/'feature_K_comparison.csv',pd.concat([read(root/'benchmark/feature_k_comparison.csv').assign(lineage=lin) for lin,root in roots.items()],ignore_index=True))
        write_json(dest/'current_14models_figure_manifest.json',dict(revision=ctx.cfg.get('execution_revision'),paths=[str(p) for p in outputs],lineage_signatures={lin:json.loads((root/'reports/run_identity.json').read_text())['signature'] for lin,root in roots.items()}))
    return outputs
