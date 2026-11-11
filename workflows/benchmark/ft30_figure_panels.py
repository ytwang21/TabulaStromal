"""Execution-only adapter of the original measured figure panels.

Original package remains unchanged to retain trained-model cache signatures.
Selection differs only in using the explicit best INTERNAL Valid fold; external
benchmark errors use MedAE while age-label MAD remains in CSV and supplement.
Per-GSE plotting is not scheduled by the caller; source results remain untouched.
"""
from tabulastromal.figures import *

def figures_from_saved(ctx):
    style();outputs=[];ranking=read(ctx.root/'benchmark/model_comparison.csv');folds=read(ctx.root/'benchmark/cv_fold_metrics.csv');external=read(ctx.root/'external/external_combined_no_adapt.csv')
    if ranking.empty or folds.empty:raise RuntimeError('Benchmark results required for figures')
    # Best K for each model is chosen using internal validation only.
    selected=ranking.sort_values('internal_rank').groupby('model',sort=False).head(1)[['model','K']]
    bm=folds.merge(selected,on=['model','K']);ext=external.merge(selected,on=['model','K']) if len(external) else pd.DataFrame()
    if len(ext):ext['split']='External-valid';bm=pd.concat([bm,ext],ignore_index=True)
    models=selected.model.tolist();representatives=[]
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for ax,metric in zip(axes,['R2','MAE']):
        for name,sub in folds[folds.split=='Valid'].groupby('model'):
            stats=sub.groupby('K')[metric].agg(['mean','std']);ax.errorbar(stats.index,stats['mean'],yerr=stats['std'].fillna(0),marker='.',lw=.8,label=name)
        ax.set_xlabel('Requested total Celltype_Gene K');ax.set_ylabel('Validation '+metric)
    axes[-1].legend(fontsize=5,frameon=False);fig.tight_layout();outputs+=save(ctx,fig,'feature_K_comparison')
    for r in selected.itertuples():
        candidates=folds[(folds.model==r.model)&(folds.K==r.K)&(folds.split=='Valid')].copy()
        row=candidates.sort_values(['R2','MAE','fold'],ascending=[False,True,True],na_position='last').iloc[0]
        representatives.append(dict(model=r.model,K=int(r.K),fold=int(row.fold),strategy='highest internal Valid R2; lowest MAE; lowest fold'))
    reps=pd.DataFrame(representatives);write_csv(ctx.root/'figures/internal_representative_folds.csv',reps)
    for label,table in [('fivefold_summary',bm),('internal_representative_fold',bm.merge(reps[['model','K','fold']],on=['model','K','fold']))]:
        fig,axes=plt.subplots(2,2,figsize=(17,10))
        for ax,metric in zip(axes.ravel(),['R2','MAE','MedAE','PCC']):plot_benchmark(ax,table,metric,models,label=='fivefold_summary')
        axes[0,0].legend(frameon=False,ncol=3,fontsize=8)
        fig.suptitle(f'{ctx.lineage}: {label.replace("_"," ")}\nExternal-valid = external combined, no adaptation; K and representative fold selected by internal CV',fontsize=12)
        fig.text(.02,.005,'MAD = mean absolute age-label deviation; it cannot rank models within the same evaluation sample set. Bars: fold means; error bars: descriptive fold SD.',fontsize=8)
        fig.tight_layout(rect=[0,.025,1,.94]);outputs+=save(ctx,fig,'benchmark_'+label)
    winner=ranking.iloc[0];K=int(winner.K);model=str(winner.model);rep=int(reps[reps.model==model].fold.iloc[0])
    cvp=read(ctx.root/'benchmark/cv_predictions.csv');cvp=cvp[(cvp.K==K)&(cvp.model==model)&(cvp.fold==rep)]
    ep=read(ctx.root/'external/external_combined_no_adapt_predictions.csv')
    if len(ep):ep=ep[(ep.K==K)&(ep.model==model)&(ep.fold==rep)].copy();ep['split']='External-valid';cvp=pd.concat([cvp,ep],ignore_index=True)
    fig,ax=plt.subplots(figsize=(8,7));scatter(ax,cvp,f'{ctx.lineage} {model}, K={K}, internal representative fold {rep}');outputs+=save(ctx,fig,'prediction_scatter')
    adapted=read(ctx.root/'external/external_combined_adaptation_predictions.csv')
    if len(adapted):
        sub=adapted[(adapted.K==K)&(adapted.model==model)&(adapted.fold==rep)&np.isclose(adapted.nominal_adaptation_fraction,.3)&(adapted['repeat']==0)]
        if len(sub):
            before=sub.copy();before['prediction']=before.baseline_prediction;before['split']='Baseline all100'
            after=sub.copy();after['split']='Updated all100'
            fig,ax=plt.subplots(figsize=(8,7));scatter(ax,pd.concat([before,after]),'Paired all100 evaluation: includes 30% adaptation samples; not independent held-out');outputs+=save(ctx,fig,'adaptation_paired_scatter')
    gse=read(ctx.root/'external/external_per_gse_no_adapt.csv')
    if len(gse):
        gs=gse[(gse.K==K)&(gse.model==model)&(gse.fold==rep)&(gse.method==ctx.cfg['external']['primary_method'])]
        if len(gs):
            fig,axes=plt.subplots(1,3,figsize=(16,5))
            for ax,metric in zip(axes,['R2','MAE','MAD']):
                ax.bar(np.arange(len(gs)),gs[metric],color='#b9d9ee');ax.set_xticks(np.arange(len(gs)));ax.set_xticklabels([f'{r.GSE}\nn={r.n_test}' for r in gs.itertuples()],rotation=65,ha='right',fontsize=6);ax.set_ylabel(metric)
            fig.suptitle('Per-GSE, no adaptation; internal representative fold; undefined metrics remain NaN');fig.tight_layout();outputs+=save(ctx,fig,'per_GSE_no_adaptation')
    ada=read(ctx.root/'external/external_per_gse_adaptation.csv')
    if len(ada):
        ada=ada[(ada.K==K)&(ada.model==model)&(ada.fold==rep)&(ada.method==ctx.cfg['external']['primary_method'])]
        if len(ada):
            fig,axes=plt.subplots(1,3,figsize=(17,5))
            for ax,metric in zip(axes,['R2','MAE','MAD']):
                for label,sub in ada.groupby('GSE'):
                    t=sub.groupby('adaptation_fraction')[metric].agg(['mean','std']);ax.plot(t.index,t['mean'],marker='.',lw=.8,label=label)
                ax.set_xlabel('Adaptation fraction');ax.set_ylabel('Held-out '+metric)
            axes[-1].legend(fontsize=5,frameon=False);fig.suptitle('Per-GSE adaptation; held-out subsets vary by fraction; compare paired deltas');fig.tight_layout();outputs+=save(ctx,fig,'per_GSE_adaptation')
            fig,ax=plt.subplots(figsize=(8,5))
            for label,sub in ada.groupby('GSE'):
                t=sub.groupby('adaptation_fraction').delta_MAE.agg(['mean','std']);ax.errorbar(t.index,t['mean'],yerr=t['std'].fillna(0),marker='.',lw=.8,label=label)
            ax.axhline(0,color='gray',ls='--',lw=.7);ax.set_xlabel('Adaptation fraction');ax.set_ylabel('Paired ΔMAE (adapted − baseline; lower improves)');ax.legend(fontsize=5,frameon=False);outputs+=save(ctx,fig,'paired_improvement_by_fraction')
    stable=read(ctx.root/'stability/stable_Celltype_Gene.csv');allf=read(ctx.root/'stability/full_train/bootstrap_all_features.csv');cell=read(ctx.root/'stability/stable_Celltype.csv');gene=read(ctx.root/'stability/stable_Gene.csv');pathways=read(ctx.root/'enrichment/pathway_enrichment.csv');overlap=read(ctx.root/'enrichment/aging_sets_overlap.csv')
    panel_specs=[]
    for metric in ['R2','MAE','MedAE','PCC']:panel_specs.append((f'Benchmark {metric}',lambda ax,m=metric:plot_benchmark(ax,bm,m,models)))
    panel_specs.append(('Predicted versus chronological age',lambda ax:scatter(ax,cvp,'Internal representative fold')))
    if len(allf):
        def frequency(ax):
            ax.hist(allf.selection_frequency,bins=25,color='#b9d9ee');ax.axvline(ctx.cfg['stability']['threshold'],color='gray',ls='--');ax.set_xlabel('Bootstrap selection frequency');ax.set_ylabel('Celltype_Gene count')
        fig,ax=plt.subplots(figsize=(7,4));frequency(ax);outputs+=save(ctx,fig,'selection_frequency');panel_specs.append(('Selection frequency',frequency))
    if len(stable):
        def coefficients(ax):
            s=stable.assign(strength=stable.coefficient_mean.abs()).nlargest(20,'strength').iloc[::-1];ax.barh(s.feature,s.coefficient_mean,color='#8ebcd9');ax.tick_params(axis='y',labelsize=5);ax.axvline(0,color='gray',lw=.6);ax.set_xlabel('Mean standardized coefficient')
        fig,ax=plt.subplots(figsize=(10,7));coefficients(ax);outputs+=save(ctx,fig,'stable_feature_coefficients');panel_specs.append(('Stable feature direction/strength',coefficients))
        topgenes=stable.assign(strength=stable.coefficient_mean.abs()).groupby('Gene').strength.sum().nlargest(25).index
        contribution=stable[stable.Gene.isin(topgenes)].pivot_table(index='Celltype',columns='Gene',values='coefficient_mean',aggfunc='sum').fillna(0)
        fig,ax=plt.subplots(figsize=(12,max(4,len(contribution)*.25)));bound=max(float(np.abs(contribution.to_numpy()).max()),1e-8)
        im=ax.imshow(contribution.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-bound,vmax=bound);ax.set_xticks(np.arange(len(contribution.columns)));ax.set_xticklabels(contribution.columns,rotation=65,ha='right',fontsize=7);ax.set_yticks(np.arange(len(contribution)));ax.set_yticklabels(contribution.index,fontsize=6);fig.colorbar(im,ax=ax,label='Mean standardized coefficient');fig.tight_layout();outputs+=save(ctx,fig,'subtype_gene_contributions')
    if len(cell):
        def subtypes(ax):
            s=cell[cell.raw_stable_features>0].sort_values('stable_fraction_of_candidates').tail(20);ax.barh(s.Celltype,s.stable_fraction_of_candidates,color='#b9d9ee');ax.set_xlabel('Stable / variance-eligible candidate pairs');ax.tick_params(axis='y',labelsize=5)
        fig,ax=plt.subplots(figsize=(10,6));subtypes(ax);outputs+=save(ctx,fig,'stable_celltypes');panel_specs.append(('Subtype coverage, candidate-normalized',subtypes))
        fig,ax=plt.subplots(figsize=(10,6));s=cell.sort_values('raw_stable_features').tail(20);ax.barh(s.Celltype,s.raw_stable_features,color='#b9d9ee');ax.set_xlabel('Raw stable pair count');fig.tight_layout();outputs+=save(ctx,fig,'stable_celltypes_raw_counts')
    if len(gene) and (gene.raw_stable_pairs>0).any():
        def genes(ax):
            s=gene[gene.raw_stable_pairs>0].sort_values(['max_frequency','raw_stable_pairs']).tail(20);ax.barh(s.Gene,s.max_frequency,color='#8ebcd9');ax.set_xlabel('Maximum pair frequency across subtypes');ax.tick_params(axis='y',labelsize=6)
        fig,ax=plt.subplots(figsize=(7,6));genes(ax);outputs+=save(ctx,fig,'stable_genes');panel_specs.append(('Top stable genes',genes))
    if len(overlap):
        fig,ax=plt.subplots(figsize=(10,7));dotplot(ax,overlap[overlap.scope=='global'],'Provided unverified aging lists');outputs+=save(ctx,fig,'aging_sets_overlap_dotplot')
        matrix=overlap.pivot_table(index='scope',columns='term',values='overlap_count',aggfunc='first').fillna(0)
        fig,ax=plt.subplots(figsize=(11,max(4,len(matrix)*.25)));im=ax.imshow(matrix.to_numpy(),aspect='auto',cmap='Blues');ax.set_xticks(np.arange(len(matrix.columns)));ax.set_xticklabels(matrix.columns,rotation=45,ha='right');ax.set_yticks(np.arange(len(matrix)));ax.set_yticklabels(matrix.index,fontsize=6);fig.colorbar(im,ax=ax,label='Overlap gene count');fig.tight_layout();outputs+=save(ctx,fig,'aging_sets_overlap_heatmap')
        panel_specs.append(('Aging sets (provided_unverified)',lambda ax:dotplot(ax,overlap[overlap.scope=='global'],'Provided aging sets')))
    if len(pathways):
        for library,sub in pathways.groupby('library'):
            fig,ax=plt.subplots(figsize=(11,7));dotplot(ax,sub[sub.scope=='global'],library);fig.tight_layout();outputs+=save(ctx,fig,'enrichment_'+library)
        panel_specs.append(('Measured pathway enrichment',lambda ax:dotplot(ax,pathways[pathways.scope=='global'],'Global pathway terms')))
    # Assemble only measured panels. No dummy data or blank scientific placeholders.
    ncols=3;nrows=math.ceil(len(panel_specs)/ncols);fig,axes=plt.subplots(nrows,ncols,figsize=(21,5.3*nrows));axes=np.asarray(axes).ravel()
    for i,(title,draw) in enumerate(panel_specs):
        draw(axes[i]);axes[i].set_title(f'{chr(65+i)}. {title}',fontsize=10,loc='left')
    for ax in axes[len(panel_specs):]:ax.remove()
    fig.suptitle(f'TabulaStromal {ctx.lineage}: measured analysis panels\nExternal combined no adaptation; all model/K/fold selection uses internal CV',fontsize=14);fig.tight_layout(rect=[0,.015,1,.96]);outputs+=save(ctx,fig,'Fig3_TabulaStromal_'+ctx.lineage)
    write_json(ctx.root/'figures/figure_manifest.json',dict(paths=[str(p) for p in outputs],included_panels=[p[0] for p in panel_specs],
                    missing_components=[name for name,df in [('stable_pairs',stable),('pathway_terms',pathways),('aging_overlap',overlap)] if df.empty],
                    note='Missing result panels are omitted and listed explicitly; no result is simulated. MAD is age-label dispersion.'))
    write_text(ctx.root/'figures/figure_legends.md','''# Figure definitions

Benchmark: each model's K chosen from internal CV; fivefold summary = mean±descriptive SD of correlated folds, not an independent-replicate CI. Representative fold = highest internal validation R2, then lowest MAE, then lowest fold. External-valid always means combined external cohort without adaptation. The externally selected-fold CSV is exploratory and is never used here or for final model choice.
MAD/Age_MAD = mean absolute deviation of true age from cohort mean; MedAE = median absolute prediction error. MAD is not a within-cohort model ranking score.
Scatter regressions and y=x are descriptive. All predictions, including negative ages, are retained. Adaptation examples compare the same held-out samples; fraction=.3, repeat=0 is fixed before evaluating outcomes.
Stability intervals are descriptive after feature selection; subtype coverage is candidate-normalized with raw counts saved separately. Prediction-associated genes are not causal aging drivers. Gene-set panels preserve provided_unverified labels; all pathway tests use training expression opportunity backgrounds.
''')
    return outputs+[ctx.root/'figures/figure_manifest.json']
