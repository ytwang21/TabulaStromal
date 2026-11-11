"""Connect FT30 CSVs to the original figure generator and its K/benchmark panels."""
from pathlib import Path
import json,shutil
import numpy as np
import pandas as pd
from tabulastromal import figures as original
from tabulastromal.storage import write_csv,write_json,write_text
import ft30

LABEL='External (FT30 → all 100%; includes adaptation)'
BASELINE='External all100 (no adaptation; best internal fold)'
NOTE='Primary external_all_post_adaptation uses all100 samples, INCLUDING the fixed30% fitted by continuation/refit. Failures use the SAME original best internal-Valid fold with explicit fallback status. Train/Valid: unchanged fivefold mean ± descriptive SD; External: one internally selected source fold, no fold SD. All100 scores mix apparent/unseen predictions, not independent external validation. Separate external_heldout_70 tables/figures report never-adapted samples.'

def prediction_scatter(ax,table,title):
    """High-contrast groups, with explicit all100 labels for paired displays."""
    if table.empty:return
    markers={'Train':'o','Valid':'^','External-valid':'s','Baseline all100':'o','Updated all100':'s'}
    for split,sub in table.groupby('split',sort=False):
        m=original.evaluate(sub.Age,sub.prediction)
        color=original.COLORS.get(split,'#009E73')
        ax.scatter(sub.Age,sub.prediction,s=12 if split=='Train' else 23,
            alpha=.5 if split=='Train' else .85,color=color,marker=markers.get(split,'s'),
            label=f'{split}: n={len(sub)}, R²={m["R2"]:.2f}, MAE={m["MAE"]:.2f}')
    v=table[['Age','prediction']].to_numpy(float);lo,hi=np.min(v),np.max(v)
    ax.plot([lo,hi],[lo,hi],'--',color='.5',lw=.8)
    ax.set_xlabel('Chronological age (years)');ax.set_ylabel('Predicted age (years)')
    ax.set_title(title);ax.legend(fontsize=6,frameon=False)

def backup_benchmarks(ctx):
    for stem in ['benchmark_fivefold_summary','feature_K_comparison']:
        for suffix in ['png','pdf']:
            p=ctx.root/'figures'/f'{stem}.{suffix}';backup=p.with_name(stem+'.pre_ft30.'+suffix)
            if p.exists() and not backup.exists():shutil.copy2(p,backup)

def frames(ctx,smoke=False):
    cv=pd.read_csv(ctx.root/'benchmark/cv_fold_metrics.csv')
    extname='external_ft30_all_sample_metrics.csv' if ft30.OPTIONS['external_eval_mode']=='adapt30_all100' else 'external_combined_no_adapt_all100.csv'
    if smoke:extname='external_combined_ft30_smoke_per_fold.csv' if ft30.OPTIONS['external_eval_mode']=='adapt30_all100' else 'external_combined_ft30_smoke_baseline.csv'
    ext=pd.read_csv(ctx.root/'external'/extname)
    return cv,ext[ext.status=='complete'].copy()

def benchmark_plots(ctx,partial=False):
    """Also used by smoke tests: the actual figure data source is saved alongside it."""
    original.style();original.COLORS.update({'Train':'#0072B2','Valid':'#D55E00','External-valid':'#009E73'});import matplotlib.pyplot as plt
    cv,ext=frames(ctx,smoke=partial)
    ranking=pd.read_csv(ctx.root/'benchmark/model_comparison.csv')
    selected=ranking.sort_values('internal_rank').groupby('model',sort=False).head(1)[['model','K']]
    label=LABEL if ft30.OPTIONS['external_eval_mode']=='adapt30_all100' else BASELINE
    # A smoke figure is clearly named and cannot overwrite a complete benchmark.
    if partial:
        pairs=ext[['model','K']].drop_duplicates();table=pd.concat([cv.merge(pairs,on=['model','K']),ext.assign(split='External-valid')],ignore_index=True)
        order=ext.model.drop_duplicates().tolist();name='ft30_smoke_benchmark'
    else:
        table=pd.concat([cv.merge(selected,on=['model','K']),ext.merge(selected,on=['model','K']).assign(split='External-valid')],ignore_index=True)
        order=selected.model.tolist();name='benchmark_fivefold_summary_ft30_all'
    write_csv(ctx.root/'figures'/(name+'_plot_data.csv'),table)
    fig,axes=plt.subplots(2,2,figsize=(17,10))
    for ax,metric in zip(axes.ravel(),['R2','MAE','MedAE','PCC']):original.plot_benchmark(ax,table,metric,order)
    for ax in axes.ravel():
        for c in ax.containers:
            if hasattr(c,'get_label') and c.get_label()=='External-valid':c.set_label(label)
    axes[0,0].legend(frameon=False,ncol=3,fontsize=7)
    fig.suptitle(f'{ctx.lineage}: '+('SMOKE TEST — ElasticNet one source fold' if partial else 'fivefold Train/Valid and best-fold External all100'),fontsize=12)
    fig.text(.02,.005,NOTE+'\nMAD = mean absolute age-label deviation; MedAE = median absolute prediction error.',fontsize=7)
    fig.tight_layout(rect=[0,.045,1,.95]);outputs=original.save(ctx,fig,name)
    if not partial:
        originals=list(outputs)
        for stem in ['benchmark_fivefold_summary','benchmark_fivefold_summary_ft30']:
            for p in originals:
                target=p.with_name(stem+p.suffix);shutil.copy2(p,target);outputs.append(target)
            dp=ctx.root/'figures'/(stem+'_plot_data.csv');write_csv(dp,table);outputs.append(dp)
        outputs.append(ctx.root/'figures'/(name+'_plot_data.csv'))
        outputs+=k_plot(ctx,cv,ext,label)
        delta=pd.read_csv(ctx.root/'external/ElasticNet_ft30_paired_comparison.csv')
        held=pd.read_csv(ctx.root/'external/external_heldout_70_paired_metrics.csv')
        held=held[held.model=='ElasticNet']
        fig,axes=plt.subplots(1,2,figsize=(12,5))
        for ax,metric in zip(axes,['R2','MAE']):
            for frame,col,title,color,ls in [
                (delta,metric,'EN refit → all100','#009E73','-'),
                (delta,'baseline_'+metric,'EN original → all100','#0072B2','--'),
                (held,metric,'EN refit → heldout70','#D55E00','-'),
                (held,'baseline_'+metric,'EN original → heldout70','#666666','--')]:
                t=frame.sort_values('K');ax.plot(t.K,t[col],marker='o',label=title,color=color,ls=ls)
            ax.set_xlabel('K (source-selected Celltype_Gene features)');ax.set_ylabel(metric);ax.legend(fontsize=7)
        fig.suptitle(f'{ctx.lineage}: ElasticNet reference refit; all100 and separate heldout70')
        fig.tight_layout();outputs+=original.save(ctx,fig,'ElasticNet_ft30_paired_comparison')
        outputs+=heldout_plot(ctx)

    return outputs

def k_plot(ctx,cv,ext,label):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(15,13))
    tables=[('Train',cv[cv.split=='Train']),('Valid',cv[cv.split=='Valid']),(label,ext)]
    colors=plt.get_cmap('tab20')
    for row,(title,table) in enumerate(tables):
        for col,metric in enumerate(['R2','MAE']):
            ax=axes[row,col]
            for i,name in enumerate(ft30.models.NAMES):
                sub=table[table.model==name]
                if sub.empty:continue
                s=sub.groupby('K')[metric].agg(['mean','std']);ax.errorbar(s.index,s['mean'],yerr=s['std'].fillna(0),marker='.',lw=.8,color=colors(i),label=name)
            ax.set_title(title+' — '+metric);ax.set_xlabel('Requested Celltype_Gene K');ax.set_ylabel(metric)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=5,fontsize=7)
    fig.suptitle(f'{ctx.lineage}: training, validation and fixed FT30 K comparison');fig.tight_layout(rect=[0,.08,1,.97])
    outputs=original.save(ctx,fig,'feature_K_comparison_ft30_all')
    originals=list(outputs)
    for stem in ['feature_K_comparison','feature_K_comparison_ft30']:
        for p in originals:
            target=p.with_name(stem+p.suffix);shutil.copy2(p,target);outputs.append(target)
    for stem in ['feature_K_comparison','feature_K_comparison_ft30','feature_K_comparison_ft30_all']:
        dp=ctx.root/'figures'/(stem+'_plot_data.csv');write_csv(dp,pd.concat([cv,ext.assign(split=label)],ignore_index=True));outputs.append(dp)
    return outputs

def heldout_plot(ctx):
    import matplotlib.pyplot as plt
    table=pd.read_csv(ctx.root/'external/external_heldout_70_paired_metrics.csv')
    rank=pd.read_csv(ctx.root/'benchmark/model_comparison.csv').sort_values('internal_rank')
    pairs=rank.groupby('model',sort=False).head(1)[['model','K']]
    table=pairs.merge(table,on=['model','K'])
    fig,axes=plt.subplots(2,2,figsize=(17,10))
    for ax,metric in zip(axes.ravel(),['R2','MAE','MedAE','PCC']):
        ix=np.arange(len(table));ax.bar(ix-.18,table['baseline_'+metric],.36,color='#0072B2',label='Original best fold: same heldout70')
        ax.bar(ix+.18,table[metric],.36,color='#D55E00',label='After adaptation/refit: heldout70')
        ax.set_xticks(ix);ax.set_xticklabels(table.model,rotation=55,ha='right',fontsize=7);ax.set_ylabel(metric)
        ax.axhline(0,color='.7',lw=.5)
    axes[0,0].legend(fontsize=7,frameon=False)
    fig.suptitle(f'{ctx.lineage}: external_heldout_70 — NEVER used for adaptation; separate from all100 primary figure')
    fig.tight_layout(rect=[0,.03,1,.95]);outputs=original.save(ctx,fig,'external_heldout_70_comparison')
    p=ctx.root/'external/external_heldout_70_comparison.csv';write_csv(p,table);outputs.append(p)
    p=ctx.root/'figures/external_heldout_70_comparison_plot_data.csv';write_csv(p,table);outputs.append(p)
    return outputs


def render(ctx):
    ft30.validate_predictions(ctx,json.loads((ctx.root/'external/ft30_protocol.json').read_text())['ft30_protocol'])
    cv,ext=frames(ctx)
    expected=len(ft30.selected_folds(ctx))
    if len(ext)!=expected:raise RuntimeError('Cannot render complete benchmark with missing FT30 tasks')
    # Preserve original output references without creating another result directory.
    for stem in ['benchmark_fivefold_summary','feature_K_comparison','Fig3_TabulaStromal_'+ctx.lineage]:
        for suffix in ['png','pdf']:
            p=ctx.root/'figures'/f'{stem}.{suffix}';backup=p.with_name(stem+'.pre_ft30.'+suffix)
            if p.exists() and not backup.exists():shutil.copy2(p,backup)
    import ft30_figure_panels as panels
    panels.scatter=prediction_scatter
    original.COLORS.update({'Train':'#0072B2','Valid':'#D55E00','External-valid':'#009E73','Baseline all100':'#8D6E63','Updated all100':'#009E73'})
    oldread,oldsave,oldtext,oldjson=panels.read,panels.save,panels.write_text,panels.write_json
    pred=pd.read_csv(ctx.root/'external/external_ft30_all_sample_predictions.csv')
    label=LABEL if ft30.OPTIONS['external_eval_mode']=='adapt30_all100' else BASELINE
    if ft30.OPTIONS['external_eval_mode']=='no_adapt_all100':pred=pred.assign(prediction=pred.baseline_prediction)
    def read(path):
        name=Path(path).name
        if name=='external_combined_no_adapt.csv':return ext.copy()
        if name in ['external_combined_no_adapt_predictions.csv','external_combined_adaptation_predictions.csv']:return pred.copy()
        if name in ['external_per_gse_adaptation.csv','external_per_gse_no_adapt.csv']:return pd.DataFrame()
        return oldread(path)
    def save(context,fig,name):
        import matplotlib.pyplot as plt
        if name=='feature_K_comparison':plt.close(fig);return []
        for text in fig.findobj(match=__import__('matplotlib').text.Text):
            value=text.get_text().replace('External-valid = external combined, no adaptation',label)
            value=value.replace('External combined no adaptation;',label+';')
            value=value.replace('External-valid',label)
            value=value.replace('External combined cohort without adaptation',label)
            text.set_text(value)
        return oldsave(context,fig,name)
    def text(path,value):
        if Path(path).name=='figure_legends.md':
            value=value.replace('External-valid always means combined external cohort without adaptation.',NOTE)
            value=value.replace('fraction=.3, repeat=0','fixed fraction=.3, one seed-42 split')
        return oldtext(path,value)
    def metadata(path,value):
        if Path(path).name=='figure_manifest.json':value.update(external_protocol='FT30',external_eval_mode=ft30.OPTIONS['external_eval_mode'],external_label=label)
        return oldjson(path,value)
    panels.read,panels.save,panels.write_text,panels.write_json=read,save,text,metadata
    try:outputs=panels.figures_from_saved(ctx)
    finally:panels.read,panels.save,panels.write_text,panels.write_json=oldread,oldsave,oldtext,oldjson
    outputs+=benchmark_plots(ctx)
    write_text(ctx.root/'figures/figure_legends.md',NOTE+'\nMAD/Age_MAD = mean absolute age-label dispersion. MedAE = median absolute prediction error. Benchmark CSVs unchanged.\nScatter uses Train/Valid from one best internal fold and all100 external from the same source fold.\nBiological panels remain the original training-only EN stability/enrichment analysis. Per-GSE outputs are unchanged.\n')
    manifest=ctx.root/'figures/figure_manifest.json';d=json.loads(manifest.read_text());d['paths']=list(dict.fromkeys(d['paths']+[str(p) for p in outputs]));d['external_note']=NOTE
    write_json(manifest,d)
    return list(dict.fromkeys(outputs+[manifest,ctx.root/'figures/figure_legends.md']))
