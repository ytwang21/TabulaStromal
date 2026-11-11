#!/usr/bin/env python3
"""One-page A–D signature figure; global training-selected lists and real scores.

A/B: Fibroblast/Endothelial functional and subtype UMAPs.
C/D: all Fibroblast/Endothelial Leiden clusters, respectively.
The marker search uses training GSE only. Cell scores use expression alone.
"""
from __future__ import annotations
import logging
import math
from pathlib import Path
import time
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import Normalize, to_rgb, to_hex
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import signature_pipeline as sp
import marker_combination_search as search

BAR_ORDER = {'fibro':['Myfibro','iFibro','eFibro','qFibro'],
             'endo':['MCs','LEC','VenEC','CapEC','ArtEC','EndoMT']}
UMAP_LAYOUT = {'fibro':['ECM','AP','Myfibro','eFibro','iFibro','qFibro'],
               'endo':['AP','Tip','Stalk','ArtEC','CapEC','VenEC','EndoMT','MCs','LEC']}


def pale(color, alpha=.13):
    return to_hex(np.asarray(to_rgb(color))*alpha+(1-alpha))


def load_inputs(cfg, lin):
    c=sp.cache(cfg,lin);audit=sp.read_json(sp.root(cfg)/f'audit/data_audit_{lin}.json')
    fp=sp.fingerprint(cfg['inputs'][lin])
    if any(audit[k]!=v for k,v in fp.items()):
        raise ValueError('Source fingerprint changed; audit/rescore before plotting')
    obs=pd.read_parquet(c/'metadata.parquet')
    scores=np.load(c/'figure_mean.npy',mmap_mode='r')
    names=np.load(c/'figure_score_names.npy').tolist()
    if len(obs)!=len(scores):raise ValueError('Score/metadata row alignment failed')
    meta,means=sp.aggregate_scores(scores,obs,['leiden'])
    compositions=obs.groupby(['leiden','Celltype','Celltype_general'],observed=True).size().rename('n_cells').reset_index()
    bar_order=BAR_ORDER[lin];means=means[:,[names.index(s) for s in bar_order]]
    records=[];long=[]
    for i,r in meta.iterrows():
        cluster=str(r['leiden']);part=compositions[compositions['leiden'].astype(str).eq(cluster)]
        part=part.sort_values('n_cells',ascending=False);dominant=part.iloc[0]
        fine=str(dominant['Celltype']);ref=str(dominant['Celltype_general'])
        nref=part['Celltype_general'].nunique();nfine=part['Celltype'].nunique()
        v=means[i];best=float(np.nanmax(v))
        top=[s for s,value in zip(bar_order,v) if np.isfinite(value) and np.isclose(value,best,atol=1e-10,rtol=0)]
        unique=len(top)==1;match=unique and top[0]==ref and nref==1
        title=f'{fine}_C{cluster}';own=float(v[bar_order.index(ref)])
        other=float(np.nanmax(np.delete(v,bar_order.index(ref))))
        record={'lineage':lin,'leiden':cluster,'Celltype':fine,'Celltype_general':ref,
                'panel_title':title,'n_cells':int(r['n_cells']),'n_Celltype':int(nfine),
                'n_Celltype_general':int(nref),'dominant_reference_fraction':float(dominant['n_cells']/r['n_cells']),
                'top_signature':';'.join(top),'is_unique_highest':unique,
                'annotation_matches_highest':match,'annotation_score':own,
                'highest_other_score':other,'annotation_minus_highest_other':own-other,
                'plot_method':'all-cell mean of training-selected global mean-log signatures',
                'score_source':str(c/'figure_mean.npy'),
                'reference_composition':part[['Celltype','Celltype_general','n_cells']].to_json(orient='records')}
        record.update({f'mean_{s}':float(v[j]) for j,s in enumerate(bar_order)})
        records.append(record)
        for j,s in enumerate(bar_order):
            long.append({'lineage':lin,'leiden':cluster,'panel_title':title,'Celltype':fine,
                         'Celltype_general':ref,'signature':s,'mean_score':float(v[j]),
                         'n_cells':int(r['n_cells']),'annotation_matches_highest':match,
                         'scoring_method':'mean_log_expression','aggregation':'all cells in Leiden; equal weight per cell',
                         'annotation_used_for_list_selection':'training GSE only',
                         'annotation_used_to_change_score':False})
    table=pd.DataFrame(records);table['subtype_order']=table['Celltype_general'].map({s:j for j,s in enumerate(bar_order)})
    table['numeric_cluster']=pd.to_numeric(table['leiden'],errors='coerce')
    table=table.sort_values(['subtype_order','numeric_cluster','leiden']).reset_index(drop=True)
    table=table.drop(columns=['subtype_order','numeric_cluster'])
    table.to_csv(sp.root(cfg)/f'data/plot_reference_cluster_check_{lin}.csv',index=False)
    pd.DataFrame(long).to_csv(sp.root(cfg)/f'data/plot_reference_cluster_scores_{lin}.csv',index=False)
    logging.info('%s: %d/%d Leiden match; all clusters retained',lin,int(table.annotation_matches_highest.sum()),len(table))
    return obs,table


def select_umap_signatures(cfg,lin,rows):
    c=sp.cache(cfg,lin);specs=sp.read_json(c/'score_specs_indexed.json')
    candidates=np.load(c/'candidate_mean.npy',mmap_mode='r')
    values=np.load(c/'figure_mean.npy',mmap_mode='r');names=np.load(c/'figure_score_names.npy').tolist()
    doc=sp.read_json(sp.root(cfg)/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_figure.json')
    selected={};records=[]
    for key in UMAP_LAYOUT[lin]:
        if key in names:
            selected[key]=np.asarray(values[rows,names.index(key)])
            genes=doc['marker_sets'][key]
            rec={'candidate_id':doc['version'],'n_genes_listed':len(genes),'n_genes_matched':len(genes),
                 'genes_used':';'.join(genes),'source':str(c/'figure_mean.npy'),'role':'subtype; same lists/scores as bars'}
        else:
            dic='Fibro_marker_sub_ya' if lin=='fibro' else 'score_dict'
            ss=[s for s in specs if s['mapping']=='original' and s['dictionary_name']==dic and s['key']==key]
            if not ss:raise ValueError(f'Missing recovered signature {lin}:{dic}:{key}')
            ss.sort(key=lambda s:(not s['candidate_id'].startswith(lin+':'),s['candidate_id']))
            s=ss[0];selected[key]=np.asarray(candidates[rows,s['score_index']])
            rec={'candidate_id':s['candidate_id'],'n_genes_listed':len(s['genes_listed']),
                 'n_genes_matched':len(s['genes']),'genes_used':';'.join(s['genes']),
                 'source':str(c/'candidate_mean.npy'),'role':'historical functional signature; not tuned'}
        records.append({'lineage':lin,'signature':key,'method':'mean_log_expression',**rec})
    pd.DataFrame(records).to_csv(sp.root(cfg)/f'data/plot_reference_umap_signatures_{lin}.csv',index=False)
    return selected


def umap_panel(fig,cfg,lin,obs,rect):
    """Rect uses (left,right,bottom,top) figure fractions."""
    rows=np.load(sp.cache(cfg,lin)/'plot_rows.npy')
    xy=np.load(sp.cache(cfg,lin)/'umap.npy',mmap_mode='r')[rows]
    vals=select_umap_signatures(cfg,lin,rows)
    finite=np.concatenate([v[np.isfinite(v)] for v in vals.values()])
    vmax=max(1.,float(math.ceil(np.quantile(finite,.995)))) if len(finite) else 1.
    norm=Normalize(0,vmax);nrows=2 if lin=='fibro' else 3
    grid=GridSpec(nrows,3,figure=fig,left=rect[0],right=rect[1],bottom=rect[2],top=rect[3],wspace=.18,hspace=.20)
    # Equal plotted limits for every signature of the same lineage; no per-map crop.
    xmin,ymin=np.nanmin(xy,axis=0);xmax,ymax=np.nanmax(xy,axis=0)
    dx=xmax-xmin;dy=ymax-ymin
    for i,key in enumerate(UMAP_LAYOUT[lin]):
        nested=GridSpecFromSubplotSpec(1,2,subplot_spec=grid[i//3,i%3],width_ratios=[1,.035],wspace=.025)
        ax=fig.add_subplot(nested[0]);cax=fig.add_subplot(nested[1])
        scatter=ax.scatter(*xy.T,c=vals[key],s=.65,cmap='viridis',norm=norm,linewidths=0,rasterized=True)
        cb=fig.colorbar(scatter,cax=cax,extend='max');cb.ax.tick_params(labelsize=9,length=2,pad=2)
        cb.locator=MaxNLocator(nbins=4);cb.update_ticks()
        if cb.solids is not None:cb.solids.set_rasterized(False)
        ax.set_title(key,fontsize=15,pad=4,fontweight='normal')
        ax.set_xlim(xmin-dx*.03,xmax+dx*.03);ax.set_ylim(ymin-dy*.03,ymax+dy*.03)
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
    sp.dump(sp.root(cfg)/f'data/plot_reference_umap_display_{lin}.json',
            {'color_vmin':0,'color_vmax':vmax,'color_strategy':'shared within lineage, ceil pooled99.5percentile; display saturation only',
             'n_sampled_cells':len(rows),'n_atlas_cells':len(obs),'panel_layout':UMAP_LAYOUT[lin],
             'subtype_score_source':str(sp.cache(cfg,lin)/'figure_mean.npy'),
             'method':'mean log expression; no clipping or annotation-dependent score changes'})


def format_value(value):
    if not np.isfinite(value):return 'NA'
    if value==0:return '0'
    if abs(value)<.001:return f'{value:.1e}'
    if abs(value)<.1:return f'{value:.3f}'
    return f'{value:.2f}'


def bar_panel(fig,cfg,lin,table,rect,ncols=6,nrows=None,legend_y=.025):
    order=BAR_ORDER[lin];pal=sp.read_json(sp.root(cfg)/'audit/fig1_v2_colors.json')['general_colors']
    nrows=nrows or math.ceil(len(table)/ncols)
    grid=GridSpec(nrows,ncols,figure=fig,left=rect[0],right=rect[1],bottom=rect[2],top=rect[3],wspace=.30,hspace=.64)
    values=table[[f'mean_{s}' for s in order]].to_numpy()
    if np.nanmin(values)<0:raise ValueError('Negative mean-log score; refuse to hide negative scores')
    ylim=float(np.nanmax(values))*1.25
    for i,row in table.iterrows():
        ax=fig.add_subplot(grid[i//ncols,i%ncols]);v=values[i]
        ax.set_facecolor(pale(pal[row['Celltype_general']]))
        bars=ax.bar(np.arange(len(order)),v,color=[pal[s] for s in order],width=.69,edgecolor='#333333',linewidth=.45)
        best=float(np.nanmax(v))
        for j,(bar,val) in enumerate(zip(bars,v)):
            top=np.isfinite(val) and np.isclose(val,best,atol=1e-10,rtol=0)
            if top:bar.set_linewidth(1.1)
            if np.isfinite(val):
                text=format_value(val)
                ax.text(j,val+ylim*.022,text,ha='center',va='bottom',fontsize=10 if top else (7 if len(text)>5 else 8),color='black')
            else:ax.text(j,ylim*.02,'NA',ha='center',fontsize=8)
        suffix='' if row['annotation_matches_highest'] else ' †'
        if row['n_Celltype_general']>1:suffix+=' [mixed]'
        title=row['panel_title']+suffix
        ax.set_title(title,fontsize=9 if len(title)<36 else 8,pad=5,fontweight='normal')
        ax.set_ylim(0,ylim);ax.set_xlim(-.65,len(order)-.35)
        ax.set_xticks([]);ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis='y',labelsize=8,length=2,pad=2)
        if i%ncols==0:ax.set_ylabel('Signature score',fontsize=9,labelpad=3)
        ax.grid(axis='y',lw=.45,alpha=.27);ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(.6);ax.spines['bottom'].set_linewidth(.6)
    mid=(rect[0]+rect[1])/2;width=rect[1]-rect[0]
    for handles,x,title in [([Patch(facecolor=pal[s],label=s) for s in order],mid-width*.25,'marker_list'),
                            ([Patch(facecolor=pale(pal[s]),label=s) for s in order],mid+width*.25,'Celltype')]:
        legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(x,legend_y),ncol=2 if len(order)==4 else 3,
                          fontsize=12,title=title,title_fontsize=13,frameon=True,fancybox=True,
                          handlelength=1.1,columnspacing=1.0,framealpha=1,borderpad=.6)
        legend.get_frame().set_edgecolor('#777777');legend.get_frame().set_linewidth(.7)


def make_figures(cfg,loaded):
    plt=sp.plotting_setup();r=sp.root(cfg)
    # One page, four panels. All 99 clusters remain, with space for full titles.
    fig=plt.figure(figsize=(36,27))
    for lin,left,right,a,b in [('fibro',.040,.477,'A','C'),('endo',.547,.984,'B','D')]:
        obs,table=loaded[lin]
        fig.text(left-.025,.978,a,fontsize=31,va='center',fontweight='normal')
        fig.text((left+right)/2,.978,sp.LABELS[lin],fontsize=21,ha='center',va='center')
        umap_panel(fig,cfg,lin,obs,(left,right,.665,.948))
        fig.text(left-.025,.617,b,fontsize=31,va='center',fontweight='normal')
        fig.text((left+right)/2,.617,'Signature score of '+sp.LABELS[lin]+' subtype markers',fontsize=19,ha='center',va='center')
        bar_panel(fig,cfg,lin,table,(left,right,.105,.585),ncols=6 if lin=='fibro' else 5,nrows=10,legend_y=.020)
        fig.text((left+right)/2,.075,
            f"All Leiden clusters: {int(table.annotation_matches_highest.sum())}/{len(table)} match reference; † discordant",
            fontsize=12,ha='center',va='center')
    dest=r/'figures/integrated/ExtFig1_SignatureScore_Fibro_Endo'
    sp.save_figure(fig,dest,cfg);plt.close(fig)
    # Separate downloads inherit the same A/C and B/D panel letters.
    for lin,a,b in [('fibro','A','C'),('endo','B','D')]:
        obs,table=loaded[lin];fig=plt.figure(figsize=(20,24));left=.060;right=.975
        fig.text(.018,.979,a,fontsize=29,va='center')
        fig.text(.52,.979,sp.LABELS[lin],fontsize=21,ha='center',va='center')
        umap_panel(fig,cfg,lin,obs,(left,right,.668,.945))
        fig.text(.018,.620,b,fontsize=29,va='center')
        fig.text(.52,.620,'Signature score of '+sp.LABELS[lin]+' subtype markers',fontsize=19,ha='center',va='center')
        bar_panel(fig,cfg,lin,table,(left,right,.105,.585),ncols=6 if lin=='fibro' else 5,nrows=10,legend_y=.020)
        fig.text(.52,.075,f"All Leiden clusters: {int(table.annotation_matches_highest.sum())}/{len(table)} match reference; † discordant",fontsize=12,ha='center',va='center')
        sp.save_figure(fig,r/f'figures/panels/ExtFig1_SignatureScore_{sp.LABELS[lin]}',cfg);plt.close(fig)


def write_figure_report(cfg,loaded):
    r=sp.root(cfg);table=pd.concat([loaded[lin][1] for lin in ['fibro','endo']],ignore_index=True)
    table.to_csv(r/'data/plot_reference_cluster_check.csv',index=False)
    lines=['# A–D layout and global marker comparison — 2026-11-06','',
        '原路径PNG/PDF已覆盖；主图为一页A–D，按上排A/B、下排C/D排列。分图沿用主图字母A/C与B/D。','',
        'A：Fibroblast ECM、AP、Myfibro、eFibro、iFibro、qFibro，2×3 UMAP。',
        'B：Endothelial AP、Tip、Stalk、ArtEC、CapEC、VenEC、EndoMT、MCs、LEC，3×3 UMAP。',
        'C/D：分别为Fibroblast/Endothelial全部Leiden；每个标题严格使用`Celltype + "_C" + leiden`。', '',
        '所有亚型UMAP与柱图使用同一套新marker、同一批真实单细胞评分。柱值为该Leiden全部细胞的未加权mean-log表达均值；柱色继承Fig1 v2；浅色背景为原Celltype_general。黑色轮廓与较大数值标出真实最高值；†标记与参考注释不一致。', '',
        '## Marker组合比较方法与边界','',
        '从真实notebook、用户候选表、Fig1 v2和既有训练gene-QC构建有限候选池；每个亚型候选先接受训练GSE的表达覆盖、配对效应与稳定性QC；之前的基线明确保留为对照。以训练GSE等权一致率及亚型平衡项选择全局组合，每个亚型一份名单用于所有细胞。Fibroblast穷举固定候选池；Endothelial使用有停止条件的多起点坐标搜索，不声称穷举全部组合。', '',
        '搜索没有使用Leiden编号，不以全Leiden达到100%为停止条件。训练时使用了Celltype_general参考注释，因此是监督式名单选择；原标签又是marker-informed，存在循环性。冻结JSON之后才评价保留GSE，再绘制全部Leiden。原始final JSON与原评分保留作基线，新图使用`*_subtype_markers_figure.json`与`figure_mean.npy`。', '',
        '不同基因集的原始均值没有跨signature校准；最高柱仅是一致性诊断，不是独立身份金标准。没有可信donor字段，保留验证以GSE代替。既有保留GSE在此前基线分析中已被查看，本次属于回顾性分组评价，不是新的盲测。不能为了全部最高柱吻合而逐cluster改名单、给参考标签加分或隐藏反例。所有不一致均保留。', '']
    for lin in ['fibro','endo']:
        doc=sp.read_json(r/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_figure.json');t=loaded[lin][1];bad=t[~t.annotation_matches_highest]
        baseline=np.load(sp.cache(cfg,lin)/'final_mean.npy',mmap_mode='r');names=np.load(sp.cache(cfg,lin)/'final_score_names.npy').tolist()
        obs=loaded[lin][0];bm,bv=sp.aggregate_scores(baseline,obs,['leiden'])
        refs=obs.groupby('leiden',observed=True).Celltype_general.first()
        oldmatch=sum(names[int(np.argmax(bv[i]))]==str(refs.loc[row.leiden]) for i,row in bm.iterrows())
        lines += [f'## {sp.LABELS[lin]}',f"实际比较{doc['selection']['trials_evaluated']}组；全部Leiden一致数：旧名单{oldmatch}/{len(t)} → 新名单{int(t.annotation_matches_highest.sum())}/{len(t)}。",'',
                  sp.mdtable(pd.DataFrame(doc['comparison'])),'',
                  '### 本次绘图marker与保留集状态','',sp.mdtable(pd.DataFrame([
                      {'subtype':s,'genes':', '.join(gg),**doc['holdout_validation'][s]} for s,gg in doc['marker_sets'].items()])),'',
                  '### 保留的不一致Leiden','',sp.mdtable(bad[['leiden','Celltype','Celltype_general','top_signature','annotation_score','highest_other_score']]),'']
    lines += ['## 可复查文件','',
        '- `marker_sets/*_subtype_markers_figure.json`：本次使用的全局名单、来源、冻结及保留验证信息。',
        '- `data/marker_combination_candidate_QC_*.csv`：所有候选组合组件的训练QC。',
        '- `data/marker_combination_search_*.csv`：每次全局组合试验、训练目标与没有使用保留集/Leiden的记录。',
        '- `data/marker_combination_train_holdout_*.csv`：冻结后新旧名单比较及GSE bootstrap95%CI。',
        '- `data/plot_reference_cluster_scores_*.csv`：每个Leiden每根柱的精确数值。',
        '- `data/plot_reference_cluster_check.csv`：全部cluster、原标签、最高signature、差值及完整标题。',
        '- `data/plot_reference_umap_signatures_*.csv`：A/B每个UMAP真实来源与基因列表。',
        '- `scripts/marker_combination_search.py`、`scripts/reference_style_figures.py`：可复现选择与绘图代码。', '',
        'PNG 300dpi；主PDF一页，文字、数值、轴、柱形、背景、图例与色条可编辑，仅UMAP百万级点云局部栅格化。']
    (r/'reports/reference_layout_report.md').write_text('\n'.join(lines)+'\n')
    (r/'reports/figure_legend.md').write_text('\n'.join(lines[:13])+'\n\nMarkers and complete validation: reports/reference_layout_report.md. Training-supervised selection; discordant clusters retained.\n')
    summary={'layout':'one_page_A_B_C_D','panel_A':UMAP_LAYOUT['fibro'],'panel_B':UMAP_LAYOUT['endo'],
             'panel_C':'all Fibroblast Leiden bars','panel_D':'all Endothelial Leiden bars',
             'scoring_method':'unweighted mean log expression; global training-selected lists',
             'annotation_changed':False,'label_dependent_score_boost':False,'pdf_pages':1,
             'panels':{lin:{'n_clusters':len(loaded[lin][1]),'n_matching_clusters':int(loaded[lin][1].annotation_matches_highest.sum())} for lin in ['fibro','endo']},
             'updated_at':time.strftime('%Y-%m-%d %H:%M:%S')}
    sp.dump(r/'logs/reference_layout_update.json',summary)
    refresh_current_documentation(cfg)


def refresh_current_documentation(cfg):
    """Keep historical baseline statistics clearly separate from current figures."""
    r=sp.root(cfg)
    if not (r/'reports/reference_layout_report.md').exists():return
    report=(r/'reports/reference_layout_report.md').read_text()
    (r/'reports/figure_legend.md').write_text('\n'.join(report.splitlines()[:13])+\
        '\n\nCurrent global lists, held-out results and all discordant clusters: reference_layout_report.md. '
        'A/B use the same subtype scores as C/D. PDF: one page, vector bars and text; point clouds locally rasterized.\n')
    summary=sp.read_json(r/'logs/reference_layout_update.json')
    figures={lin:sp.read_json(r/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_figure.json') for lin in ['fibro','endo']}
    rates={lin:next(x['GSE_equal_top1'] for x in figures[lin]['comparison'] if x['split']=='holdout' and x['marker_version']=='selected_global_combination') for lin in figures}
    matched={lin:f"{summary['panels'][lin]['n_matching_clusters']}/{summary['panels'][lin]['n_clusters']}" for lin in figures}
    note=('<!-- CURRENT_FIGURE_START -->\n'
          '## 当前绘图版本：2026-11-06，单页 A–D\n\n'
          'A为Fibroblast ECM/AP/4亚型UMAP；B为Endothelial AP/Tip/Stalk/6亚型UMAP；C/D分别为全部Leiden柱状图。'
          '本图使用新的训练GSE组合选择名单 `*_subtype_markers_figure.json` 与 `figure_mean.npy`，'
          '与下方历史基线 `*_subtype_markers_final.json` / mean/module统计区分。'
          f"新名单保留GSE一致率：Fibroblast {rates['fibro']:.2%}，Endothelial {rates['endo']:.2%}；全部Leiden一致数{matched['fibro']}、{matched['endo']}。"
          '†反例保留，没有按标签加分。详见 [当前图形与组合验证报告](reference_layout_report.md)。\n\n'
          '以下正文保留2026-11-04基线Phase0–5分析；其分数及状态不应解释为本次新名单的统计。\n'
          '<!-- CURRENT_FIGURE_END -->\n\n')
    for name in ['run_summary.md','marker_selection_report.md','annotation_validation_report.md']:
        path=r/'reports'/name;text=path.read_text()
        if '<!-- CURRENT_FIGURE_END -->' in text:text=text.split('<!-- CURRENT_FIGURE_END -->',1)[1].lstrip()
        # Retire the superseded, manually prepended C/D update block.
        if name=='run_summary.md':
            text='\n'.join(line for line in text.splitlines() if not line.startswith('**当前图形已于2026-11-05'))
        path.write_text(note+text)


def render_reference_figures(cfg,lineages):
    # A–D is a joint deliverable; prepare both lineages to avoid mixed versions.
    loaded={}
    for lin in ['fibro','endo']:
        search.prepare(cfg,lin,resume=True)
        loaded[lin]=load_inputs(cfg,lin)
    make_figures(cfg,loaded);write_figure_report(cfg,loaded)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.yaml'))
    p.add_argument('--lineage',choices=['fibro','endo','both'],default='both');a=p.parse_args()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    render_reference_figures(sp.read_config(a.config),['fibro','endo'] if a.lineage=='both' else [a.lineage])
