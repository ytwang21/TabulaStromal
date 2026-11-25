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


def inch_axes(fig,left,bottom,width,height):
    fw,fh=fig.get_size_inches()
    return fig.add_axes([left/fw,bottom/fh,width/fw,height/fh])


def inch_text(fig,x,y,text,**kwargs):
    fw,fh=fig.get_size_inches()
    return fig.text(x/fw,y/fh,text,**kwargs)


def geometry():
    """Shared physical dimensions; panels consume space according to content."""
    g={'width':32.,'height':36.8,'left':.72,'right':.46,'lineage_gap':.95,
       'umap_gap_x':.36,'umap_gap_y':.42,'colorbar_width':.11,
       'colorbar_gap':.05,'colorbar_tick_space':.28,
       'bar_gap_x':.35,'bar_gap_y':.48,'bar_bottom':1.12,'legend_y':.10}
    # 2+3 UMAP columns and 7+6 bar-chart columns on independent horizontal grids.
    usable=g['width']-g['left']-g['right']-g['lineage_gap']
    g['umap_side']=(usable-3*g['umap_gap_x'])/5-g['colorbar_width']-g['colorbar_gap']-g['colorbar_tick_space']
    g['umap_tile_width']=g['umap_side']+g['colorbar_width']+g['colorbar_gap']+g['colorbar_tick_space']
    g['x_unit_inches']=(usable-11*g['bar_gap_x'])/(7*4+6*6)
    g['umap_widths']={lin:n*g['umap_tile_width']+(n-1)*g['umap_gap_x'] for lin,n in [('fibro',2),('endo',3)]}
    g['bar_widths']={lin:n*k*g['x_unit_inches']+(n-1)*g['bar_gap_x'] for lin,n,k in [('fibro',7,4),('endo',6,6)]}
    g['umap_region_height']=3*g['umap_side']+2*g['umap_gap_y']
    g['bar_height']=(g['umap_region_height']-7*g['bar_gap_y'])/8
    g['umap_top']=g['height']-.95
    g['umap_header_y']=g['height']-.33
    g['bar_header_y']=g['bar_bottom']+g['umap_region_height']+.47
    return g


def umap_panel(fig,cfg,lin,obs,left,top,g):
    rows=np.load(sp.cache(cfg,lin)/'plot_rows.npy')
    xy=np.load(sp.cache(cfg,lin)/'umap.npy',mmap_mode='r')[rows]
    vals=select_umap_signatures(cfg,lin,rows)
    finite=np.concatenate([v[np.isfinite(v)] for v in vals.values()])
    vmax=max(1.,float(math.ceil(np.quantile(finite,.995)))) if len(finite) else 1.
    norm=Normalize(0,vmax);ncols=2 if lin=='fibro' else 3;side=g['umap_side']
    # Square padding preserves coordinates and isotropic aspect, while making
    # the physical axes boxes identical between the two differently shaped UMAPs.
    low=np.nanmin(xy,axis=0);high=np.nanmax(xy,axis=0)
    center=(low+high)/2;radius=float(np.max(high-low))*.53
    for i,key in enumerate(UMAP_LAYOUT[lin]):
        x=left+(i%ncols)*(g['umap_tile_width']+g['umap_gap_x'])
        y=top-(i//ncols+1)*side-(i//ncols)*g['umap_gap_y']
        ax=inch_axes(fig,x,y,side,side)
        cax=inch_axes(fig,x+side+g['colorbar_gap'],y,g['colorbar_width'],side)
        scatter=ax.scatter(*xy.T,c=vals[key],s=.65,cmap='viridis',norm=norm,linewidths=0,rasterized=True)
        cb=fig.colorbar(scatter,cax=cax,extend='max');cb.ax.tick_params(labelsize=9,length=2,pad=2)
        cb.locator=MaxNLocator(nbins=4);cb.update_ticks()
        if cb.solids is not None:cb.solids.set_rasterized(False)
        ax.set_title(key,fontsize=15,pad=4,fontweight='normal')
        ax.set_xlim(center[0]-radius,center[0]+radius);ax.set_ylim(center[1]-radius,center[1]+radius)
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
        ax._signature_geometry={'kind':'umap','lineage':lin,'signature':key}
    sp.dump(sp.root(cfg)/f'data/plot_reference_umap_display_{lin}.json',
            {'color_vmin':0,'color_vmax':vmax,'color_strategy':'shared within lineage, ceil pooled99.5percentile; display saturation only',
             'n_sampled_cells':len(rows),'n_atlas_cells':len(obs),'panel_layout':UMAP_LAYOUT[lin],
             'grid_rows':3,'grid_columns':ncols,'physical_axes_inches':[side,side],
             'coordinate_display':'original coordinates, isotropic square padding; no deformation',
             'subtype_score_source':str(sp.cache(cfg,lin)/'figure_mean.npy'),
             'method':'mean log expression; no clipping or annotation-dependent score changes'})


def format_value(value):
    if not np.isfinite(value):return 'NA'
    if value==0:return '0'
    if abs(value)<.001:return f'{value:.1e}'
    if abs(value)<.1:return f'{value:.3f}'
    return f'{value:.2f}'


def bar_panel(fig,cfg,lin,table,left,top,g):
    order=BAR_ORDER[lin];pal=sp.read_json(sp.root(cfg)/'audit/fig1_v2_colors.json')['general_colors']
    ncols=7 if lin=='fibro' else 6;axwidth=len(order)*g['x_unit_inches']
    values=table[[f'mean_{s}' for s in order]].to_numpy()
    if np.nanmin(values)<0:raise ValueError('Negative mean-log score; refuse to hide negative scores')
    ylim=float(np.nanmax(values))*1.25
    for i,row in table.iterrows():
        x=left+(i%ncols)*(axwidth+g['bar_gap_x'])
        y=top-(i//ncols+1)*g['bar_height']-(i//ncols)*g['bar_gap_y']
        ax=inch_axes(fig,x,y,axwidth,g['bar_height']);v=values[i]
        ax.set_facecolor(pale(pal[row['Celltype_general']]))
        bars=ax.bar(np.arange(len(order)),v,color=[pal[s] for s in order],width=.69,edgecolor='#333333',linewidth=.45)
        best=float(np.nanmax(v))
        for j,(bar,val) in enumerate(zip(bars,v)):
            top_value=np.isfinite(val) and np.isclose(val,best,atol=1e-10,rtol=0)
            if top_value:bar.set_linewidth(1.1)
            if np.isfinite(val):
                text=format_value(val)
                ax.text(j,val+ylim*.022,text,ha='center',va='bottom',fontsize=11 if top_value else (7 if len(text)>5 else 9),color='black')
            else:ax.text(j,ylim*.02,'NA',ha='center',fontsize=8)
        suffix='' if row['annotation_matches_highest'] else ' †'
        if row['n_Celltype_general']>1:suffix+=' [mixed]'
        title=row['panel_title']+suffix
        ax.set_title(title,fontsize=9 if len(title)<36 else 8,pad=5,fontweight='normal')
        # Identical physical x-unit size gives exactly equal bar width across
        # four- and six-signature panels. D axes are 6/4 times as wide as C axes.
        ax.set_ylim(0,ylim);ax.set_xlim(-.5,len(order)-.5)
        ax.set_xticks([]);ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis='y',labelsize=9,length=2,pad=2)
        if i%ncols==0:ax.set_ylabel('Signature score',fontsize=10,labelpad=3)
        ax.grid(axis='y',lw=.45,alpha=.27);ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(.6);ax.spines['bottom'].set_linewidth(.6)
        ax._signature_geometry={'kind':'barplot','lineage':lin,'leiden':str(row['leiden'])}
    fw,fh=fig.get_size_inches();width=g['bar_widths'][lin]
    for handles,x,title in [([Patch(facecolor=pal[s],label=s) for s in order],left+width*.25,'marker_list'),
                            ([Patch(facecolor=pale(pal[s]),label=s) for s in order],left+width*.75,'Celltype')]:
        legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(x/fw,g['legend_y']/fh),ncol=len(order),
                          fontsize=12,title=title,title_fontsize=13,frameon=True,fancybox=True,
                          handlelength=1.1,columnspacing=1.0,framealpha=1,borderpad=.6)
        legend.get_frame().set_edgecolor('#777777');legend.get_frame().set_linewidth(.7)
        legend._signature_geometry={'kind':'legend','lineage':lin,'title':title,'n_items':len(order)}


def validate_geometry(fig,label):
    """Measure final rendered axes and bar transforms, not just grid settings."""
    fig.canvas.draw();fw,fh=fig.get_size_inches();umaps=[];charts=[];legends=[]
    for ax in fig.axes:
        meta=getattr(ax,'_signature_geometry',None)
        if not meta:continue
        pos=ax.get_position();row={**meta,'width_inches':pos.width*fw,'height_inches':pos.height*fh,
                                  'left_inches':pos.x0*fw,'bottom_inches':pos.y0*fh}
        if meta['kind']=='umap':umaps.append(row)
        else:
            row['bar_width_inches']=[abs(ax.transData.transform((b.get_x()+b.get_width(),0))[0]-ax.transData.transform((b.get_x(),0))[0])/fig.dpi for b in ax.patches]
            charts.append(row)
    for legend in fig.legends:
        meta=getattr(legend,'_signature_geometry',None)
        if not meta:continue
        labels=legend.get_texts();centers=[t.get_window_extent(fig.canvas.get_renderer()).y0 for t in labels]
        assert np.ptp(centers)<1e-7,'Legend entries wrap to multiple rows'
        box=legend.get_window_extent(fig.canvas.get_renderer())
        assert box.x0>=0 and box.x1<=fw*fig.dpi and box.y0>=0,'Legend outside canvas'
        legends.append({**meta,'n_columns':legend._ncols,'single_row':True})
    assert umaps and charts and legends
    np.testing.assert_allclose([r['width_inches'] for r in umaps],umaps[0]['width_inches'],rtol=0,atol=1e-9)
    np.testing.assert_allclose([r['height_inches'] for r in umaps],umaps[0]['height_inches'],rtol=0,atol=1e-9)
    np.testing.assert_allclose([r['height_inches'] for r in charts],charts[0]['height_inches'],rtol=0,atol=1e-9)
    bars=[v for r in charts for v in r['bar_width_inches']]
    np.testing.assert_allclose(bars,bars[0],rtol=0,atol=1e-9)
    dims={lin:next(r['width_inches'] for r in charts if r['lineage']==lin) for lin in dict.fromkeys(r['lineage'] for r in charts)}
    if len(dims)==2:np.testing.assert_allclose(dims['endo']/dims['fibro'],1.5,rtol=0,atol=1e-9)
    umap_span=max(r['bottom_inches']+r['height_inches'] for r in umaps)-min(r['bottom_inches'] for r in umaps)
    bar_span=max(r['bottom_inches']+r['height_inches'] for r in charts)-min(r['bottom_inches'] for r in charts)
    np.testing.assert_allclose(bar_span,umap_span,rtol=0,atol=1e-9)
    logging.info('%s geometry verified: UMAP %.4f × %.4f in; bars %.4f in; subplot height %.4f in',label,umaps[0]['width_inches'],umaps[0]['height_inches'],bars[0],charts[0]['height_inches'])
    return {'passed':True,'figure_inches':[fw,fh],'n_UMAPs':len(umaps),'n_Leiden_panels':len(charts),
            'UMAP_width_inches':umaps[0]['width_inches'],'UMAP_height_inches':umaps[0]['height_inches'],
            'bar_width_inches':bars[0],'Leiden_panel_height_inches':charts[0]['height_inches'],
            'UMAP_region_height_inches':umap_span,'Leiden_region_height_inches':bar_span,
            'Leiden_to_UMAP_region_height_ratio':bar_span/umap_span,
            'Leiden_panel_width_inches':dims,'legends':legends,'UMAPs':umaps,'Leiden_panels':charts}


def make_figures(cfg,loaded):
    plt=sp.plotting_setup();r=sp.root(cfg);g=geometry();metrics={}
    nr=max(math.ceil(len(loaded[lin][1])/(7 if lin=='fibro' else 6)) for lin in loaded)
    bar_top=g['bar_bottom']+nr*g['bar_height']+(nr-1)*g['bar_gap_y']
    for which in ['integrated','fibro','endo']:
        lineages=['fibro','endo'] if which=='integrated' else [which]
        width=g['width'] if which=='integrated' else max(g['umap_widths'][which],g['bar_widths'][which])+g['left']+g['right']
        fig=plt.figure(figsize=(width,g['height']))
        for lin in lineages:
            obs,table=loaded[lin]
            if which=='integrated':
                ux=g['left'] if lin=='fibro' else g['left']+g['umap_widths']['fibro']+g['lineage_gap']
                bx=g['left'] if lin=='fibro' else g['left']+g['bar_widths']['fibro']+g['lineage_gap']
            else:
                ux=(width-g['umap_widths'][lin])/2;bx=(width-g['bar_widths'][lin])/2
            a,b=('A','C') if lin=='fibro' else ('B','D')
            inch_text(fig,ux-.40,g['umap_header_y'],a,fontsize=31,va='center',fontweight='normal')
            inch_text(fig,ux+g['umap_widths'][lin]/2,g['umap_header_y'],sp.LABELS[lin],fontsize=21,ha='center',va='center')
            umap_panel(fig,cfg,lin,obs,ux,g['umap_top'],g)
            inch_text(fig,bx-.40,g['bar_header_y'],b,fontsize=31,va='center',fontweight='normal')
            inch_text(fig,bx+g['bar_widths'][lin]/2,g['bar_header_y'],'Signature score of '+sp.LABELS[lin]+' subtype markers',fontsize=19,ha='center',va='center')
            bar_panel(fig,cfg,lin,table,bx,bar_top,g)
            inch_text(fig,bx+g['bar_widths'][lin]/2,.84,
                f"All Leiden clusters: {int(table.annotation_matches_highest.sum())}/{len(table)} match reference; † discordant",
                fontsize=12,ha='center',va='center')
        metrics[which]=validate_geometry(fig,which)
        path=r/'figures/integrated/ExtFig1_SignatureScore_Fibro_Endo' if which=='integrated' else r/f'figures/panels/ExtFig1_SignatureScore_{sp.LABELS[which]}'
        sp.save_figure(fig,path,cfg);plt.close(fig)
    return metrics


def write_figure_report(cfg,loaded,geometry_metrics=None):
    r=sp.root(cfg);table=pd.concat([loaded[lin][1] for lin in ['fibro','endo']],ignore_index=True)
    table.to_csv(r/'data/plot_reference_cluster_check.csv',index=False)
    lines=['# A–D layout and global marker comparison — 2026-11-24','',
        '原路径PNG/PDF已覆盖；主图为一页A–D，按上排A/B、下排C/D排列。分图沿用主图字母A/C与B/D。','',
        'A：Fibroblast ECM、AP、Myfibro、eFibro、iFibro、qFibro，3×2 UMAP。',
        'B：Endothelial AP、Tip、Stalk、ArtEC、CapEC、VenEC、EndoMT、MCs、LEC，3×3 UMAP。',
        'C/D：分别为Fibroblast/Endothelial全部Leiden；每个标题严格使用`Celltype + "_C" + leiden`。', '',
        '2026-11-24版式：A为3行×2列、B为3行×3列；所有UMAP的绘图区宽高严格相同，原坐标只作等比例方形留白。上方按2:3列数分配空间。C每行7个Leiden、D每行6个；共用同一物理x单位和小图高度，D的6柱小图宽度为C的4柱小图的1.5倍。C/D小图加高约61%，柱图区与UMAP区的总高度均为16.794英寸；画布纵向延长至36.8英寸，柱图数值与y轴标签适当放大。底部marker_list和Celltype图例中的4/6个亚型各排成单行。实测轴框、全部柱宽和图例行数保存在logs/reference_layout_update.json的physical_geometry中。', '',
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
    if geometry_metrics is None:
        previous=r/'logs/reference_layout_update.json'
        geometry_metrics=sp.read_json(previous).get('physical_geometry',{}) if previous.exists() else {}
    summary['physical_geometry']=geometry_metrics
    summary['layout_updated_date']='2026-11-24'
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
          '## 当前绘图版本：2026-11-24，单页 A–D\n\n'
          'A为Fibroblast ECM/AP/4亚型UMAP；B为Endothelial AP/Tip/Stalk/6亚型UMAP；C/D分别为全部Leiden柱状图。'
          '本图使用新的训练GSE组合选择名单 `*_subtype_markers_figure.json` 与 `figure_mean.npy`，'
          '与下方历史基线 `*_subtype_markers_final.json` / mean/module统计区分。'
          f"新名单保留GSE一致率：Fibroblast {rates['fibro']:.2%}，Endothelial {rates['endo']:.2%}；全部Leiden一致数{matched['fibro']}、{matched['endo']}。"
          '†反例保留，没有按标签加分。详见 [当前图形与组合验证报告](reference_layout_report.md)。\n\n'
          '以下正文为按原Phase0–5方法、在2026-11-21新注释对象上重新汇总验证的基线；其分数及状态与组合选择名单分开报告。\n'
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
    metrics=make_figures(cfg,loaded);write_figure_report(cfg,loaded,metrics)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.yaml'))
    p.add_argument('--lineage',choices=['fibro','endo','both'],default='both');a=p.parse_args()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    render_reference_figures(sp.read_config(a.config),['fibro','endo'] if a.lineage=='both' else [a.lineage])
