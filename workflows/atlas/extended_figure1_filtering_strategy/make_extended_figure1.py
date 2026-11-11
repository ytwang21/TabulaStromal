#!/usr/bin/env python
"""Extended Figure 1: real stored AveExp filtering evaluations and UMAPs.

Four dataset blocks in a 2x2 grid: TISCH 6 datasets/block; HUSCH 8/block.
HUSCH width / TISCH width is exactly 8/6; both figures have the same height.
"""
import argparse
import colorsys
import gc
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import pickle
import shutil
import sys
import textwrap
import zipfile

HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'tmp'/'matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR',str(HERE/'tmp'/'numba'))
Path(os.environ['MPLCONFIGDIR']).mkdir(parents=True,exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from pickle_metadata import metadata_records

SOURCES={
    'TISCH':Path('/mnt/Mercury/home/wangyuting/project/LLM/LLM_Aging/small_file/tisch_24dataset_marker_deg_aveExp.pkl'),
    'HUSCH':Path('/mnt/Mercury/home/wangyuting/project/LLM/LLM_Aging/small_file/husch_32dataset_marker_deg_aveExp.pkl'),
}
COHORTS={'TISCH':{'datasets':24,'per_block':6,'annotation':'Celltype_curated'},
         'HUSCH':{'datasets':32,'per_block':8,'annotation':'Celltype_level3'}}
FILTERS={'Marker':'Fib_marker','DEG(only)':'FibOnly_deg','DEG(like)':'FibLike_deg'}
METHODS=['DEG(only)','DEG(like)','Marker']
METHOD_COLORS={'DEG(only)':'#1F77B4','DEG(like)':'#FF7F0E','Marker':'#2CA02C'}
METHOD_STYLES={'DEG(only)':'-','DEG(like)':'--','Marker':'-.'}
BINARY_COLORS={'Fibroblast':'#1F77B4','None':'#FF7F0E'}
METRICS=['Recall','Accuracy','Precision','F1']
SEED=42


def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False))


def source_fingerprints(paths):
    return {k:dict(path=str(v.resolve()),bytes=v.stat().st_size,mtime_ns=v.stat().st_mtime_ns)
            for k,v in paths.items()}


def clean(series):
    return series.astype('string').fillna('Unknown').astype(str)


def sample_rows(obs,maximum,seed):
    """Same stratified IDs in all four UMAP rows; retain every label stratum."""
    n=len(obs)
    if maximum==0 or maximum>=n:return np.arange(n)
    codes,_=pd.factorize(pd.MultiIndex.from_frame(obs),sort=True)
    counts=np.bincount(codes)
    if maximum<len(counts):raise ValueError('Display cap is smaller than the number of annotation/filter strata')
    raw=counts*maximum/n
    quota=np.minimum(counts,np.maximum(1,np.floor(raw).astype(int)))
    while quota.sum()>maximum:
        candidates=np.where(quota>1)[0]
        quota[candidates[np.argmax((quota-raw)[candidates])]]-=1
    while quota.sum()<maximum:
        candidates=np.where(quota<counts)[0]
        quota[candidates[np.argmax((raw-quota)[candidates])]]+=1
    rng=np.random.default_rng(seed)
    rows=np.concatenate([rng.choice(np.flatnonzero(codes==i),size=int(q),replace=False)
                         for i,q in enumerate(quota)])
    return np.sort(rows)


def prepare(paths,out,maximum,full_pickle=False,refresh=False):
    cache=out/'metadata_cache.pkl'
    key=dict(sources=source_fingerprints(paths),plot_cells=maximum,seed=SEED,version=1)
    if cache.exists() and not refresh:
        with cache.open('rb') as f:saved=pickle.load(f)
        if saved['key']!=key:raise ValueError('Input/cap changed; rerun with --refresh to regenerate the metadata cache')
        logging.info('Reusing verified metadata-only plotting cache')
        return saved['cohorts'],key
    cohorts={};audit=[];metrics=[];reference_counts=[];selection_counts=[]
    cells_dir=out/'plot_cells';cells_dir.mkdir(exist_ok=True)
    for cohort,path in paths.items():
        logging.info('Extracting %s plotting metadata from %s',cohort,path)
        cfg=COHORTS[cohort];records=[]
        with metadata_records(path,full=full_pickle) as (raw_records,reader):
            # The actual files are lists of dictionaries containing recall_df,
            # adata and proportion. DataFrame/dictionary adapters are in
            # normalize_records() in pickle_metadata.py for other inputs.
            if len(raw_records)!=cfg['datasets']:raise ValueError(f'{cohort}: expected {cfg["datasets"]} records')
            for i,raw in enumerate(raw_records):
                if not isinstance(raw,dict) or not {'recall_df','adata'}<=set(raw):
                    raise ValueError('Each record must provide recall_df and adata; inspect the input structure here')
                frame=pd.DataFrame(raw['recall_df']).copy()
                required={'Celltype','Marker','dataset','thresh',*METRICS}
                if not required<=set(frame):raise ValueError(f'Missing evaluation columns: {required-set(frame)}')
                names=frame['dataset'].dropna().astype(str).unique()
                if len(names)!=1:raise ValueError('Each dataset record must have one evaluation dataset identifier')
                dataset=names[0]
                adata=raw['adata'];source_obs=adata.obs
                cols=[cfg['annotation'],*FILTERS.values()]
                if not set(cols)<=set(source_obs):raise ValueError(f'{dataset}: missing obs fields {set(cols)-set(source_obs)}')
                obs=source_obs[cols].copy()
                for col in cols:obs[col]=clean(obs[col])
                obs=obs.rename(columns={cfg['annotation']:'Reference'})
                for col in FILTERS.values():
                    if set(obs[col])-set(BINARY_COLORS):raise ValueError(f'{dataset}: unexpected binary labels in {col}: {set(obs[col])}')
                evaluation=frame.loc[frame.Celltype=='Fibroblast'].copy()
                for col in [*METRICS,'thresh']:evaluation[col]=pd.to_numeric(evaluation[col],errors='raise')
                if evaluation.duplicated(['Marker','thresh']).any():raise ValueError('Repeated method/threshold rows require review')
                if set(evaluation.Marker)!=set(METHODS):raise ValueError(f'{dataset}: expected the three stored Fibroblast methods')
                values=evaluation[METRICS].to_numpy()
                if np.nanmin(values)<-1e-8 or np.nanmax(values)>1+1e-8:raise ValueError('Evaluation metric outside [0,1]')
                if 'X_umap' not in adata.obsm:raise ValueError(f'{dataset}: no saved X_umap')
                xy=reader.array(adata.obsm['X_umap']) if reader else np.asarray(adata.obsm['X_umap'])
                if xy.shape!=(len(obs),2):raise ValueError(f'{dataset}: UMAP/obs shape mismatch {xy.shape}')
                valid=np.isfinite(xy).all(axis=1)
                if not valid.all():raise ValueError(f'{dataset}: non-finite original UMAP coordinates require review')
                rows=sample_rows(obs,maximum,SEED+i)
                displayed=obs.iloc[rows].copy()
                displayed.insert(0,'source_row',rows)
                displayed.insert(1,'cell_id',obs.index.to_numpy()[rows].astype(str))
                displayed['UMAP1']=xy[rows,0];displayed['UMAP2']=xy[rows,1]
                displayed.to_csv(cells_dir/f'{cohort}_{i+1:02d}.csv.gz',index=False)
                counts={method:int(obs[field].eq('Fibroblast').sum()) for method,field in FILTERS.items()}
                rc=obs.Reference.value_counts().to_dict()
                record=dict(dataset=dataset,index=i,n_cells=len(obs),n_plot_cells=len(rows),
                    obs=displayed,xy=np.array(xy[rows],copy=True),
                    bounds=np.array([xy.min(axis=0),xy.max(axis=0)]),
                    evaluation=evaluation,selected_counts=counts,reference_counts=rc,
                    annotation_field=cfg['annotation'],original_obs_columns=list(source_obs.columns))
                records.append(record)
                audit.append(dict(cohort=cohort,source_position=i+1,block='ABCD'[i//cfg['per_block']],
                    dataset=dataset,reference_field=cfg['annotation'],n_cells=len(obs),n_plot_cells=len(rows),
                    n_reference_types=len(rc),**{f'selected_{m}':v for m,v in counts.items()},
                    thresholds=';'.join(map(str,sorted(evaluation.thresh.unique())))))
                metrics.append(frame.assign(cohort=cohort,source_position=i+1))
                reference_counts.extend(dict(cohort=cohort,dataset=dataset,celltype=k,n_cells=int(v)) for k,v in rc.items())
                selection_counts.extend(dict(cohort=cohort,dataset=dataset,method=k,n_selected=v,n_cells=len(obs),
                                              fraction_selected=v/len(obs)) for k,v in counts.items())
                logging.info('%s %02d %s: %s cells; %s displayed',cohort,i+1,dataset,len(obs),len(rows))
            if len(set(r['dataset'] for r in records))!=len(records):raise ValueError('Duplicate dataset identifiers')
            cohorts[cohort]=records
        del raw_records,raw,adata,source_obs,obs,xy,reader
        gc.collect()
    pd.DataFrame(audit).to_csv(out/'dataset_audit.csv',index=False)
    pd.concat(metrics,ignore_index=True).to_csv(out/'stored_evaluation_all_lineages.csv',index=False)
    pd.DataFrame(reference_counts).to_csv(out/'full_reference_cell_counts.csv',index=False)
    pd.DataFrame(selection_counts).to_csv(out/'full_saved_filter_counts.csv',index=False)
    with cache.open('wb') as f:pickle.dump(dict(key=key,cohorts=cohorts),f,protocol=4)
    dump(out/'data_provenance.json',key)
    return cohorts,key


def reference_palette(cohorts):
    labels=sorted({label for records in cohorts.values() for r in records for label in r['reference_counts']})
    colors={label:matplotlib.colors.to_hex(colorsys.hls_to_rgb((i*.61803398875)%1,.56,.57))
            for i,label in enumerate(labels)}
    for label in labels:
        canonical=label.lower().replace(' ','')
        if canonical in ['fibroblast','fibroblasts']:colors[label]='#D64141'
        elif canonical in ['endothelial','endothelialcells']:colors[label]='#A7D641'
        elif canonical=='mono/macro':colors[label]='#4147D6'
        elif canonical in ['unknown','nan','none']:colors[label]='#A7ABB1'
    return colors


def curve(ax,record,metric,column):
    table=record['evaluation']
    for method in METHODS:
        part=table[table.Marker==method].sort_values('thresh')
        ax.plot(part.thresh,part[metric],color=METHOD_COLORS[method],ls=METHOD_STYLES[method],
                lw=1.4,marker='o',ms=2.4,rasterized=False)
    thresholds=sorted(table.thresh.unique())
    ax.set_xlim(min(thresholds)-.025,max(thresholds)+.025);ax.set_ylim(-.025,1.05)
    ax.set_xticks([thresholds[0],thresholds[len(thresholds)//2],thresholds[-1]])
    ax.set_yticks([0,.5,1]);ax.tick_params(labelsize=7,length=2,width=.6,pad=2)
    ax.set_xlabel('Threshold',fontsize=8,labelpad=2)
    ax.set_ylabel('F1 score' if metric=='F1' else metric,fontsize=9,labelpad=4)
    ax.grid(color='#DDDFE1',lw=.5,zorder=0)
    ax.spines[['top','right']].set_visible(False)
    for spine in ax.spines.values():spine.set_linewidth(.6)
    if column==0:ax.yaxis.label.set_fontweight('normal')


def umap(ax,record,field,palette,seed):
    xy=record['xy'];values=record['obs'][field].to_numpy()
    order=np.random.default_rng(seed).permutation(len(xy))
    if field=='Reference':
        ax.scatter(xy[order,0],xy[order,1],c=[palette[v] for v in values[order]],
                   s=.60,alpha=.93,edgecolors='none',rasterized=False)
        ax.set_title('Reference annotation',fontsize=9,pad=3)
    else:
        for category in ['None','Fibroblast']:
            selected=order[values[order]==category]
            ax.scatter(xy[selected,0],xy[selected,1],color=BINARY_COLORS[category],s=.60,
                       alpha=.93,edgecolors='none',rasterized=False)
        method=next(k for k,v in FILTERS.items() if v==field)
        ax.set_title(f'{method}: {record["selected_counts"][method]:,} selected',fontsize=9,pad=3)
    low,high=record['bounds'];span=np.maximum(high-low,1e-8)
    ax.set_xlim(low[0]-.03*span[0],high[0]+.03*span[0])
    ax.set_ylim(low[1]-.03*span[1],high[1]+.03*span[1])
    ax.set_aspect('equal',adjustable='datalim')
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_xlabel('UMAP-1',fontsize=7.5,labelpad=2);ax.set_ylabel('UMAP-2',fontsize=7.5,labelpad=2)
    for spine in ax.spines.values():spine.set_linewidth(.55);spine.set_color('#6F7479')


def section_band(fig,axes,rows,label,color):
    positions=[axes[row,0].get_position() for row in rows]
    lower=min(p.y0 for p in positions);upper=max(p.y1 for p in positions)
    x=positions[0].x0-.65/fig.get_figwidth();width=.28/fig.get_figwidth()
    fig.add_artist(FancyBboxPatch((x,lower),width,upper-lower,boxstyle='round,pad=0.001',
        transform=fig.transFigure,facecolor=color,edgecolor='none',zorder=-1))
    fig.text(x+width/2,(lower+upper)/2,label,ha='center',va='center',rotation=90,fontsize=10)


def make_figure(cohort,records,palette,width,height,maximum):
    per_block=COHORTS[cohort]['per_block']
    # Four macro panels; each macro panel contains 8 rows of genuine results.
    fig,macro_axes=plt.subplots(2,2,figsize=(width,height),squeeze=False)
    fig.subplots_adjust(left=.04,right=.985,bottom=.025,top=.94,wspace=.10,hspace=.15)
    fig.text(.04,.992,f'Extended Figure 1 | {cohort} | Fibroblast filtering strategy',
             fontsize=24,va='top')
    cap='all cells' if maximum==0 else f'up to {maximum:,} stratified cells/dataset'
    fig.text(.04,.969,f'{len(records)} datasets; 4 blocks x {per_block} datasets; stored AveExp evaluation; '
             f'UMAP: {cap}; metrics/counts: all source cells',fontsize=11,va='top')
    handles=[Line2D([],[],color=METHOD_COLORS[m],ls=METHOD_STYLES[m],lw=1.8,label=f'{m} (curve)') for m in METHODS]
    handles += [Line2D([],[],ls='',marker='o',ms=6,color=c,label=f'{k} (UMAP)') for k,c in BINARY_COLORS.items()]
    fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.985,.994),ncol=5,
               frameon=False,fontsize=10,handlelength=1.5,columnspacing=1.0)
    group_reference_counts=[]
    for block,macro in enumerate(macro_axes.flat):
        macro.set_axis_off()
        nested=macro.get_subplotspec().subgridspec(10,per_block,
            height_ratios=[.55,1,1,1,1,1.15,1.15,1.15,1.15,1.20],hspace=.30,wspace=.30)
        header=fig.add_subplot(nested[0,:]);header.set_axis_off()
        start=block*per_block;batch=records[start:start+per_block]
        header.text(0,.60,'ABCD'[block],fontsize=28,ha='left',va='center')
        header.text(.055,.60,f'Datasets {start+1:02d}–{start+per_block:02d}',fontsize=16,ha='left',va='center')
        axes=np.empty((8,per_block),dtype=object)
        for row in range(8):
            for col,record in enumerate(batch):
                ax=fig.add_subplot(nested[row+1,col]);axes[row,col]=ax
                if row<4:
                    curve(ax,record,METRICS[row],col)
                    if row==0:ax.set_title(textwrap.fill(record['dataset'],width=24,break_long_words=True),fontsize=10,pad=6)
                elif row==4:umap(ax,record,'Reference',palette,SEED+record['index'])
                else:umap(ax,record,list(FILTERS.values())[row-5],palette,SEED+record['index'])
        section_band(fig,axes,range(4),'AveExp evaluation','#E5E3E5')
        section_band(fig,axes,[4],'Reference','#E7EDB0')
        section_band(fig,axes,range(5,8),'Saved filtering','#BFDFDC')
        legend_ax=fig.add_subplot(nested[9,:]);legend_ax.set_axis_off()
        labels=sorted({v for record in batch for v in record['reference_counts']})
        group_reference_counts.append(len(labels))
        items=[Line2D([],[],ls='',marker='o',ms=4,color=palette[v],label=v) for v in labels]
        legend_ax.legend(handles=items,loc='center',ncol=per_block,frameon=False,fontsize=8,
            title='Reference cell types',title_fontsize=10,handletextpad=.4,handlelength=.8,
            columnspacing=1.1,labelspacing=.35,borderaxespad=0)
    return fig,group_reference_counts


def save(fig,stem,dpi):
    for artist in fig.findobj():
        if hasattr(artist,'set_rasterized'):artist.set_rasterized(False)
    logging.info('Saving %s (%s x %s inches)',stem.name,*fig.get_size_inches())
    # Fixed full canvas: no tight cropping, so width ratio and equal heights
    # are preserved exactly in both PNG pixels and PDF page dimensions.
    # Agg truncates floating-point pixel sizes: 53.333... * 300 can become
    # 15999.999... . Round to the intended pixel grid and nudge only the
    # raster export upward by one floating-point step, preserving the exact
    # physical figure size for the vector PDF.
    physical_size=fig.get_size_inches().copy()
    pixel_size=np.rint(physical_size*dpi).astype(int)
    fig.set_size_inches(np.nextafter(pixel_size/dpi,np.inf),forward=False)
    try:
        fig.savefig(stem.with_suffix('.png'),dpi=dpi,bbox_inches=None,facecolor='white')
    finally:
        fig.set_size_inches(physical_size,forward=False)
    fig.savefig(stem.with_suffix('.pdf'),dpi=dpi,bbox_inches=None,facecolor='white')
    plt.close(fig)


def validate(out,dpi,width,height):
    from PIL import Image
    from pypdf import PdfReader
    Image.MAX_IMAGE_PIXELS=int(width*8/6*height*dpi*dpi*1.25)
    results={}
    for cohort,factor in [('TISCH',1),('HUSCH',8/6)]:
        stem=out/f'Extended_Figure1_{cohort}'
        with Image.open(stem.with_suffix('.png')) as im:
            size=im.size;im.verify()
        expected=(round(width*factor*dpi),round(height*dpi))
        if size!=expected:raise ValueError(f'{cohort}: unexpected PNG dimensions {size}, expected {expected}')
        reader=PdfReader(stem.with_suffix('.pdf'));images=[];fonts=[];embedded=[]
        if len(reader.pages)!=1:raise ValueError('Expected one page per cohort')
        page=reader.pages[0]
        def inspect(resources):
            resources=resources.get_object()
            for ref in resources.get('/Font',{}).get_object().values() if '/Font' in resources else []:
                font=ref.get_object();fonts.append(str(font.get('/Subtype')))
                for child in font.get('/DescendantFonts',[font]):
                    child=child.get_object();descriptor=child.get('/FontDescriptor')
                    if descriptor:embedded.append(any(k in descriptor.get_object() for k in ['/FontFile','/FontFile2','/FontFile3']))
            for ref in resources.get('/XObject',{}).get_object().values() if '/XObject' in resources else []:
                obj=ref.get_object()
                if obj.get('/Subtype')=='/Image':images.append(True)
                if '/Resources' in obj:inspect(obj['/Resources'])
        inspect(page['/Resources'])
        if images or '/Type3' in fonts or not embedded or not all(embedded):raise ValueError('PDF vector/font check failed')
        page_size=(float(page.mediabox.width),float(page.mediabox.height))
        if not np.allclose(page_size,[width*factor*72,height*72],atol=1e-5):raise ValueError('Unexpected PDF page size')
        results[cohort]=dict(png_pixels=size,pdf_points=page_size,image_xobjects=0,fonts_embedded=True)
    ratio=results['HUSCH']['pdf_points'][0]/results['TISCH']['pdf_points'][0]
    if not np.isclose(ratio,8/6,atol=1e-10):raise ValueError('Width ratio differs from 8/6')
    report=dict(status='passed',cohorts=results,width_ratio=ratio,equal_heights=True,
                blocks_per_cohort=4,grid=[2,2],datasets_per_block={'TISCH':6,'HUSCH':8})
    dump(out/'validation_report.json',report)
    return report


def documents(out,key,cohorts,width,height,dpi,maximum):
    import nbformat as nbf
    from matplotlib import font_manager
    fonts=out/'fonts';fonts.mkdir(exist_ok=True)
    font=Path(font_manager.findfont('DejaVu Sans'));shutil.copy2(font,fonts/font.name)
    env={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','anndata','matplotlib','pypdf','Pillow','nbformat']}
    dump(out/'environment.json',dict(python=sys.version,packages=env))
    (HERE/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in env.items())+'\n')
    (out/'README.md').write_text(f'''# Extended Figure 1: filtering strategy

实际完成 TISCH 24 套与 HUSCH 32 套数据的绘图。两张图都分 A–D 四个大面板，采用 2×2 网格；TISCH 每面板 6 套，HUSCH 每面板 8 套。原始记录顺序保留，见 dataset_audit.csv。

TISCH 画布 {width:g}×{height:g} 英寸；HUSCH {width*8/6:g}×{height:g} 英寸。PNG 为 {dpi} dpi，PDF 完全矢量并嵌入 TrueType 字体，未将 PNG 嵌入 PDF。固定画布、不使用独立 tight 裁切，因而实际宽度比例精确为 8/6、高度相同。验证结果见 validation_report.json。

## 实际数据内容与图行

两份 pkl 的实际结构都是 list[dict]，每条记录含 recall_df（DataFrame）、adata（AnnData）和 proportion（DataFrame）。每套 Fibroblast 评估包含 Marker、DEG(only)、DEG(like) 三种方法和 0.1–0.9 的阈值。每个面板从上到下为八行：Recall、Accuracy、Precision、F1、原始参考注释 UMAP、保存的 Marker UMAP、保存的 DEG(only) UMAP、保存的 DEG(like) UMAP。

这些指定文件提供 AveExp 结果。参考图片中的 logFC、100 DEG、150 DEG 不是这两份文件中已验证的结果字段，因此没有虚构这些分析。四项曲线直接读取完整原始 recall_df，不从展示子集重新估算；三种过滤图直接使用原 obs 的 Fib_marker、FibOnly_deg、FibLike_deg，不重设阈值。图中 selected 数量为对应标签在全体源细胞中的数量。保存标签的具体操作阈值未在源对象中明确标注，图中标为 Saved filtering。

参考注释：TISCH 使用 Celltype_curated；HUSCH 使用 Celltype_level3（含 Fibroblast），该字段与该队列的精细细胞注释一致。分类标签保持原文，未重命名或合并；全局颜色字典见 reference_palette.json。相同名称在两个队列中使用同色，Fibroblast/Fibroblasts 均为红色。曲线颜色和二分类 UMAP 颜色有各自明示图例。

## 显示与统计

UMAP 每套最多显示 {maximum:,} 个固定种子 42 的分层代表细胞（--plot-cells 0 为显示全部）。分层同时考虑参考类别和三种二分类结果；同一套数据四张 UMAP 使用完全相同的 ID 和原始坐标。图中的指标、参考注释计数及过滤 selected 数量使用全部源细胞。plot_cells/*.csv.gz 给出实际绘图 ID、源行号、注释、保存标签和 UMAP 坐标。完整统计表为 full_reference_cell_counts.csv、full_saved_filter_counts.csv；原评估表（含其他原始 lineage）为 stored_evaluation_all_lineages.csv。

## 大文件读取

原始 pkl 分别约 51 GB 和 44 GB。默认使用 pickle_metadata.py 跳过绘图不需要的大型表达数组，仅提取真实表格、注释和 UMAP；不修改源文件。此过程创建的内存 AnnData 是绘图专用的不完整对象，不能用于表达分析或另存为完整 h5ad；保存的 metadata_cache.pkl 只包含普通绘图记录，不含这些不完整 AnnData。读取器通过与标准 pickle.load 的合成对象对照检查，见 loader_validation.json。

如需标准完整读取，可使用 --full-pickle --refresh（需要足够内存）。若换用其他 pkl，先确认 normalize_records() 与 required obs 字段；兼容 DataFrame 行记录、单记录字典、按数据集索引的字典和记录列表。缺少评估表、UMAP、注释或过滤字段时明确报错，不生成假数据。

## 运行

```bash
cd {HERE}
external_inputs/python -u make_extended_figure1.py
```

后台运行：`bash run_extended_fig1.sh`。绘图缓存经过源文件路径、字节数、mtime、抽样规模核验后复用。输入路径可用 --tisch 和 --husch 修改；画布可用 --tisch-width 和 --height 修改，HUSCH 始终自动取 8/6 倍宽度。--prepare-only 只提取与核查数据，不画图。原提示词代码已备份为 reference_code.txt。

Notebook 仅查看已有结果，默认不重新读取大型文件。Adobe Illustrator 编辑时如提示缺字体，可安装 fonts/DejaVuSans.ttf；本环境未实际启动 Illustrator。
''')
    nb=nbf.v4.new_notebook(metadata={'kernelspec':{'display_name':'Python 3','name':'python3','language':'python'}})
    nb.cells=[nbf.v4.new_markdown_cell('# Extended Figure 1 — filtering strategy\n\n真实 TISCH/HUSCH 数据；四个 2×2 大面板。图中数据字段、八行定义和显示抽样详见 README。'),
        nbf.v4.new_code_cell(f'from pathlib import Path\nimport pandas as pd,json\nfrom IPython.display import display,Image,Markdown\nOUT=Path({str(out)!r})\nCODE=Path({str(HERE)!r})\n'),
        nbf.v4.new_code_cell('display(pd.read_csv(OUT/"dataset_audit.csv"))\ndisplay(json.loads((OUT/"validation_report.json").read_text()))'),
        nbf.v4.new_code_cell('for cohort in ["TISCH","HUSCH"]:\n    display(Image(filename=str(OUT/f"Extended_Figure1_{cohort}.png"),width=1800))'),
        nbf.v4.new_code_cell('display(Markdown((OUT/"README.md").read_text()))'),
        nbf.v4.new_code_cell('REBUILD=False\nif REBUILD:\n    import subprocess,sys\n    subprocess.run([sys.executable,str(CODE/"make_extended_figure1.py")],check=True,cwd=CODE)')]
    nbf.validate(nb);nbf.write(nb,HERE/'Extended_Figure1.ipynb')
    src=out/'source_code';src.mkdir(exist_ok=True)
    for name in ['make_extended_figure1.py','pickle_metadata.py','verify_loader.py','run_extended_fig1.sh','requirements.txt','Extended_Figure1.ipynb']:
        shutil.copy2(HERE/name,src/name)
    archive=out/'Extended_Figure1_filtering_strategy.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for path in out.rglob('*'):
            if path.is_file() and path!=archive and path.name!='metadata_cache.pkl':z.write(path,path.relative_to(out))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tisch',type=Path,default=SOURCES['TISCH'])
    parser.add_argument('--husch',type=Path,default=SOURCES['HUSCH'])
    parser.add_argument('--output-dir',type=Path,default=HERE/'results')
    parser.add_argument('--tisch-width',type=float,default=40)
    parser.add_argument('--height',type=float,default=42)
    parser.add_argument('--dpi',type=int,default=300)
    parser.add_argument('--plot-cells',type=int,default=12000)
    parser.add_argument('--full-pickle',action='store_true')
    parser.add_argument('--refresh',action='store_true')
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    if args.dpi<1 or args.tisch_width<=0 or args.height<=0 or args.plot_cells<0:parser.error('Invalid size/dpi/display cap')
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(out/'run.log')],force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'text.color':'black',
        'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
        'pdf.fonttype':42,'ps.fonttype':42,'savefig.bbox':None,'savefig.facecolor':'white'})
    cohorts,key=prepare({'TISCH':args.tisch,'HUSCH':args.husch},out,args.plot_cells,args.full_pickle,args.refresh)
    palette=reference_palette(cohorts);dump(out/'reference_palette.json',palette)
    if args.prepare_only:
        logging.info('Preparation complete: 24 TISCH and 32 HUSCH datasets')
        return
    for cohort,factor in [('TISCH',1),('HUSCH',8/6)]:
        fig,counts=make_figure(cohort,cohorts[cohort],palette,args.tisch_width*factor,args.height,args.plot_cells)
        logging.info('%s reference categories per block: %s',cohort,counts)
        save(fig,out/f'Extended_Figure1_{cohort}',args.dpi)
    validate(out,args.dpi,args.tisch_width,args.height)
    documents(out,key,cohorts,args.tisch_width,args.height,args.dpi,args.plot_cells)
    logging.info('SUCCESS: %s',out)


if __name__=='__main__':main()
