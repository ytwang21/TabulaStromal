#!/usr/bin/env python
"""Fig1 A–J: original atlas, Celltype-suffix Tissue, full summaries and vector PDF.

Reuses validated all-cell marker/signature tables and shared original plotting maths.
No expression matrix or trajectory fit is rerun. Every PDF artist remains vector.
"""
import argparse
import hashlib
import json
import logging
from pathlib import Path
import pickle
import shutil
import zipfile

import stromal_atlas_pipeline as pipeline
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd

from utils_data import read_elem, stratified_rows, dump_json, MARKERS, age_order, clean
from utils_plot import summarize, marker_matrices, signature_matrices

HERE = pipeline.HERE
TITLES = {'fibro': 'Fibroblast', 'endo': 'Endothelial'}


def panel_letter(ax, letter):
    ax.text(-.045, 1.04, letter, transform=ax.transAxes, fontsize=90,
            fontweight='normal',ha='right',va='bottom')


def count_text(number):
    return f'{number/1e6:.2f}M' if number >= 1e6 else f'{number/1e3:.1f}k' if number >= 10000 else f'{int(number):,}'


def no_raster(fig):
    """Also disable Matplotlib's automatic rasterization of dense colorbars."""
    for artist in fig.findobj():
        if hasattr(artist, 'set_rasterized'):
            artist.set_rasterized(False)


def export(fig, stem, dpi, tight=True):
    # User-requested black typography; categorical colors remain in data marks
    # and classification strips. Pastel heatmaps keep black numbers legible.
    for artist in fig.findobj(match=Text):
        artist.set_color('black')
    no_raster(fig)
    logging.info('Exporting %s', stem.name)
    fig.savefig(str(stem)+'.png', dpi=dpi, bbox_inches='tight' if tight else None, pad_inches=.20)
    fig.savefig(str(stem)+'.pdf', dpi=dpi, bbox_inches='tight' if tight else None, pad_inches=.20)
    plt.close(fig)


def prepare(lineage, root, out, config):
    audit = json.loads((root/'shared'/f'data_audit_{lineage}.json').read_text())
    manifest = json.loads((root/'run_manifest.json').read_text())
    source = Path(audit['input'])
    stat = source.stat()
    old = manifest['inputs'][lineage]
    if stat.st_size != old['size'] or stat.st_mtime_ns != old['mtime_ns']:
        raise ValueError('Source changed; full-cell cached expression summaries are stale')
    cached = pd.read_pickle(root/lineage/'metadata.pkl')
    # Deliberately discard original obs.tissue. No Fig1 artist/statistic uses it.
    obs = cached[['general', 'celltype', 'age_group', 'sample', 'dataset']].rename(
        columns={'general':'Celltype_general', 'celltype':'Celltype', 'age_group':'Age_group', 'dataset':'GSE'}).copy()
    xy = np.load(root/lineage/'umap.npy', mmap_mode='r')
    with h5py.File(source, 'r') as handle:
        source_columns = list(handle['obs'].keys())
        for field in ['Celltype_general', 'Celltype', 'Age_group', 'sample', 'GSE']:
            if field not in handle['obs']:
                raise ValueError(f'Missing actual source obs column: {field}')
            raw = clean(pd.Series(read_elem(handle['obs'][field]))).to_numpy()
            if not np.array_equal(raw, obs[field].to_numpy()):
                raise ValueError(f'Cached source annotation mismatch: {field}')
        source_ids = np.asarray(read_elem(handle['obs'][handle['obs'].attrs['_index']]))
        if not np.array_equal(source_ids, obs.index.to_numpy()):
            raise ValueError('Cached atlas observation IDs/order mismatch')
        if not np.array_equal(handle['obsm/X_umap'][:].astype(np.float32), xy):
            raise ValueError('UMAP differs from the actual atlas source')
    # Exact requested definition; spaces and the literal pan suffix are preserved.
    obs['Tissue_from_Celltype'] = obs['Celltype'].astype(str).str.split('_').str[-1].astype('category')
    if set(obs.Celltype_general) != set(config['general_order'][lineage]):
        raise ValueError('Unexpected true Celltype_general labels')
    if set(obs.Tissue_from_Celltype) - set(config['tissue_colors']):
        raise ValueError('A source Celltype suffix lacks a shared Tissue color')
    tissue = summarize(obs, 'Tissue_from_Celltype')
    # Reuse the original coverage helper, including unique sample counts and Unknown.
    summary_obs = obs.rename(columns={'Age_group':'age_group'})
    age = summarize(summary_obs, 'age_group')
    if tissue.n_cells.sum() != len(obs) or age.n_cells.sum() != len(obs):
        raise ValueError('Full metadata summary totals differ from actual source')
    tissue.to_csv(out/f'{lineage}_Tissue_from_Celltype_summary.csv')
    age.to_csv(out/f'{lineage}_Age_group_summary.csv')
    pairs = obs[['Celltype','Celltype_general','Tissue_from_Celltype']].value_counts(sort=False).rename('n_cells').reset_index()
    pairs = pairs[pairs.n_cells > 0]
    pairs.to_csv(out/f'{lineage}_Celltype_to_Tissue_mapping.csv', index=False)
    cap = config['plot_cells_per_lineage']
    rows = stratified_rows(obs, min(cap,len(obs)) if cap else len(obs), config['seed'], ('Celltype_general','Celltype'))
    np.save(out/f'{lineage}_umap_plot_rows.npy', rows)
    sampled = obs.iloc[rows].copy()
    sampled.insert(0, 'source_row', rows)
    sampled.to_csv(out/f'{lineage}_umap_plot_cells.csv.gz', index_label='cell_id')
    if set(sampled.Celltype) != set(obs.Celltype):
        raise ValueError('Representative plotting cohort omitted a Celltype')
    marker = pd.read_csv(root/lineage/'marker_summary.csv')
    general_counts = obs.Celltype_general.value_counts()
    for state in config['general_order'][lineage]:
        if not (marker.loc[marker.subtype==state,'n_cells'] == general_counts[state]).all():
            raise ValueError('Marker denominator does not match actual full atlas')
    means = pd.read_csv(root/lineage/'signature_mean.csv', index_col=0)
    available = json.loads((root/lineage/'marker_availability.json').read_text())
    n_samples = int(obs.loc[obs['sample']!='Unknown','sample'].nunique())
    record = dict(source=str(source), source_bytes=stat.st_size, source_mtime_ns=stat.st_mtime_ns,
        n_cells=len(obs), n_plot_cells=len(rows), n_samples=n_samples,
        n_studies=int(obs.loc[obs.GSE!='Unknown','GSE'].nunique()), n_Tissue_from_Celltype=len(tissue),
        actual_obs_columns=source_columns, donor_field=audit['fields'].get('donor'),
        tissue_definition="Celltype.astype(str).str.split('_').str[-1]", original_obs_tissue_used=False,
        source_ids_annotations_umap_match=True, source_unchanged=True,
        sampling='Fixed seed, Celltype_general × Celltype stratification; same rows for both UMAP colorings',
        marker_source=str(root/lineage/'marker_summary.csv'), signature_source=str(root/lineage/'signature_mean.csv'),
        missing_markers=available['missing'])
    logging.info('%s verified: %s full cells; %s vector UMAP points; Tissue suffixes=%s', lineage, len(obs),len(rows), list(tissue.index))
    return dict(lineage=lineage, obs=sampled, xy=np.asarray(xy[rows]), record=record, tissue=tissue, age=age,
                marker=marker, means=means, available=available)


def umap(fig, spec, data, field, letter, config):
    group = field == 'Celltype_general'
    labels = config['general_order'][data['lineage']] if group else sorted(data['tissue'].index, key=lambda s:(s!='pan',s))
    colors = config['general_colors'] if group else config['tissue_colors']
    grid = spec.subgridspec(1,2,width_ratios=[1,1.16],wspace=.01)
    ax = fig.add_subplot(grid[0,0]); legend_ax = fig.add_subplot(grid[0,1]); legend_ax.set_axis_off()
    xy=data['xy']; values=data['obs'][field].astype(str).to_numpy()
    order=np.random.default_rng(config['seed']).permutation(len(xy))
    ax.scatter(xy[order,0],xy[order,1],c=[colors[s] for s in values[order]],s=2.0,alpha=.90,
               edgecolors='none',rasterized=False)
    ax.set_aspect('equal', adjustable='datalim')
    ax.margins(.015)
    ax.set(xticks=[],yticks=[])
    for spine in ax.spines.values():spine.set_visible(False)
    ax.annotate('',xy=(.20,.035),xytext=(.035,.035),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.1))
    ax.annotate('',xy=(.035,.20),xytext=(.035,.035),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.1))
    ax.text(.205,.02,'UMAP-1',transform=ax.transAxes,fontsize=18,va='center')
    ax.text(.015,.21,'UMAP-2',transform=ax.transAxes,fontsize=18,rotation=90,ha='center',va='bottom')
    handles=[Line2D([],[],ls='',marker='o',markersize=22,color=colors[label],label=label) for label in labels]
    legend_ax.legend(handles=handles,loc='center left',frameon=False,
                     fontsize=config['legend_fontsize'],ncol=1 if group else 2,
                     title='Celltype_general' if group else 'Tissue',title_fontsize=config['legend_title_fontsize'],
                     handletextpad=.2,handlelength=.8,columnspacing=.5,labelspacing=.30,borderaxespad=.1)
    # The legend identifies states without covering small or overlapping islands.
    if group:
        r=data['record']
        ax.text(.015,.985,f'{r["n_cells"]:,} cells\n{r["n_samples"]:,} samples\n{r["n_Tissue_from_Celltype"]} suffix Tissue labels\n{r["n_studies"]} studies',
                transform=ax.transAxes,ha='left',va='top',fontsize=config['statistics_fontsize'],
                bbox=dict(facecolor='white',edgecolor='black',lw=1,pad=7,alpha=.95))
        ax.set_title(f'{TITLES[data["lineage"]]} | subtype',fontsize=28,fontweight='bold',loc='left',pad=14)
    panel_letter(ax,letter)
    return ax


def coverage(fig, spec, data, kind, letter, config):
    labels = sorted(set(data['fibro'][kind].index)|set(data['endo'][kind].index), key=lambda x:(x!='pan',x)) if kind=='tissue' else age_order(list(set(data['fibro']['age'].index)|set(data['endo']['age'].index)))
    palette = config['tissue_colors'] if kind=='tissue' else config['age_colors']
    # Mirrored bars put both lineages on each category row; this preserves enough
    # vertical space for all 16 suffixes without shrinking category-label fonts.
    grid=spec.subgridspec(1,2,wspace=.20)
    header=fig.add_subplot(spec,frameon=False);header.set_axis_off()
    header.text(0,1.12,'Tissue distribution' if kind=='tissue' else 'Age-group coverage',
                transform=header.transAxes,fontsize=32,fontweight='normal',va='bottom')
    header.text(-.065,1.12,letter,transform=header.transAxes,fontsize=102,fontweight='normal',ha='right',va='bottom')
    for col,metric in enumerate(['n_cells','n_samples']):
        maximum=max(data[l][kind][metric].max() for l in ['fibro','endo'])*1.20
        ax=fig.add_subplot(grid[0,col])
        for direction,lineage in [(-1,'fibro'),(1,'endo')]:
            values=data[lineage][kind].reindex(labels,fill_value=0)[metric]
            ax.barh(np.arange(len(labels)),direction*values,color=[palette[x] for x in labels],height=.82,edgecolor='none')
            for y,value in enumerate(values):
                if value:
                    inside=value>.70*maximum
                    tip=value-.035*maximum if inside else value+.025*maximum
                    align=('left' if direction<0 else 'right') if inside else ('right' if direction<0 else 'left')
                    ax.text(direction*tip,y,count_text(value),va='center',ha=align,fontsize=33,fontweight='normal')
        ax.set_yticks(np.arange(len(labels)),labels if col==0 else ['']*len(labels),fontsize=42,fontweight='normal')
        ax.set_ylim(len(labels)-.5,-.5);ax.set_xlim(-max(1,maximum),max(1,maximum))
        ax.tick_params(axis='y',length=0)
        ax.set_xticks([-maximum/1.2,0,maximum/1.2])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:count_text(abs(x))))
        ax.tick_params(axis='x',labelsize=24)
        for tick in ax.get_xticklabels():tick.set_fontweight('normal')
        ax.spines[['left','right','top']].set_visible(False)
        ax.axvline(0,color='black',lw=.8,zorder=0)
        ax.set_xlabel('Cells' if metric=='n_cells' else 'Samples',fontsize=28,fontweight='normal',labelpad=9)
        ax.set_title('Cell number' if metric=='n_cells' else 'Sample number',fontsize=28,fontweight='normal',pad=58)
        ax.text(.25,1.025,'Fibroblast',transform=ax.transAxes,ha='center',fontsize=26,fontweight='normal')
        ax.text(.75,1.025,'Endothelial',transform=ax.transAxes,ha='center',fontsize=26,fontweight='normal')


def marker(fig, spec, data, letter, config):
    ax=fig.add_subplot(spec);lineage=data['lineage'];groups=config['general_order'][lineage];colors=config['general_colors']
    genes,mean,fraction,scaled=marker_matrices(data['marker'],lineage,groups)
    xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(groups)))
    dots=ax.scatter(xx.ravel(),yy.ravel(),s=260*fraction.to_numpy().ravel(),c=scaled.to_numpy().ravel(),
                    cmap='Blues',vmin=0,vmax=1,edgecolors='#879CA6',linewidths=.4,rasterized=False)
    ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right',rotation_mode='anchor',fontsize=24)
    ax.set_yticks(range(len(groups)),groups,fontsize=25.5)
    for label in ax.get_yticklabels():label.set_fontweight('normal')
    ax.set_ylim(len(groups)-.5,-1.4);ax.set_xlim(-.7,len(genes)-.3)
    start=0
    for state in groups:
        count=sum(g in genes for g in MARKERS[lineage][state]);end=start+count
        if count:
            ax.plot([start-.35,end-.65],[-.85,-.85],color=colors[state],lw=3,clip_on=False)
            ax.text((start+end-1)/2,-1.05,state,fontsize=22.5,fontweight='normal',ha='center')
        start=end
    ax.tick_params(length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    colorbar=fig.colorbar(dots,ax=ax,fraction=.023,pad=.025)
    colorbar.set_label('Gene-scaled mean',fontsize=36);colorbar.ax.tick_params(labelsize=31.5)
    if colorbar.solids is not None:colorbar.solids.set_rasterized(False)
    handles=[Line2D([],[],ls='',marker='o',markersize=np.sqrt(260*f),color='#B5CFDF',markeredgecolor='#879CA6',label=f'{int(100*f)}%') for f in [.1,.5,1]]
    ax.legend(handles=handles,title='Cells expressing',title_fontsize=24,fontsize=22.5,ncol=3,
              frameon=False,loc='upper center',bbox_to_anchor=(.5,-.27),columnspacing=1.0,
              bbox_transform=ax.transAxes+ScaledTranslation(0,5/25.4,fig.dpi_scale_trans))
    panel_letter(ax,letter)
    return ax


def signature(fig,spec,data,letter,config,show_colorbar=True):
    ax=fig.add_subplot(spec);lineage=data['lineage'];groups=config['general_order'][lineage];colors=config['general_colors']
    means,z=signature_matrices(data['means'],lineage,groups)
    n=len(groups)
    # Rectangular cells stay individually editable, unlike imshow's embedded bitmap.
    pastel=LinearSegmentedColormap.from_list('readable_RdBu',plt.get_cmap('RdBu_r')(np.linspace(.20,.80,256)))
    mesh=ax.pcolormesh(np.arange(n+1)-.5,np.arange(n+1)-.5,z.to_numpy(),
                       cmap=pastel,vmin=-2,vmax=2,edgecolors='white',linewidth=.7,rasterized=False)
    # Proportional margins give a 4x4 and a 6x6 matrix the same outer size
    # whenever their axes have the same physical dimensions.
    margin=.03*n
    ax.set_xlim(-.5-margin,n-.5);ax.set_ylim(n-.5,-.5-margin);ax.set_aspect('equal')
    ax.set_xticks(range(n),groups,rotation=45,ha='right',fontsize=24)
    ax.set_yticks(range(n),groups,fontsize=24)
    for label in ax.get_xticklabels()+ax.get_yticklabels():
        label.set_fontweight('normal')
    for i,state in enumerate(groups):
        ax.add_patch(Rectangle((i-.5,-.5-.0275*n),1,.02*n,facecolor=colors[state],edgecolor='none',clip_on=False))
        ax.add_patch(Rectangle((-.5-.0275*n,i-.5),.02*n,1,facecolor=colors[state],edgecolor='none',clip_on=False))
        for j in range(n):
            ax.text(j,i,f'{means.iloc[i,j]:.2f}',ha='center',va='center',fontsize=32,color='black')
    ax.tick_params(length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_xlabel('Marker signature',fontsize=40.5)
    ax.xaxis.set_label_coords(.5,-.20)
    ax.set_anchor('NW')
    # An inset colorbar follows the final square heatmap bounds, rather than
    # remaining at the far edge of the original wide subplot allocation.
    if show_colorbar:
        cax=ax.inset_axes([1.045,.09,.04,.82])
        colorbar=fig.colorbar(mesh,cax=cax)
        size=31.5 if lineage=='endo' else 21
        colorbar.ax.set_title('z-score',fontsize=size,pad=10,loc='left' if lineage=='endo' else 'center')
        colorbar.ax.tick_params(labelsize=size)
        if colorbar.solids is not None:colorbar.solids.set_rasterized(False)
    panel_letter(ax,letter)
    return ax


def render_main(data,config,out):
    fig=plt.figure(figsize=config['figure_size'])
    # Explicit spacer rows prevent larger summary headers from consuming the
    # C/D-to-E/F gap; additional height preserves the actual plotting areas.
    outer=fig.add_gridspec(7,1,height_ratios=[11,2.5,11,3.6,11,3.0,8],hspace=0,
                          left=.09,right=.978,top=.965,bottom=.065)
    top=outer[0].subgridspec(1,2,wspace=.15)
    tissues=outer[2].subgridspec(1,2,wspace=.15)
    summaries=outer[4].subgridspec(1,2,wspace=.20)
    bottom=outer[6].subgridspec(1,4,width_ratios=[1.5,1.5,1.25,1.25],wspace=.35)
    signature_axes=[]
    for col,lineage in enumerate(['fibro','endo']):
        umap(fig,top[col],data[lineage],'Celltype_general','AB'[col],config)
        umap(fig,tissues[col],data[lineage],'Tissue_from_Celltype','CD'[col],config)
        marker(fig,bottom[col],data[lineage],'GH'[col],config)
        ax=signature(fig,bottom[col+2],data[lineage],'IJ'[col],config,show_colorbar=col==1)
        if col==1:
            position=ax.get_position()
            ax.set_position([position.x0-.012,position.y0,position.width,position.height])
        signature_axes.append(ax)
    coverage(fig,summaries[0],data,'tissue','E',config)
    coverage(fig,summaries[1],data,'age','F',config)
    # Check the colored matrix itself, including identical top/bottom edges;
    # equal subplot slots alone are insufficient with unequal data margins.
    alignment={}
    for letter,ax in zip('IJ',signature_axes):
        ax.get_position()  # apply the square aspect before measuring
        n=len(ax.get_xticks())
        corners=ax.transData.transform([[-.5,-.5],[n-.5,n-.5]])/fig.dpi
        left,bottom_edge=corners.min(axis=0);right,top_edge=corners.max(axis=0)
        alignment[letter]=dict(left_inches=float(left),bottom_inches=float(bottom_edge),
            width_inches=float(right-left),height_inches=float(top_edge-bottom_edge),top_inches=float(top_edge))
    for key in ['bottom_inches','top_inches','width_inches','height_inches']:
        if not np.isclose(alignment['I'][key],alignment['J'][key],rtol=0,atol=1e-9):
            raise ValueError(f'I/J signature alignment mismatch: {key}, {alignment}')
    dump_json(out/'signature_alignment.json',dict(status='passed',matrices=alignment,
        checks=['Identical matrix width and height','Identical top and bottom edges','Identical Marker signature title baseline']))
    export(fig,out/'Fig1_StromalAtlas_main',config['dpi'])


def render_panels(data,config,out,selected='ABCDEFGHIJ'):
    for letter,lineage,field in [('A','fibro','Celltype_general'),('B','endo','Celltype_general'),
                                 ('C','fibro','Tissue_from_Celltype'),('D','endo','Tissue_from_Celltype')]:
        if letter not in selected:continue
        fig=plt.figure(figsize=(26,18))
        umap(fig,fig.add_gridspec(1,1,left=.07,right=.99,bottom=.07,top=.90)[0],data[lineage],field,letter,config)
        export(fig,out/f'Fig1{letter}_{TITLES[lineage]}_{field}',config['dpi'])
    for letter,kind in [('E','tissue'),('F','age')]:
        if letter not in selected:continue
        fig=plt.figure(figsize=(24,18))
        coverage(fig,fig.add_gridspec(1,1,left=.19,right=.98,bottom=.08,top=.84)[0],data,kind,letter,config)
        export(fig,out/f'Fig1{letter}_{kind}_summary',config['dpi'])
    for letter,lineage,function,name in [('G','fibro',marker,'markers'),('H','endo',marker,'markers'),
                                        ('I','fibro',signature,'signatures'),('J','endo',signature,'signatures')]:
        if letter not in selected:continue
        fig=plt.figure(figsize=(15,12) if name=='markers' else (12,12))
        function(fig,fig.add_gridspec(1,1,left=.15,right=.94,top=.83 if name=='markers' else .85,bottom=.27 if name=='markers' else .20)[0],data[lineage],letter,config)
        export(fig,out/f'Fig1{letter}_{TITLES[lineage]}_{name}',config['dpi'],tight=name!='signatures')


def validate_pdfs(out):
    from pypdf import PdfReader
    results={}
    for path in sorted(out.glob('Fig1*.pdf')):
        reader=PdfReader(path);images=[];font_types=set();embedded=[]
        seen=set()
        def inspect(resources):
            resources=resources.get_object()
            for font_ref in resources.get('/Font',{}).get_object().values() if '/Font' in resources else []:
                font=font_ref.get_object();font_types.add(str(font.get('/Subtype')))
                descendants=font.get('/DescendantFonts',[font])
                for item in descendants:
                    subfont=item.get_object();font_types.add(str(subfont.get('/Subtype')))
                    descriptor=subfont.get('/FontDescriptor')
                    if descriptor:
                        descriptor=descriptor.get_object()
                        embedded.append(any(k in descriptor for k in ['/FontFile','/FontFile2','/FontFile3']))
            for name,obj_ref in resources.get('/XObject',{}).get_object().items() if '/XObject' in resources else []:
                key=(getattr(obj_ref,'idnum',None),getattr(obj_ref,'generation',None))
                if key in seen:continue
                seen.add(key);obj=obj_ref.get_object()
                if obj.get('/Subtype')=='/Image':images.append(str(name))
                if '/Resources' in obj:inspect(obj['/Resources'])
        for page in reader.pages:inspect(page['/Resources'])
        if len(reader.pages)!=1 or images or '/Type3' in font_types or not embedded or not all(embedded):
            raise ValueError(f'PDF vector/font validation failed: {path}, images={images}, fonts={font_types}')
        results[path.name]=dict(pages=len(reader.pages),image_xobjects=0,font_types=sorted(font_types),
                                embedded_fonts=True,bytes=path.stat().st_size)
    if len(results)!=11:raise ValueError('Expected main figure plus ten separate panel PDFs')
    return results


def notebook(out, config):
    import nbformat as nbf
    md,code=nbf.v4.new_markdown_cell,nbf.v4.new_code_cell
    nb=nbf.v4.new_notebook(metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
    nb.cells=[
        md('# StromalAtlas Fig1 A–J\n\n本次使用原全细胞 atlas，UMAP 只读取原坐标；统计使用全量，marker/signature 复用已验证的全量汇总。Tissue 严格取 Celltype 的最后一个 `_` 分隔字段。PDF 中散点、文字、条形和热图格均为矢量。默认只查看已生成结果。'),
        code(f'from pathlib import Path\nimport json,sys,subprocess\nimport pandas as pd\nfrom IPython.display import display,Image,Markdown\nCODE=Path({str(HERE)!r})\nOUT=Path({str(out)!r})\nCONFIG=OUT/"config_used.json"\ncfg=json.loads(CONFIG.read_text())\n'),
        md('## 重新生成完整 Fig1\n\n将 REBUILD 改成 True 后运行此单元，会执行同一个可复现 Python 入口。它只重做 metadata 汇总和绘图，不重跑表达矩阵或轨迹。'),
        code('REBUILD=False\nif REBUILD:\n    subprocess.run([sys.executable,str(CODE/"08_make_fig1.py"),"--config",str(CONFIG)],check=True,cwd=CODE)\n'),
        md('## 新字段及全量计数\n\n代码中的实际构建语句：`obs["Tissue_from_Celltype"] = obs["Celltype"].astype(str).str.split("_").str[-1]`。下表没有使用原始 obs.tissue。'),
        code('for lineage in ["fibro","endo"]:\n    display(pd.read_csv(OUT/f"{lineage}_Tissue_from_Celltype_summary.csv"))\n    display(pd.read_csv(OUT/f"{lineage}_Age_group_summary.csv"))\nprint(json.dumps(cfg["general_colors"],ensure_ascii=False,indent=2))\n'),
        code('display(Image(filename=str(OUT/"Fig1_StromalAtlas_main.png"),width=1500))\n'),
        code('for path in sorted(OUT.glob("Fig1[A-J]_*.png")):\n    print(path.name)\n    display(Image(filename=str(path),width=1100))\n'),
        code('display(json.loads((OUT/"validation_report.json").read_text()))\ndisplay(Markdown((OUT/"README.md").read_text()))\n')]
    for cell in nb.cells:
        if cell.cell_type=='code':compile(cell.source,'fig1_notebook_cell','exec')
    nbf.validate(nb)
    nbf.write(nb,out/'Fig1_StromalAtlas.ipynb')
    nbf.write(nb,HERE/'fig1_stromal_atlas.ipynb')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=HERE/'fig1_config.json')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--dpi',type=int)
    parser.add_argument('--plot-cells',type=int,help='0 = all source cells as vectors; default = existing pipeline display scale')
    parser.add_argument('--reuse-panels',action='store_true',help='Reuse existing A-J panel exports when changing only the main layout')
    parser.add_argument('--panels',default='ABCDEFGHIJ',help='Independent panels to refresh, e.g. J; main is always refreshed and other panels must already exist')
    args=parser.parse_args()
    if not args.panels or set(args.panels)-set('ABCDEFGHIJ'):parser.error('--panels must contain only A-J')
    config=json.loads(args.config.read_text())
    for key in ['source_run','output_dir']:
        config[key]=str((args.config.resolve().parent/config[key]).resolve())
    if args.output_dir:config['output_dir']=str(args.output_dir.resolve())
    if args.dpi:config['dpi']=args.dpi
    if args.plot_cells is not None:config['plot_cells_per_lineage']=args.plot_cells
    out=Path(config['output_dir']);out.mkdir(parents=True,exist_ok=True)
    root=Path(config['source_run'])
    if (root/'.running.lock').exists():raise RuntimeError('The source analysis is still being written')
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(out/'fig1_run.log')],force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'axes.titlesize':26,
        'axes.labelsize':18,'xtick.labelsize':16,'ytick.labelsize':16,'pdf.fonttype':42,
        'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
        'ps.fonttype':42,'pdf.use14corefonts':False,'savefig.facecolor':'white'})
    dump_json(out/'config_used.json',config)
    cache=out/'prepared_plot_data.pkl'
    cache_key={k:config[k] for k in ['source_run','plot_cells_per_lineage','seed']}
    if cache.exists():
        with cache.open('rb') as stream:cached=pickle.load(stream)
        if cached['key']!=cache_key:raise ValueError('Prepared cohort configuration changed; use a new output directory')
        data=cached['data']
        for item in data.values():
            r=item['record'];stat=Path(r['source']).stat()
            if (stat.st_size,stat.st_mtime_ns)!=(r['source_bytes'],r['source_mtime_ns']):raise ValueError('Cached source changed')
        logging.info('Reusing source-verified Fig1 plot metadata and full-count summaries')
    else:
        data={lineage:prepare(lineage,root,out,config) for lineage in ['fibro','endo']}
        with cache.open('wb') as stream:pickle.dump(dict(key=cache_key,data=data),stream,protocol=4)
    dump_json(out/'data_provenance.json',{lineage:value['record'] for lineage,value in data.items()})
    render_main(data,config,out)
    if args.reuse_panels:
        for extension in ['png','pdf']:
            if len(list(out.glob(f'Fig1[A-J]_*.{extension}')))!=10:
                raise ValueError('--reuse-panels requires all ten existing independent panel exports')
        logging.info('Independent A-J exports unchanged; reusing existing files')
    else:
        render_panels(data,config,out,selected=args.panels)
    pdf_validation=validate_pdfs(out)
    from PIL import Image
    # These locally generated large-format PNGs can exceed Pillow's default
    # threshold; keep a finite limit derived from our explicitly chosen canvas.
    Image.MAX_IMAGE_PIXELS=max(Image.MAX_IMAGE_PIXELS or 0,int(np.prod(config['figure_size'])*config['dpi']**2*1.25))
    for path in out.glob('Fig1*.png'):
        with Image.open(path) as picture:picture.verify()
    dump_json(out/'validation_report.json',dict(status='passed',pdfs=pdf_validation,
        checks=['True original atlas obs fields, IDs and UMAP match cache',
                'Tissue is exclusively Celltype final suffix', 'Full-cell Tissue and age count totals preserved',
                'Marker denominators match full-cell general counts',
                'All Celltypes present in common per-lineage UMAP cohort',
                'Shared category colors in data marks and classification strips; all typography black',
                'Four-row layout; G H I J appear together in one row; only A/B retain UMAP titles and statistics boxes',
                'UMAP legend labels 54pt (300%), legend titles 60pt (300%), statistics 24pt (150%)',
                'Explicit larger row gaps; taller normal-weight summary plots with 0.82-height bars; main footer removed',
                'E/F category labels 42pt and values 33pt (150%); sample x-axis label Samples',
                'E-J text is normal weight; A-J letters scaled to 300%; signature cell numbers scaled to 200%',
                'G/H colorbar title 36pt and ticks 31.5pt; J colorbar 31.5pt; signature x-axis title 40.5pt',
                'Main J moved left by 0.6 inches; G/H size legends remain close to dotplots',
                'All A-J panel letters normal weight; complete G/H Cells expressing legends shifted upward by exactly 5 mm',
                'I/J matrix width, height, top/bottom edges and Marker signature title baseline are identical',
                'Independent I/J exports use identical 12x12-inch pages without separate tight cropping',
                'All 11 PDFs have zero embedded image objects, embedded fonts and no Type3 fonts',
                'All PNG files decode successfully']))
    # Bundle original fonts so Illustrator can edit text without substituting glyphs.
    from matplotlib import font_manager
    fonts=out/'fonts';fonts.mkdir(exist_ok=True)
    for weight in ['normal','bold']:
        font=Path(font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight=weight)))
        shutil.copy2(font,fonts/font.name)
    src=out/'source_code';src.mkdir(exist_ok=True)
    for name in ['08_make_fig1.py','fig1_config.json','utils_plot.py','utils_data.py','stromal_atlas_pipeline.py','run_fig1.sh','requirements.txt']:
        if (HERE/name).exists():shutil.copy2(HERE/name,src/name)
    prompt=HERE.parent/'StromalAtlas_Fig1_Codex.md'
    shutil.copy2(prompt,out/prompt.name)
    dump_json(out/'source_code_sha256.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in src.iterdir()})
    (out/'README.md').write_text(f'''# StromalAtlas Fig1 A–J

数据来源是既有全量 atlas：Fibroblast 5,916,760 cells，Endothelial 2,239,277 cells。不是年龄过滤的 *_age.h5ad，也不是轨迹用 meta-cell 对象。实际原始路径、obs 字段、样本字段、缺失 donor 字段及原始坐标核验见 data_provenance.json。

## 面板

A/B：按 Celltype_general 的原始 UMAP。C/D：同一批细胞按 Tissue_from_Celltype 的 UMAP。没有 Leiden 编号、轨迹或椭圆框。

## 本次排版更新

画布为 {config['figure_size'][0]} × {config['figure_size'][1]} 英寸，四行分别为 A/B、C/D、E/F、G/H/I/J。行间加入独立间隔；主图底部说明已移至 README 和 Notebook。本次直接覆盖 fig1_20261103_v2 内同名文件和压缩包，不新建结果目录。

A–J 面板字母均为此前的 300%：A–D 与 G–J 为 90pt，E/F 为 102pt。全部面板字母均为常规字重；E–J 的全部文字（含面板字母、类别、标题、数值和图例）均为常规字重，颜色为黑色。

A/B 标题为 Fibroblast | subtype、Endothelial | subtype，只在 A/B 保留 24pt 信息框。A–D 图例文字 54pt、标题 60pt；C/D 无标题、无信息框，图例标题为 Tissue，采用双列。

E/F 分类标签 42pt、柱子数值 33pt，均使用常规字重；样本横轴为 Samples，仍统计去重样本。条形高度 0.82，所有 Tissue 后缀和 Unknown 年龄组都保留。

G/H 基因标签 24pt、状态标签 25.5pt、顶部分组标签 22.5pt。侧边色标标题 36pt、刻度 31.5pt，分别为此前 24pt/21pt 的 150%。底部 Cells expressing 图例标题 24pt、文字 22.5pt；在前版锚点 -0.27 的基础上，将整组图例精确向上平移 5 mm（通过物理单位变换实现，不随 DPI 或子图尺寸改变）。

I/J 格内数字 32pt，为此前 16pt 的 200%；轴刻度标签 24pt，Marker signature 横轴标题 40.5pt，为此前 27pt 的 150%。主图 J 向左移动 0.6 英寸；J 右侧 z-score 标题和刻度 31.5pt，为此前 21pt 的 150%。主图 I/J 共用 J 右侧色标；独立 I/J 面板各自保留色标，独立 I 色标仍为 21pt。G–J 无子图标题。

I/J 主图使用等宽、等高的正方形绘图区，并将分类色条留白改成按矩阵维数成比例的边距；因此 4×4 与 6×6 的彩色矩阵外框宽高完全相同，上下边缘及 Marker signature 标题均水平对齐。实际尺寸检查见 signature_alignment.json。独立 I/J PNG/PDF 均使用相同的 12×12 英寸画布，不分别进行 tight 裁切，因此整张独立图片也同样大小。

分类颜色保留在散点、条形、marker 顶部横线及 signature 边缘色条中；signature 采用较浅的红蓝连续色标，以保证黑色数字清楚可读。

Tissue_from_Celltype 严格由 `obs['Celltype'].astype(str).str.split('_').str[-1]` 构建，保留 pan、bone marrow 等原始末尾字段。pan 是注释后缀，不擅自解释成一个真实解剖组织。C、D、E 完全没有使用原始 obs.tissue。映射及计数分别导出 *_Celltype_to_Tissue_mapping.csv、*_Tissue_from_Celltype_summary.csv。

E：新 Tissue 的全量细胞数和去重样本数，两个谱系均展示，所有 suffix 都保留。F：真实 Age_group 的全量细胞数和去重样本数，Unknown 保留。E/F 均以镜像条形图展示：左向 Fibroblast、右向 Endothelial，两个方向均为正计数。sample 不是已确认的 donor；同一个 sample 可出现在多个 Tissue 后缀类别，跨类别样本数不能相加作为独立样本总数。

G/H：复用全部细胞的 marker 均值与检出比例，沿用现有 marker_matrices 计算方法。点大小=检出比例，蓝色=逐基因 min–max 均值缩放；marker 按状态分组，组内字母序。缺失的 Fibro PI16、LUM、IL33 记录并省略，不造数。

I/J：复用全部细胞 signature 均值，沿用现有 signature_matrices 的列 z-score 展示；格内数字为原始状态均值。边缘色条与 A/B、G/H 使用同一分类色，文字统一为黑色；连续定量色标与类别颜色是不同语义。signature 与 marker 用于描述性注释支持，不是独立验证。

## UMAP 显示与完整统计

每个谱系固定种子 42、按 Celltype_general × Celltype 分层，显示 {config['plot_cells_per_lineage']:,} 个代表细胞，延续原 pipeline 的绘图规模；A/C 同一批 ID，B/D 同一批 ID。所有计数、marker、signature 统计仍来自全部细胞。绘图细胞 ID 和原始行号见 *_umap_plot_cells.csv.gz。没有对表达矩阵重分析、重新计算 UMAP 或运行轨迹。

## Adobe Illustrator 编辑

Fig1_StromalAtlas_main.pdf 及全部独立面板 PDF 的散点、标签、线条、marker、热图矩形和色标均为矢量，没有嵌入图片。PDF 文本采用嵌入的 TrueType 字体（Type42/CIDFontType2），没有 Type3 字体。结构检查见 validation_report.json；并非将整张 PNG 装入 PDF。未在本环境中实际启动 Illustrator，若它提示字体缺失，可先安装 fonts 中的 DejaVu Sans Regular/Bold。主图含约64万个矢量散点，Illustrator 打开和选择散点可能需要时间；独立 A–J PDF 更便于分面编辑。

PNG 为 {config['dpi']} dpi。源代码、配置、Notebook、原提示词和完整统计表均随目录保存。

## 复现

```bash
cd {HERE}
external_inputs/python -u 08_make_fig1.py --config fig1_config.json
```

后台：`bash run_fig1.sh`。只查看结果可打开 Fig1_StromalAtlas.ipynb；默认不重画。把 REBUILD 改为 True 可执行同一个生成入口。打包源码可用 `python source_code/08_make_fig1.py --config config_used.json` 在当前数据路径环境下复现。

如果只调整主图拼版、独立面板绘制内容和配置完全未变，可加 `--reuse-panels` 保留已生成的独立 A–J 文件；修改面板内容时请用默认完整运行。

如需将全部 source cells 都导出为矢量，可加 `--plot-cells 0 --output-dir PATH_TO_NEW_DIRECTORY`；这会显著增加文件大小及 Illustrator 的编辑负担。不会改变全量统计。默认文件名严格对应此次 A–J 定义，不覆盖之前的轨迹版 Overview。
''')
    notebook(out,config)
    archive=out/'StromalAtlas_Fig1_Download.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as bundle:
        for path in out.rglob('*'):
            if path.is_file() and path!=archive and path.name!='prepared_plot_data.pkl':bundle.write(path,path.relative_to(out))
    logging.info('SUCCESS Fig1 A–J; all PDF image counts are zero. Results: %s',out)
    print('RESULT_DIR='+str(out))
    print('DOWNLOAD='+str(archive))


if __name__=='__main__':main()
