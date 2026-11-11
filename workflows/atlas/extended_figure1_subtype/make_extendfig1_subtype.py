#!/usr/bin/env python
"""Extend Fig1 A–F: verified atlas UMAP facets and full-cell composition.

Reuse the frozen Fig1 v2 palette, displayed cells and metadata helpers.
Read source H5AD annotations/UMAP only; never load the expression matrices.
"""
import argparse
import gc
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import pickle
import shutil
import sys
import zipfile

HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'tmp'/'matplotlib'))
Path(os.environ['MPLCONFIGDIR']).mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(HERE/'upstream'))

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd
from utils_data import read_elem,clean,age_order,stratified_rows

FIG1=HERE.parent/'results'/'fig1_20261103_v2'
TITLES={'fibro':'Fibroblast','endo':'Endothelial'}
STEMS={'A':'ExtendFig1_A_Fibro_subtype_umap','B':'ExtendFig1_B_Endo_subtype_umap',
       'C':'ExtendFig1_C_Fibro_tissue_proportion','D':'ExtendFig1_D_Endo_tissue_proportion',
       'E':'ExtendFig1_E_Fibro_agegroup_proportion','F':'ExtendFig1_F_Endo_agegroup_proportion'}
MAIN='ExtendFig1_subtype_umap_composition'
FONT_SCALE=1.5
BASE_FONT_SIZES={'base':20,'title':27,'panel':42,'strip':24,'celltype':13,
                 'tick':20,'axis':22,'legend':21,'note':16}
FONT_SIZES={key:value*FONT_SCALE for key,value in BASE_FONT_SIZES.items()}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,
        default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))


def stat_record(path):
    p=Path(path);s=p.stat()
    return dict(path=str(p.resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns)


def composition(obs,field,order,states):
    """Pooled-cell proportions: denominator is all lineage cells in one group."""
    counts=pd.crosstab(obs[field],obs.Celltype_general,dropna=False).reindex(
        index=order,columns=states,fill_value=0).astype('int64')
    totals=counts.sum(axis=1)
    if (totals==0).any():raise ValueError('A displayed group has an empty denominator')
    fractions=counts.div(totals,axis=0)
    if counts.to_numpy().sum()!=len(obs):raise ValueError('Composition lost source cells')
    if not np.allclose(fractions.sum(axis=1),1,atol=1e-12):raise ValueError('Proportions do not sum to one')
    summary=obs.groupby(field,observed=True).size().rename('n_cells').to_frame().reindex(order)
    for source,target in [('sample','n_samples'),('GSE','n_studies')]:
        summary[target]=obs.loc[obs[source]!='Unknown'].groupby(field,observed=True)[source].nunique().reindex(order,fill_value=0)
    summary.index.name=field;counts.index.name=field;fractions.index.name=field
    return dict(counts=counts,fractions=fractions,summary=summary)


def prepare(fig1,out,plot_cells=None,refresh=False):
    cfg=json.loads((fig1/'config_used.json').read_text())
    provenance=json.loads((fig1/'data_provenance.json').read_text())
    cache_root=Path(cfg['source_run'])
    maximum=cfg['plot_cells_per_lineage'] if plot_cells is None else plot_cells
    key=dict(version=1,fig1=str(fig1.resolve()),config_sha256=sha(fig1/'config_used.json'),
        sources={k:stat_record(v['source']) for k,v in provenance.items()},plot_cells=maximum,
        original_plot_rows={k:sha(fig1/f'{k}_umap_plot_rows.npy') for k in TITLES})
    cache=out/'prepared_subtype_plot_data.pkl'
    if cache.exists() and not refresh:
        with cache.open('rb') as f:saved=pickle.load(f)
        if saved['key']!=key:raise ValueError('Source/config changed; rerun with --refresh')
        logging.info('Using previously source-verified plot cache')
        return saved['data'],cfg,key
    data={}
    for lineage,title in TITLES.items():
        original=provenance[lineage];source=Path(original['source'])
        before=stat_record(source)
        if before['size']!=original['source_bytes'] or before['mtime_ns']!=original['source_mtime_ns']:
            raise ValueError(f'{lineage}: original atlas changed since Fig1 v2')
        logging.info('%s: reading metadata and verifying actual H5AD',title)
        full=pd.read_pickle(cache_root/lineage/'metadata.pkl')
        obs=full[['general','celltype','age_group','sample','dataset']].rename(columns={
            'general':'Celltype_general','celltype':'Celltype','age_group':'Age_group','dataset':'GSE'}).copy()
        del full
        xy=np.load(cache_root/lineage/'umap.npy',mmap_mode='r')
        with h5py.File(source,'r') as h:
            fields=list(h['obs'].keys())
            required=['Celltype_general','Celltype','Age_group','sample','GSE']
            if set(required)-set(fields):raise ValueError('Missing required actual H5AD obs columns')
            if 'X_umap' not in h['obsm'] or h['obsm/X_umap'].shape!=(len(obs),2):
                raise ValueError('No aligned source UMAP')
            for field in required:
                actual=clean(pd.Series(read_elem(h['obs'][field])))
                if not np.array_equal(actual.to_numpy(),obs[field].to_numpy()):
                    raise ValueError(f'{lineage}: cached/source annotation mismatch in {field}')
                del actual
            idx=h['obs'][h['obs'].attrs['_index']]
            ids=obs.index.to_numpy()
            for start in range(0,len(obs),250000):
                stop=min(start+250000,len(obs))
                if not isinstance(idx,h5py.Dataset):raise ValueError('Expected an HDF5 cell-ID dataset')
                if not np.array_equal(idx.asstr()[start:stop],ids[start:stop]):
                    raise ValueError('Source and cached observation IDs/order differ')
                source_xy=h['obsm/X_umap'][start:stop].astype(np.float32)
                if not np.isfinite(source_xy).all() or not np.array_equal(source_xy,xy[start:stop]):
                    raise ValueError('Source UMAP is not the Fig1 atlas UMAP')
            donor_fields=[name for name in ['donor','donor_id','patient_id','individual'] if name in fields]
            # Report observed numeric ages without inventing boundaries for A–H.
            age_ranges=None
            if 'Age' in h['obs'] and isinstance(h['obs/Age'],h5py.Dataset):
                ages=np.asarray(h['obs/Age'][:],dtype=float)
                age_ranges=pd.DataFrame({'Age_group':obs.Age_group.to_numpy(),'Age':ages})
                age_ranges=age_ranges.loc[np.isfinite(ages)&(ages>=0)].groupby('Age_group').Age.agg(['count','min','median','max'])
                age_ranges.to_csv(out/'tables'/f'{lineage}_observed_numeric_age_by_group.csv')
        states=cfg['general_order'][lineage]
        if set(obs.Celltype_general)!=set(states):raise ValueError('Unexpected source subtype categories')
        if set(states)-set(cfg['general_colors']):raise ValueError('Fig1 v2 palette is incomplete')
        # The strict requested Tissue definition; original obs.tissue is unused.
        obs['Tissue_from_Celltype']=obs.Celltype.astype(str).str.split('_').str[-1].astype('category')
        tissue_order=sorted(obs.Tissue_from_Celltype.unique(),key=lambda s:(s!='pan',s))
        ages_order=age_order(obs.Age_group)
        tissue=composition(obs,'Tissue_from_Celltype',tissue_order,states)
        age=composition(obs,'Age_group',ages_order,states)
        for kind,result in [('tissue',tissue),('agegroup',age)]:
            for metric,table in result.items():table.to_csv(out/'tables'/f'{lineage}_{kind}_{metric}.csv')
            long=result['counts'].stack().rename('n_cells').reset_index()
            long['fraction']=result['fractions'].stack().to_numpy()
            long['percent']=100*long.fraction
            long.to_csv(out/'tables'/f'{lineage}_{kind}_composition_long.csv',index=False)
        if plot_cells is None:
            rows=np.load(fig1/f'{lineage}_umap_plot_rows.npy')
            reuse=True
        else:
            rows=stratified_rows(obs,len(obs) if maximum==0 else min(maximum,len(obs)),cfg['seed'],
                                 ('Celltype_general','Celltype'))
            reuse=False
        if len(np.unique(rows))!=len(rows) or rows.min()<0 or rows.max()>=len(obs):raise ValueError('Invalid plotting rows')
        plotted=obs.iloc[rows].copy()
        if set(plotted.Celltype)!=set(obs.Celltype):raise ValueError('Plot cohort lost a fine Celltype')
        plotted.insert(0,'source_row',rows)
        plotted['UMAP1']=xy[rows,0];plotted['UMAP2']=xy[rows,1]
        plotted.to_csv(out/'tables'/f'{lineage}_umap_plot_cells.csv.gz',index_label='cell_id')
        np.save(out/'tables'/f'{lineage}_umap_plot_rows.npy',rows)
        fine_counts=obs.groupby(['Celltype_general','Celltype','Tissue_from_Celltype'],observed=True).size().rename('n_cells').reset_index()
        fine_counts.to_csv(out/'tables'/f'{lineage}_Celltype_to_Tissue_mapping.csv',index=False)
        fine_medians=plotted.groupby('Celltype',observed=True)[['UMAP1','UMAP2']].median()
        anchors=[]
        for _,record in fine_counts.iterrows():
            mask=plotted.Celltype.eq(record.Celltype).to_numpy()
            coords=np.asarray(xy[rows[mask]])
            med=fine_medians.loc[record.Celltype].to_numpy(dtype=float)
            point=coords[np.argmin(((coords-med)**2).sum(axis=1))]
            anchors.append(dict(Celltype_general=record.Celltype_general,Celltype=record.Celltype,
                Tissue_from_Celltype=record.Tissue_from_Celltype,n_cells=int(record.n_cells),
                x=float(point[0]),y=float(point[1]),label=str(record.Celltype).split('_',1)[1].replace('_',' · ')))
        anchors=pd.DataFrame(anchors)
        anchors.to_csv(out/'tables'/f'{lineage}_Celltype_label_anchors.csv',index=False)
        counts=obs.Celltype_general.value_counts().reindex(states).astype(int)
        counts.rename('n_cells').to_csv(out/'tables'/f'{lineage}_subtype_cell_counts.csv')
        known=obs.Age_group.ne('Unknown')
        note=dict(source=before,fields_checked=required,actual_obs_columns=fields,
            n_cells=len(obs),n_plot_cells=len(rows),n_samples=int(obs.loc[obs['sample']!='Unknown','sample'].nunique()),
            donor_fields=donor_fields,sample_is_verified_donor=False,source_annotations_ids_umap_match=True,
            source_unchanged=stat_record(source)==before,reused_exact_fig1_plot_rows=reuse,
            n_Celltypes=int(obs.Celltype.nunique()),subtypes=states,tissue_order=list(tissue_order),age_order=ages_order,
            n_known_age_cells=int(known.sum()),n_unknown_age_cells=int((~known).sum()),
            unknown_age_percent=100*(~known).mean(),original_obs_tissue_used=False,
            tissue_definition="Celltype.astype(str).str.split('_').str[-1]",
            proportions='Pooled source-cell counts within Tissue/age group, summing to 1; not donor/sample means',
            minimum_tissue_cells=int(tissue['summary'].n_cells.min()),source_umap_recomputed=False)
        if not note['source_unchanged']:raise ValueError('Source changed during preparation')
        data[lineage]=dict(lineage=lineage,obs=plotted,xy=np.array(xy[rows]),
            bounds=np.array([xy.min(axis=0),xy.max(axis=0)]),counts=counts,anchors=anchors,
            tissue=tissue,agegroup=age,record=note)
        logging.info('%s verified: %s full cells, %s UMAP points; %s Tissue suffixes; ages %s',
                     title,len(obs),len(rows),len(tissue_order),ages_order)
        del obs,xy,plotted,ids
        gc.collect()
    dump(out/'data_provenance.json',{k:v['record'] for k,v in data.items()})
    dump(out/'color_map.json',dict(source=str(fig1/'config_used.json'),source_sha256=key['config_sha256'],
                                  general_colors=cfg['general_colors'],general_order=cfg['general_order']))
    with cache.open('wb') as f:pickle.dump(dict(key=key,data=data),f,protocol=4)
    return data,cfg,key


def load_celltype_colors(data,cfg,out):
    """Use the actual category-aligned palettes saved in each original H5AD."""
    palettes={};sources={}
    for lineage,item in data.items():
        path=item['record']['source']['path']
        with h5py.File(path,'r') as h:
            labels=list(map(str,read_elem(h['obs/Celltype/categories'])))
            values=list(map(str,read_elem(h['uns/Celltype_colors'])))
        if len(labels)!=len(values):raise ValueError('Source Celltype/color length mismatch')
        palette=dict(zip(labels,values))
        if set(item['obs'].Celltype)!=set(palette):raise ValueError('Missing actual Celltype colors')
        if not all(matplotlib.colors.is_color_like(v) for v in palette.values()):raise ValueError('Invalid source color')
        if len(set(matplotlib.colors.to_hex(v) for v in palette.values()))!=len(palette):
            raise ValueError('Actual Celltypes do not have distinct colors')
        palettes[lineage]=palette
        sources[lineage]=dict(h5ad=path,labels='obs/Celltype/categories',colors='uns/Celltype_colors')
    cfg['celltype_colors']=palettes
    dump(out/'color_map.json',dict(source=str(FIG1/'config_used.json'),
        general_colors=cfg['general_colors'],general_order=cfg['general_order'],
        celltype_colors=palettes,celltype_palette_sources=sources))
    dump(out/'typography_config.json',dict(scale_from_previous=FONT_SCALE,previous=BASE_FONT_SIZES,
        current=FONT_SIZES,umap_axes=13*FONT_SCALE,figure_heading=29*FONT_SCALE,
        color='#000000',weight='normal'))


def subtype_panel(fig,spec,data,cfg,letter):
    lineage=data['lineage'];states=cfg['general_order'][lineage];colors=cfg['general_colors']
    fine_colors=cfg['celltype_colors'][lineage]
    ncols=2 if lineage=='fibro' else 3
    grid=spec.subgridspec(2,1,height_ratios=[.07,1],hspace=.04)
    header=fig.add_subplot(grid[0]);header.set_axis_off()
    header.text(0,.55,letter,transform=header.transAxes,fontsize=FONT_SIZES['panel'],va='center')
    header.text(.12,.55,f'{data["record"]["n_cells"]:,} atlas cells · {len(data["obs"]):,} displayed',
                transform=header.transAxes,fontsize=FONT_SIZES['note'],va='center')
    cells=grid[1].subgridspec(2,ncols,hspace=.16,wspace=.10)
    xy=data['xy'];values=data['obs'].Celltype_general.to_numpy();axes=[]
    fine_values=data['obs'].Celltype.to_numpy()
    for pos,state in enumerate(states):
        # The former list-legend region now belongs to the enlarged UMAP.
        facet=cells[pos//ncols,pos%ncols].subgridspec(2,1,height_ratios=[.07,1],hspace=.06)
        strip=fig.add_subplot(facet[0]);strip.set_axis_off()
        strip.add_patch(Rectangle((0,.05),1,.9,transform=strip.transAxes,
                                 facecolor=colors[state],edgecolor='#43484C',linewidth=.7))
        strip.text(.5,.5,state,ha='center',va='center',transform=strip.transAxes,
                   fontsize=FONT_SIZES['strip'],color='black',fontweight='normal')
        ax=fig.add_subplot(facet[1]);axes.append(ax)
        order=np.random.default_rng(42).permutation(len(xy))
        selected=order[values[order]==state]
        points=xy[selected]
        ax.scatter(points[:,0],points[:,1],s=1.8,
                   c=[fine_colors[name] for name in fine_values[selected]],alpha=.95,
                   edgecolors='none',rasterized=False,zorder=2)
        # Equal physical coordinate scaling, with an independently zoomed
        # viewport for this subtype; original cell coordinates never move.
        low,high=points.min(axis=0),points.max(axis=0)
        span=np.maximum(high-low,1e-6);center=(low+high)/2
        position=ax.get_position()
        inches=np.array([position.width*fig.get_figwidth(),position.height*fig.get_figheight()])
        scale=np.max(span/inches)*1.32
        limits=np.array([center-scale*inches/2,center+scale*inches/2])
        ax.set_xlim(limits[:,0]);ax.set_ylim(limits[:,1]);ax.set_aspect('equal',adjustable='box')
        ax.set_axis_off()
        records=data['anchors'].loc[data['anchors'].Celltype_general==state]
        count_text=ax.text(.98,.015,f'n = {data["counts"][state]:,}\n{len(records)} Celltype'+('s' if len(records)!=1 else ''),
                transform=ax.transAxes,ha='right',va='bottom',fontsize=FONT_SIZES['note'],
                bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2),zorder=6)
        ax.annotate('',xy=(.255,.035),xytext=(.035,.035),xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='->',color='black',lw=1.0),zorder=7)
        ax.annotate('',xy=(.035,.255),xytext=(.035,.035),xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='->',color='black',lw=1.0),zorder=7)
        ax.text(.145,.003,'UMAP-1',ha='center',va='top',transform=ax.transAxes,fontsize=13*FONT_SCALE)
        ax.text(.002,.145,'UMAP-2',ha='right',va='center',rotation=90,transform=ax.transAxes,fontsize=13*FONT_SCALE)
        names=[name for name in fine_colors if name in set(records.Celltype)]
        records=records.set_index('Celltype').loc[names]
        texts=[ax.text(row.x,row.y,name,fontsize=FONT_SIZES['celltype'],
                color='black',ha='center',va='center',zorder=6,
                bbox=dict(facecolor='white',edgecolor='none',alpha=.86,pad=1))
                for name,row in records.iterrows()]
        ax._callout_texts=texts
        ax._callout_targets=records[['x','y']].to_numpy()
        ax._callout_names=names
        ax._callout_points=points
        ax._callout_obstacles=[count_text,Rectangle((0,0),.28,.28,transform=ax.transAxes,visible=False)]
        ax._target_subtype=state
        ax._drawn_subtypes=sorted(set(values[selected]))
        ax._foreground_celltypes=sorted(set(fine_values[selected]))
        ax._point_count=len(points)
        ax._zoom_bounds=limits.tolist()
        if set(ax._foreground_celltypes)!=set(names):raise ValueError('Callouts and colored cells differ')
        for name,target in zip(names,ax._callout_targets):
            matching=points[fine_values[selected]==name]
            if not np.any(np.all(matching==target,axis=1)):
                raise ValueError('Callout anchor is not an actual displayed Celltype cell')
    return axes


def count_label(value):
    return f'{int(value):,}'


def proportion_panel(fig,spec,data,cfg,kind,letter):
    lineage=data['lineage'];states=cfg['general_order'][lineage];colors=cfg['general_colors']
    grid=spec.subgridspec(2,1,height_ratios=[.15,1],hspace=.07)
    legend_ax=fig.add_subplot(grid[0]);legend_ax.set_axis_off()
    legend_ax.text(0,.5,letter,transform=legend_ax.transAxes,fontsize=FONT_SIZES['panel'],va='center')
    handles=[Line2D([],[],marker='o',ls='',ms=15,color=colors[s],label=s) for s in states]
    legend=legend_ax.legend(handles=handles,loc='center left',bbox_to_anchor=(.09,.5),
                     ncol=len(states),frameon=False,fontsize=FONT_SIZES['legend'],
                     handlelength=.65,handletextpad=.35,columnspacing=.9,labelspacing=.2,
                     borderaxespad=0,borderpad=0)
    legend_ax._checked_legend=legend;legend_ax._legend_kind='subtype'
    legend_ax._full_celltypes_expected=states
    # Explicit gutters keep long group names and counts clear of plot edges.
    body=grid[1].subgridspec(1,3,width_ratios=[.22,.60,.18],wspace=.015)
    ax=fig.add_subplot(body[1]);lab=fig.add_subplot(body[0]);count=fig.add_subplot(body[2])
    lab.set_axis_off();count.set_axis_off()
    result=data[kind];table=result['fractions'];summary=result['summary'];labels=list(table.index)
    y=np.arange(len(labels));left=np.zeros(len(labels))
    for state in states:
        width=100*table[state].to_numpy()
        ax.barh(y,width,left=left,height=.77,color=colors[state],edgecolor='none',rasterized=False)
        left+=width
    ax.set_xlim(0,100);ax.set_ylim(len(labels)-.5,-.5);ax.set_yticks([])
    ax.set_xticks([0,25,50,75,100]);ax.tick_params(axis='x',labelsize=FONT_SIZES['tick'],length=4,pad=7)
    ax.set_xlabel(f'Percentage of {TITLES[lineage]} cells (%)',fontsize=FONT_SIZES['axis'],labelpad=12)
    for side in ['left','top','right']:ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_linewidth(.8)
    ax.set_axisbelow(True);ax.grid(axis='x',color='#E9EBED',lw=.8)
    lab.set_ylim(ax.get_ylim());count.set_ylim(ax.get_ylim())
    for row,label in enumerate(labels):
        display='Unknown' if label=='Unknown' else str(label)
        lab.text(.97,row,display,ha='right',va='center',fontsize=FONT_SIZES['tick'],transform=lab.transData)
        count.text(.02,row,count_label(summary.loc[label,'n_cells']),ha='left',va='center',
                   fontsize=FONT_SIZES['note'],transform=count.transData)
    count.text(.02,1.035,'Cells',transform=count.transAxes,fontsize=FONT_SIZES['note'],va='bottom')
    ax._composition_kind=kind;ax._composition_lineage=lineage
    return ax


def place_umap_callouts(ax,renderer):
    """Pack full names near real cells, avoiding dense clouds and other text."""
    texts=ax._callout_texts;targets=ax._callout_targets
    bounds=ax.get_window_extent(renderer)
    occupied=[obj.get_window_extent(renderer).expanded(1.08,1.15).extents
              for obj in ax._callout_obstacles]
    boxes=[text.get_window_extent(renderer) for text in texts]
    # Integral histogram allows rapid density scoring of candidate text boxes.
    cloud=ax.transData.transform(ax._callout_points)
    bins=96
    hist,_,_=np.histogram2d(cloud[:,0],cloud[:,1],bins=bins,
        range=[[bounds.x0,bounds.x1],[bounds.y0,bounds.y1]])
    integral=np.pad(hist,((1,0),(1,0))).cumsum(axis=0).cumsum(axis=1)
    order=sorted(range(len(texts)),key=lambda i:boxes[i].width,reverse=True)
    connectors=[]
    for i in order:
        width=boxes[i].width+9;height=boxes[i].height+9
        if width>=bounds.width-4 or height>=bounds.height-4:
            raise ValueError('Full Celltype label does not fit the enlarged UMAP')
        xx=np.linspace(bounds.x0+width/2+2,bounds.x1-width/2-2,110)
        yy=np.linspace(bounds.y0+height/2+2,bounds.y1-height/2-2,110)
        gx,gy=np.meshgrid(xx,yy);centers=np.c_[gx.ravel(),gy.ravel()]
        candidate=np.c_[centers[:,0]-width/2,centers[:,1]-height/2,
                        centers[:,0]+width/2,centers[:,1]+height/2]
        valid=np.ones(len(candidate),dtype=bool)
        for x0,y0,x1,y1 in occupied:
            valid&=~((candidate[:,0]<x1)&(candidate[:,2]>x0)&
                     (candidate[:,1]<y1)&(candidate[:,3]>y0))
        target=ax.transData.transform(targets[i])
        # A visible leader segment must remain outside its own label box.
        valid&=~((candidate[:,0]-10<target[0])&(candidate[:,2]+10>target[0])&
                 (candidate[:,1]-10<target[1])&(candidate[:,3]+10>target[1]))
        if not valid.any():raise ValueError('No non-overlapping full-name callout position')
        indices=np.c_[np.floor((candidate[:,0]-bounds.x0)/bounds.width*bins),
                      np.floor((candidate[:,1]-bounds.y0)/bounds.height*bins),
                      np.ceil((candidate[:,2]-bounds.x0)/bounds.width*bins),
                      np.ceil((candidate[:,3]-bounds.y0)/bounds.height*bins)]
        indices=np.clip(indices.astype(int),0,bins)
        x0,y0,x1,y1=indices.T
        coverage=(integral[x1,y1]-integral[x0,y1]-integral[x1,y0]+integral[x0,y0])/len(cloud)
        distance=(((centers-target)/[bounds.width,bounds.height])**2).sum(axis=1)
        cost=distance+2.5*coverage;cost[~valid]=np.inf
        best=int(np.argmin(cost));occupied.append(candidate[best])
        position=ax.transData.inverted().transform(centers[best])
        texts[i].set_position(position)
        connector=ax.annotate('',xy=targets[i],xytext=position,
            arrowprops=dict(arrowstyle='-',color='#303030',lw=.8,shrinkA=0,shrinkB=0),zorder=4)
        connectors.append(connector)
    ax._callout_connectors=connectors


def finalize_figure(fig):
    # Apply the requested typography to every rendered text artist.
    fig.canvas.draw()
    for artist in fig.findobj(match=Text):
        artist.set_color('#000000')
        artist.set_fontweight('normal')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for ax in fig.axes:
        if hasattr(ax,'_callout_texts'):place_umap_callouts(ax,renderer)
    fig.canvas.draw()


def check_layout(fig,umap_axes,label):
    renderer=fig.canvas.get_renderer();outside=[];texts=list(fig.texts);legend_items=0
    callout_items=0;leader_lines=0;point_counts={};zoom_limits={};plot_sizes={}
    legend_overflows=[];single_rows=[]
    for ax in fig.axes:
        texts.extend(ax.texts)
        if ax.axison:
            texts.extend([ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()])
        legend=ax.get_legend()
        if legend:
            texts.extend(legend.get_texts())
            if legend.get_title().get_text():raise ValueError('Unexpected legend title')
        if hasattr(ax,'_checked_legend'):
            legend=ax._checked_legend
            names=[t.get_text() for t in legend.get_texts()]
            if names!=ax._full_celltypes_expected:raise ValueError('Incomplete or abbreviated legend')
            box=legend.get_window_extent(renderer);area=ax.get_window_extent(renderer)
            if box.x0<area.x0-2 or box.y0<area.y0-2 or box.x1>area.x1+2 or box.y1>area.y1+2:
                legend_overflows.append(names)
            if ax._legend_kind=='Celltype':legend_items+=len(names)
            else:
                y=[t.get_transform().transform(t.get_position())[1] for t in legend.get_texts()]
                single_rows.append(float(np.ptp(y))<1)
        if hasattr(ax,'_callout_texts'):
            names=[text.get_text() for text in ax._callout_texts]
            if names!=ax._callout_names:raise ValueError('Incomplete or abbreviated Celltype callouts')
            boxes=[text.get_window_extent(renderer).expanded(1.01,1.03) for text in ax._callout_texts]
            reserved=[obj.get_window_extent(renderer).expanded(1.03,1.05) for obj in ax._callout_obstacles]
            if any(box.overlaps(other) for i,box in enumerate(boxes) for other in [*boxes[:i],*reserved]):
                raise ValueError('Celltype callouts overlap another label or the UMAP annotations')
            if ax._drawn_subtypes!=[ax._target_subtype] or len(ax.collections)!=1:
                raise ValueError('UMAP contains a non-target subtype/background cloud')
            if len(ax.collections[0].get_offsets())!=ax._point_count:raise ValueError('UMAP point count mismatch')
            if len(ax._callout_connectors)!=len(names):raise ValueError('Missing Celltype leader lines')
            callout_items+=len(names);leader_lines+=len(ax._callout_connectors)
            point_counts[ax._target_subtype]=ax._point_count
            zoom_limits[ax._target_subtype]=[*ax.get_xlim(),*ax.get_ylim()]
            area=ax.get_window_extent(renderer)
            plot_sizes[ax._target_subtype]=[area.width/fig.dpi,area.height/fig.dpi]
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -2 or b.y0 < -2 or b.x1>fig.bbox.width+2 or b.y1>fig.bbox.height+2:outside.append(t.get_text())
        if matplotlib.colors.to_hex(t.get_color())!='#000000' or t.get_fontweight() not in ['normal',400]:
            raise ValueError('Text must be pure black and normal weight')
        if (' | ' in t.get_text() or 'Celltype suffix groups' in t.get_text()
                or 'Biological order' in t.get_text()):raise ValueError('Removed panel title/note remains')
    if outside:raise ValueError(f'{label}: text outside canvas: {outside}')
    if legend_overflows:raise ValueError(f'{label}: legends exceed their reserved space: {legend_overflows}')
    if not all(single_rows):raise ValueError('Subtype legend is not a single row')
    facet_limits={}
    for ax in umap_axes:
        # Diagnostic only: subtype-only facets intentionally use independent
        # zoomed viewports while preserving the original cell coordinates.
        key=tuple(np.round([*ax.get_xlim(),*ax.get_ylim()],8))
        facet_limits[key]=facet_limits.get(key,0)+1
    return dict(outside_canvas_texts=outside,celltype_label_overlaps=0,legend_overflows=[],
                full_celltype_legend_items=legend_items,subtype_legends_single_row=True,
                full_celltype_callout_items=callout_items,leader_lines=leader_lines,
                subtype_only=True,grey_background_points=0,umap_zoom='per_subtype',
                plotted_cells_per_subtype=point_counts,zoom_limits=zoom_limits,umap_size_inches=plot_sizes,
                black_normal_text=True,font_scale=FONT_SCALE,panel_titles_removed=True,
                identical_limits_per_lineage=list(facet_limits.values()),
                figsize_inches=fig.get_size_inches().tolist())


def export(fig,stem,dpi):
    for artist in fig.findobj():
        if hasattr(artist,'set_rasterized'):artist.set_rasterized(False)
    logging.info('Exporting %s',stem.name)
    inches=fig.get_size_inches().copy()
    fig.set_size_inches(np.nextafter(np.rint(inches*dpi)/dpi,np.inf),forward=False)
    try:fig.savefig(stem.with_suffix('.png'),dpi=dpi,bbox_inches=None,facecolor='white')
    finally:fig.set_size_inches(inches,forward=False)
    fig.savefig(stem.with_suffix('.pdf'),dpi=dpi,bbox_inches=None,facecolor='white')
    plt.close(fig)


def render(data,cfg,out,dpi,preview=False):
    reports={}
    fig=plt.figure(figsize=(40,42))
    fig.subplots_adjust(left=.035,right=.98,bottom=.04,top=.955)
    gs=fig.add_gridspec(3,2,height_ratios=[2.20,.95,.80],width_ratios=[1,1.20],hspace=.11,wspace=.10)
    fig.text(.035,.99,'Extend Fig1. Subtype UMAP and composition across tissue and age groups',
             fontsize=29*FONT_SCALE,ha='left',va='top')
    facets=[]
    for col,lineage in enumerate(TITLES):
        facets+=subtype_panel(fig,gs[0,col],data[lineage],cfg,'AB'[col])
        proportion_panel(fig,gs[1,col],data[lineage],cfg,'tissue','CD'[col])
        proportion_panel(fig,gs[2,col],data[lineage],cfg,'agegroup','EF'[col])
    finalize_figure(fig)
    reports['main']=check_layout(fig,facets,'main')
    if preview:
        fig.savefig(out/f'{MAIN}_preview.png',dpi=60,bbox_inches=None)
        plt.close(fig)
        dump(out/'preview_layout_validation.json',reports)
        return
    export(fig,out/MAIN,dpi)
    for letter in 'ABCDEF':
        lineage='fibro' if letter in 'ACE' else 'endo'
        if letter in 'AB':
            fig=plt.figure(figsize=(24 if letter=='A' else 30,25))
            fig.subplots_adjust(left=.045,right=.975,bottom=.05,top=.965)
            spec=fig.add_gridspec(1,1)[0]
            axes=subtype_panel(fig,spec,data[lineage],cfg,letter)
            finalize_figure(fig)
        else:
            fig=plt.figure(figsize=(19,12 if letter in 'CD' else 10))
            fig.subplots_adjust(left=.03,right=.97,bottom=.12,top=.97)
            spec=fig.add_gridspec(1,1)[0]
            proportion_panel(fig,spec,data[lineage],cfg,'tissue' if letter in 'CD' else 'agegroup',letter)
            axes=[];finalize_figure(fig)
        reports[letter]=check_layout(fig,axes,letter)
        export(fig,out/STEMS[letter],dpi)
    dump(out/'layout_validation.json',reports)


def validate(out,data,cfg,dpi):
    from PIL import Image
    from pypdf import PdfReader
    # All images here are trusted, locally generated scientific figures.
    Image.MAX_IMAGE_PIXELS=250_000_000
    reports={}
    for name in [MAIN,*STEMS.values()]:
        pdf=PdfReader(out/f'{name}.pdf');page=pdf.pages[0];fonts=[];images=[]
        def inspect(resources):
            resources=resources.get_object()
            if '/Font' in resources:
                for ref in resources['/Font'].get_object().values():
                    font=ref.get_object()
                    if font.get('/Subtype')=='/Type3':raise ValueError('Type3 PDF font')
                    for child in font.get('/DescendantFonts',[font]):
                        desc=child.get_object().get('/FontDescriptor')
                        if desc:fonts.append(any(k in desc.get_object() for k in ['/FontFile','/FontFile2','/FontFile3']))
            if '/XObject' in resources:
                for ref in resources['/XObject'].get_object().values():
                    obj=ref.get_object()
                    if obj.get('/Subtype')=='/Image':images.append(True)
                    if '/Resources' in obj:inspect(obj['/Resources'])
        inspect(page['/Resources'])
        if len(pdf.pages)!=1 or images or not fonts or not all(fonts):raise ValueError('PDF vector/font validation failed')
        with Image.open(out/f'{name}.png') as im:
            pixels=im.size;png_dpi=im.info.get('dpi');im.verify()
        expected=tuple(round(float(v)/72*dpi) for v in [page.mediabox.width,page.mediabox.height])
        if pixels!=expected or not png_dpi or not np.allclose(png_dpi,[dpi,dpi],atol=.02):raise ValueError('PNG size/dpi mismatch')
        reports[name]=dict(png_pixels=pixels,dpi=png_dpi,pdf_images=0,embedded_fonts=True)
    for lineage,item in data.items():
        for kind in ['tissue','agegroup']:
            if not np.allclose(item[kind]['fractions'].sum(axis=1),1):raise ValueError('Fractions fail sum-to-one')
            if item[kind]['counts'].to_numpy().sum()!=item['record']['n_cells']:raise ValueError('Counts fail source total')
        if cfg['general_colors']!=json.loads((out/'color_map.json').read_text())['general_colors']:
            raise ValueError('Palette differs from original Fig1 v2')
        if cfg['celltype_colors'][lineage]!=json.loads((out/'color_map.json').read_text())['celltype_colors'][lineage]:
            raise ValueError('Fine Celltype palette differs from the actual source mapping')
    dump(out/'validation_report.json',dict(status='passed',figures=reports,
        full_count_totals={k:v['record']['n_cells'] for k,v in data.items()},
        proportions_sum_to_one=True,source_fields_ids_umap_verified=True,
        original_obs_tissue_used=False,exact_fig1_palette=True,actual_h5ad_Celltype_colors=True,
        font_scale=FONT_SCALE,text_color='#000000',font_weight='normal',
        umap_subtype_only=True,umap_grey_background=False,umap_zoom='per_subtype',
        Celltype_labels='full names with leader lines, on the UMAP'))


def documents(out,data,cfg,key,dpi):
    import nbformat as nbf
    from matplotlib import font_manager
    env={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','anndata','h5py','matplotlib','Pillow','pypdf','nbformat']}
    if out!=HERE:
        (out/'upstream').mkdir(exist_ok=True)
        for name in ['make_extendfig1_subtype.py','run_extendfig1_subtype.sh','input_prompt.md',
                     'upstream/utils_data.py','upstream/08_make_fig1_reference.py']:
            shutil.copy2(HERE/name,out/name)
    dump(out/'environment.json',dict(python=sys.version,packages=env))
    (out/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in env.items())+'\n')
    rare=data['endo']['tissue']['summary'].loc['lymph']
    (out/'README.md').write_text(f'''# Extend Fig1. Subtype UMAP and composition across tissue and age groups

实际读取 Fibroblast {data['fibro']['record']['n_cells']:,} 与 Endothelial {data['endo']['record']['n_cells']:,} 个源细胞，生成 A–F 主图和六个独立 panel。PNG 为 {dpi} dpi；PDF 中点、线、条形、文字均为矢量，TrueType 字体嵌入，无栅格图片。验证见 validation_report.json；未在本环境实际启动 Adobe Illustrator。

## 真实输入与复用

- Fibroblast：`{data['fibro']['record']['source']['path']}`。
- Endothelial：`{data['endo']['record']['source']['path']}`。
- 颜色、类别顺序：`{key['fig1']}/config_used.json` 的 general_colors / general_order；原样复制到 color_map.json，未自行调色。
- 复用冻结 Fig1 v2 的 utils_data.py（read_elem、clean、age_order、stratified_rows）及其 source_run metadata.pkl/umap.npy。逐项核验实际 H5AD 的 Celltype_general、Celltype、Age_group、sample、GSE、全部细胞 ID/顺序及全部 UMAP 坐标，并核验源文件大小和 mtime。无需读取数百 GB 的表达矩阵，不运行轨迹、邻居图或 UMAP 拟合。

## A–B：Subtype UMAP

Fibroblast 按 Myfibro/eFibro/iFibro/qFibro 排列，2×2；Endothelial 按 ArtEC/CapEC/VenEC/EndoMT/MCs/LEC 排列，2×3。**每个小图只绘制对应 subtype 的展示细胞**，不绘制其他 subtype，不保留灰色背景点，也不再使用图下方的列表图例。原图例区域全部交给 UMAP，使图形占用面积明显增大。

每个 subtype 内按完整 Celltype 字段分别着色。颜色从每个实际 H5AD 的 uns['Celltype_colors'] 与 obs['Celltype'].cat.categories 按原顺序配对，Fibroblast 32 类、Endothelial 38 类；每个谱系内颜色均不同。完整 Celltype 名称直接放在 UMAP 内，通过引线连接到该 Celltype 的一个实际展示细胞；锚点沿用已核验的近中位位置细胞，不移动原始细胞坐标。标签位置避让其他文字、计数和坐标箭头，并尽量避开稠密细胞区域。每个名称均保留原始下划线字符串，不缩写、不改名。70 个标注均有对应引线。

各 subtype 根据自身全部展示细胞的坐标范围独立缩放，并保留 1:1 的 UMAP 坐标物理比例；没有裁掉离群点，没有重新拟合或非线性变形。不同小图的坐标范围因此不再相同，不能依据画面中的距离直接比较不同小图的绝对尺度。完整坐标范围、图形物理尺寸及每个 subtype 展示点数记录在 layout_validation.json。

默认仍逐行复用 Fig1 v2 的每谱系 160,000 个分层代表细胞；各小图分别取其中属于当前 subtype 的部分，合并后仍恰好 160,000 个/谱系，保留所有细分 Celltype。这是绘图抽样，不是下采样统计。`--plot-cells 0` 可画全部源细胞。UMAP 图中的 n 为全量源 subtype 细胞数。subtype 色条和 C–F 的 subtype 堆叠条、图例仍保持 Fig1 v2 的 general_colors；细分 Celltype 与总体 subtype 的两个颜色层级均写入 color_map.json。

本次按要求在原路径覆盖主图、独立 panel、代码与文档。C–F 沿用上次版式和统计，仍无 panel 标题、图例标题及额外说明，subtype 图例为一行。保留 A–F 字母与 subtype 色条名称。字号维持已确认的初版 150%，本次不额外增加字号：UMAP 全名标注 19.5pt，subtype 图例 31.5pt，坐标刻度 30pt；所有文字均为 #000000、normal 字重。绝对字号见 typography_config.json。检查覆盖：仅目标 subtype、无背景云、完整名称和引线、标注无重叠、单行 subtype 图例、画布无溢出、黑色及正常字重。

## C–F：比例的定义与顺序

C/D 的 Tissue **只取** `obs['Celltype'].astype(str).str.split('_').str[-1]`，显式建立 Tissue_from_Celltype。没有使用 obs.tissue。`pan` 是注释后缀类别，不能解读为具体解剖组织。

这些是 Celltype 命名定义的后缀组，不是按真实采样组织重新分组。不少后缀只对应一个 subtype，因此 C/D 出现 100% 单色条是该注释定义的直接结果，不能据此声称真实器官只存在这一种 subtype。

每条为一个 Tissue/年龄组的 100% 横向堆叠条：分子为该组该 subtype 的全部源细胞数；分母为同一谱系同一组的全部源细胞数。每条比例加总为 1。采用 pooled cell composition，不是个体条形图，不是供体或样本均值，也没有统计显著性检验。

Tissue 顺序：pan 首位，其余后缀按字母顺序；保留全部类别，不合并 Others。Age_group 使用已有项目 age_order，按生物学 A→H，Unknown 单独放最后；不把未知年龄并入任何已知组。A–H 的观察到的数值年龄范围单独导出，不凭组名字母推断未经核实的年龄边界。

Fibroblast Unknown 年龄：{data['fibro']['record']['n_unknown_age_cells']:,} 个（{data['fibro']['record']['unknown_age_percent']:.2f}%）；Endothelial Unknown：{data['endo']['record']['n_unknown_age_cells']:,} 个（{data['endo']['record']['unknown_age_percent']:.2f}%）。年龄组比例需结合该缺失程度解释。Endothelial 的 lymph 后缀组只有 {int(rare.n_cells)} 个细胞、{int(rare.n_samples)} 个样本，保留并通过右侧 Cells 数量明示；小组不宜过度解释。

原始数值字段 Age 也存在需要核对的条目：例如 A 组的最大数值在 Fibroblast 为 442、Endothelial 为 200，数值单位和质量尚未核实。本图直接沿用原始 Age_group，不根据 Age 自动重分组，也不将这些数值直接解读为年岁。tables/*_observed_numeric_age_by_group.csv 仅用于核查原始数值范围。

sample 字段用于核对去重样本数，但未发现已验证的 donor/patient/individual ID，不能把 sample 当作独立供体。各组样本数见 tables/*_summary.csv；同一 sample 可跨组，跨组样本数不能相加得到独立供体总数。

## 输出与复现

主图：ExtendFig1_subtype_umap_composition.png/pdf；独立 panel：ExtendFig1_A...F_*.png/pdf。tables 包含全量 counts、fractions、long tables、n_cells/n_samples/n_studies、Celltype 对照、绘图 ID/坐标及原始行号。data_provenance.json 记录数据核验；color_map.json 记录原配色；layout_validation.json / validation_report.json 记录版式和输出检查。

```bash
cd {HERE}
external_inputs/python -u make_extendfig1_subtype.py
```

后台运行：`bash run_extendfig1_subtype.sh`（内部 nohup），返回 PID 和日志位置。`--prepare-only` 仅核查并准备数据；`--preview-only` 快速预览；`--refresh` 从真实源对象重新核验。默认只重用经过指纹核验的绘图缓存，不修改任何输入或主图结果。Notebook 默认查看已产出结果，REBUILD=False。
''')
    nb=nbf.v4.new_notebook(metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
    nb.cells=[nbf.v4.new_markdown_cell('# Extend Fig1 — subtype UMAP and composition\n\nUMAP 仅显示当前 subtype 并独立放大；Celltype 全名在图内以引线标注，使用源对象配色。C–F 沿用前次版式；字号维持初版 150%，纯黑且不加粗。详见 README。'),
        nbf.v4.new_code_cell(f'from pathlib import Path\nimport pandas as pd,json\nfrom IPython.display import display,Image,Markdown\nOUT=Path({str(out)!r})'),
        nbf.v4.new_code_cell('display(json.loads((OUT/"validation_report.json").read_text()))\ndisplay(Markdown((OUT/"README.md").read_text()))'),
        nbf.v4.new_code_cell(f'display(Image(filename=str(OUT/{(MAIN+".png")!r}),width=1800))'),
        nbf.v4.new_code_cell('for lineage in ["fibro","endo"]:\n    for group in ["tissue","agegroup"]:\n        display(pd.read_csv(OUT/"tables"/f"{lineage}_{group}_composition_long.csv"))'),
        nbf.v4.new_code_cell('REBUILD=False\nif REBUILD:\n    import subprocess,sys\n    subprocess.run([sys.executable,str(OUT/"make_extendfig1_subtype.py")],cwd=OUT,check=True)')]
    nbf.validate(nb);nbf.write(nb,out/'ExtendFig1_subtype_umap_composition.ipynb')
    fonts=out/'fonts';fonts.mkdir(exist_ok=True)
    font=Path(font_manager.findfont('DejaVu Sans'));shutil.copy2(font,fonts/font.name)
    dump(out/'code_sha256.json',{str(p.relative_to(out)):sha(p) for p in [out/'make_extendfig1_subtype.py',out/'run_extendfig1_subtype.sh',out/'upstream/utils_data.py',out/'input_prompt.md']})
    bundle=out/'ExtendFig1_subtype_download.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
        for path in out.rglob('*'):
            if not path.is_file() or path==bundle or path.name=='prepared_subtype_plot_data.pkl':continue
            if {'tmp','__pycache__','logs'}&set(path.relative_to(out).parts):continue
            z.write(path,path.relative_to(out))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fig1-dir',type=Path,default=FIG1)
    p.add_argument('--output-dir',type=Path,default=HERE)
    p.add_argument('--dpi',type=int,default=300)
    p.add_argument('--plot-cells',type=int,default=None)
    p.add_argument('--refresh',action='store_true')
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--preview-only',action='store_true')
    args=p.parse_args()
    if args.dpi<1 or (args.plot_cells is not None and args.plot_cells<0):p.error('Invalid dpi/cell limit')
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);(out/'tables').mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(out/'run.log')],force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':FONT_SIZES['base'],
        'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
        'font.weight':'normal','axes.labelweight':'normal','axes.titleweight':'normal',
        'pdf.fonttype':42,'ps.fonttype':42,'savefig.bbox':None,'savefig.facecolor':'white'})
    data,cfg,key=prepare(args.fig1_dir.resolve(),out,args.plot_cells,args.refresh)
    load_celltype_colors(data,cfg,out)
    if args.prepare_only:logging.info('Preparation complete');return
    render(data,cfg,out,args.dpi,preview=args.preview_only)
    if args.preview_only:logging.info('Preview complete');return
    validate(out,data,cfg,args.dpi)
    documents(out,data,cfg,key,args.dpi)
    logging.info('SUCCESS: %s',out)


if __name__=='__main__':main()
