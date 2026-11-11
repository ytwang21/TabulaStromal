"""Common visual grammar; vector labels/lines, rasterized million-cell scatter layers."""
import json
import logging
from pathlib import Path
import string

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from utils_data import ORDER, COLORS, MARKERS, AGE_LABELS, age_order, stratified_rows, dump_json

logging.getLogger('fontTools').setLevel(logging.WARNING)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,
                     'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
TITLES={'fibro':'Fibroblast','endo':'Endothelial'}


def save(fig, stem, dpi=300):
    stem=Path(stem);stem.parent.mkdir(exist_ok=True,parents=True)
    fig.savefig(str(stem)+'.png',dpi=dpi,bbox_inches='tight',pad_inches=.12)
    fig.savefig(str(stem)+'.pdf',dpi=dpi,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)


def palette(labels, role, info):
    labels=list(labels)
    if role=='general':return {k:COLORS[k] for k in labels}
    previous=info.get('palettes',{}).get(role,{})
    cmap=plt.get_cmap('tab20' if len(labels)<=20 else 'gist_ncar')
    return {k:previous.get(k,matplotlib.colors.to_hex(cmap(i/max(1,len(labels)-1)))) for i,k in enumerate(labels)}


def outlines(ax, xy, values, colors, gridsize=220, mass=.90):
    """Smoothed density contours retain disjoint islands, with one shared color/state.

    A 90% highest-density region is an annotation aid, not a confidence interval.
    Small groups (<15 sampled cells) receive no unreliable contour.
    """
    limits=np.quantile(xy,[0,.9999],axis=0)
    span=np.maximum(limits[1]-limits[0],1e-4)
    low=limits[0]-.05*span; high=limits[1]+.05*span
    edges=[np.linspace(low[d],high[d],gridsize+1) for d in [0,1]]
    centers=[(a[:-1]+a[1:])/2 for a in edges]
    paths={}
    for group,color in colors.items():
        pts=xy[np.asarray(values)==group]
        if len(pts)<15:continue
        h,_,_=np.histogram2d(pts[:,0],pts[:,1],bins=edges)
        z=gaussian_filter(h.T,1.8)
        zz=np.sort(z.ravel())[::-1];cum=np.cumsum(zz)
        level=zz[min(np.searchsorted(cum,mass*cum[-1]),len(zz)-1)]
        if not 0<level<z.max():continue
        c=ax.contour(centers[0],centers[1],z,levels=[level],colors=[color],
                     linewidths=.9,linestyles='--',alpha=.95)
        paths[group]=[s.tolist() for s in c.allsegs[0] if len(s)>=4]
    return paths


def atlas_panel(ax, obs, xy, lineage, info, role='general', dashed=True, legend=True):
    labels=ORDER[lineage] if role=='general' else sorted(obs[role].unique(),key=lambda x:(ORDER[lineage].index(x.split('_')[0]) if x.split('_')[0] in ORDER[lineage] else 99,x))
    colors=palette(labels,role,info)
    # Deterministic shuffle prevents late-drawn populations systematically covering others.
    rng=np.random.default_rng(42);idx=rng.permutation(len(obs))
    ax.scatter(xy[idx,0],xy[idx,1],c=[colors[x] for x in obs[role].to_numpy()[idx]],
               s=.7 if len(obs)>30000 else 2,alpha=.55,lw=0,rasterized=True)
    paths=outlines(ax,xy,obs[role].to_numpy(),colors) if dashed else {}
    if role=='general':
        for k in labels:
            pts=xy[obs[role].to_numpy()==k]
            if len(pts):
                med=np.median(pts,axis=0)
                ax.text(*med,k,fontsize=8,fontweight='bold',ha='center',va='center',
                        bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=.8))
    if legend:
        ax.legend(handles=[Line2D([],[],marker='o',ls='',color=colors[k],label=k,markersize=4) for k in labels],
                  loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=6.5,
                  ncol=2 if len(labels)>22 else 1,handletextpad=.25,labelspacing=.4)
    ax.set(xticks=[],yticks=[],xlabel='',ylabel='')
    for s in ax.spines.values():s.set_visible(False)
    ax.set_aspect('equal',adjustable='datalim')
    ax.annotate('',xy=(.16,.04),xytext=(.04,.04),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.7))
    ax.annotate('',xy=(.04,.16),xytext=(.04,.04),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.7))
    ax.text(.17,.04,'UMAP1',transform=ax.transAxes,fontsize=6,va='center')
    ax.text(.04,.18,'UMAP2',transform=ax.transAxes,fontsize=6,ha='center')
    ax.set_title(f'{TITLES[lineage]} | {info["n_cells"]:,} cells',loc='left')
    return paths


def summarize(obs, key, top=None):
    labels=obs[key].astype(str)
    if top:
        keep=labels.value_counts().head(top).index
        labels=labels.where(labels.isin(keep),'Others')
    frame=pd.DataFrame({'group':labels.to_numpy(),'sample':obs['sample'].to_numpy()})
    n=frame.groupby('group',observed=True).size().rename('n_cells')
    ns=frame[frame['sample']!='Unknown'].groupby('group',observed=True)['sample'].nunique().rename('n_samples')
    table=pd.concat([n,ns],axis=1).fillna(0).astype(int)
    order=age_order(table.index) if key=='age_group' else table.sort_values('n_cells',ascending=False).index.tolist()
    if 'Others' in order:order=[x for x in order if x!='Others']+['Others']
    return table.reindex(order)


def summary_panel(fig, spec, tissue, age, lineage, xmax=None):
    grid=spec.subgridspec(2,2,width_ratios=[1,.8],height_ratios=[1.4,1],wspace=.65,hspace=.65)
    axes=[]
    for i,(table,kind) in enumerate([(tissue,'Tissue'),(age,'Age group')]):
        for j,(col,title,color) in enumerate([('n_cells','Cells','#72A9BB'),('n_samples','Samples','#C4A06F')]):
            ax=fig.add_subplot(grid[i,j]);axes.append(ax)
            ax.barh(np.arange(len(table)),table[col],color=color,height=.7)
            labs=[AGE_LABELS.get(v,v) for v in table.index] if kind=='Age group' else list(table.index)
            ax.set_yticks(np.arange(len(table)),labs if j==0 else ['']*len(labs))
            ax.invert_yaxis();ax.set_xlabel(title);ax.tick_params(axis='y',length=0)
            ax.spines['left'].set_visible(False)
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,4),useMathText=True)
            if xmax is not None:ax.set_xlim(0,max(1,xmax[i,j])*1.18)
            else:ax.set_xlim(0,max(1,table[col].max())*1.23)
            for y,v in enumerate(table[col]):
                ax.text(v+.01*ax.get_xlim()[1],y,f'{v/1e6:.2f}M' if v>=1e6 else f'{v/1e3:.1f}k' if v>=10000 else str(v),va='center',fontsize=5.7)
            if i==0:ax.set_title(f'{TITLES[lineage]} · {title.lower()}',fontsize=9,loc='left')
    return axes


def atlas_stage(lineage, config, root):
    out=root/lineage
    obs=pd.read_pickle(out/'metadata.pkl');xy=np.load(out/'umap.npy',mmap_mode='r')
    info=json.loads((root/'shared'/f'data_audit_{lineage}.json').read_text())
    rows=stratified_rows(obs,config['plot_cells'],config['seed'],('subtype','tissue'))
    np.save(out/'plot_rows.npy',rows)
    for role,suffix,dashed in [('celltype','celltype',False),('subtype','subtype_dashed',True),('general','general_dashed',True),('tissue','tissue',False)]:
        fig,ax=plt.subplots(figsize=(9.5 if role in ['subtype','celltype'] else 7.5,5.5))
        paths=atlas_panel(ax,obs.iloc[rows],xy[rows],lineage,info,role,dashed)
        if dashed:dump_json(out/f'{lineage}_{suffix}_outlines.json',dict(method='90% smoothed highest-density contours on deterministic stratified plot sample; disconnected contours share colors; no confidence/trajectory interpretation',sample_n=len(rows),paths=paths))
        save(fig,out/f'{lineage}_umap_{suffix}',config['dpi'])
    tissue=summarize(obs,'tissue');age=summarize(obs,'age_group')
    tissue.to_csv(out/'summary_tissue.csv');age.to_csv(out/'summary_age_group.csv')
    top=summarize(obs,'tissue',config['top_tissues']);top.to_csv(out/'summary_tissue_top.csv')
    fig=plt.figure(figsize=(9,9));summary_panel(fig,fig.add_gridspec(1,1)[0],top,age,lineage)
    fig.text(.05,.005,'Sample counts exclude unknown IDs; a sample can occur in multiple categories.\nAge A/H retain historical mixed numeric/text bins; Unknown is displayed explicitly.',fontsize=7)
    save(fig,out/f'{lineage}_dataset_summary',config['dpi'])


def marker_matrices(table, lineage, groups=None):
    """Shared all-cell marker calculations, used by the original atlas and Fig1."""
    groups = groups if groups is not None else ORDER[lineage]
    genes=[g for k in groups for g in sorted(MARKERS[lineage][k]) if g in set(table.gene)]
    mean=table.pivot(index='subtype',columns='gene',values='mean_expression').reindex(index=groups,columns=genes)
    frac=table.pivot(index='subtype',columns='gene',values='fraction_expressing').reindex(index=groups,columns=genes)
    color=(mean-mean.min(axis=0))/(mean.max(axis=0)-mean.min(axis=0)).replace(0,1)
    return genes, mean, frac, color


def signature_matrices(means, lineage, groups=None):
    """Return unchanged raw state means and the existing column-z-score display."""
    groups = groups if groups is not None else ORDER[lineage]
    means=means.reindex(index=groups,columns=groups)
    z=(means-means.mean(axis=0))/means.std(axis=0,ddof=0).replace(0,1)
    return means, z


def marker_panel(ax, table, lineage, fig, colorbar=True):
    groups=ORDER[lineage]
    genes,mean,frac,color=marker_matrices(table,lineage,groups)
    xs,ys=np.meshgrid(np.arange(len(genes)),np.arange(len(groups)))
    dots=ax.scatter(xs.ravel(),ys.ravel(),s=frac.to_numpy().ravel()*95,
                    c=color.to_numpy().ravel(),cmap='Blues',vmin=0,vmax=1,lw=.3,edgecolors='#A0B3BE')
    ax.set_xticks(np.arange(len(genes)),genes,rotation=60,ha='right',rotation_mode='anchor')
    ax.set_yticks(np.arange(len(groups)),groups);ax.set_ylim(len(groups)-.5,-1.5)
    start=0
    for k in groups:
        ng=sum(g in genes for g in MARKERS[lineage][k]);end=start+ng
        if ng:
            ax.plot([start-.35,end-.65],[-.8,-.8],color=COLORS[k],lw=2,clip_on=False)
            ax.text((start+end-1)/2,-1.05,k,color=COLORS[k],ha='center',fontsize=7)
        start=end
    ax.set_title(f'{TITLES[lineage]} | marker expression',loc='left',pad=18)
    ax.tick_params(length=0)
    for s in ax.spines.values():s.set_visible(False)
    if colorbar:
        fig.colorbar(dots,ax=ax,fraction=.025,pad=.02,label='Mean log expression (gene-scaled)')
    handles=[plt.scatter([],[],s=f*95,color='#C2D8E5',edgecolors='#A0B3BE',lw=.3,label=f'{int(f*100)}%') for f in [.1,.5,1]]
    ax.legend(handles=handles,title='Expressing',frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.48),fontsize=6,title_fontsize=7,handletextpad=.3,columnspacing=.8)


def signature_panel(ax, means, lineage, fig):
    means,z=signature_matrices(means,lineage)
    im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_xticks(range(len(means.columns)),means.columns,rotation=30,ha='right')
    ax.set_yticks(range(len(means)),means.index)
    ax.set_xlabel('Marker signature');ax.set_ylabel('Annotated state')
    ax.set_title(f'{TITLES[lineage]} | signature enrichment',loc='left')
    for i in range(len(means)):
        for j in range(len(means.columns)):
            ax.text(j,i,f'{means.iloc[i,j]:.2f}',fontsize=6,ha='center',va='center',color='white' if abs(z.iloc[i,j])>1.4 else '#222222')
    fig.colorbar(im,ax=ax,fraction=.035,pad=.025,label='Column z-score of state means')


def marker_stage(lineage,config,root):
    out=root/lineage
    table=pd.read_csv(out/'marker_summary.csv');means=pd.read_csv(out/'signature_mean.csv',index_col=0)
    fig,ax=plt.subplots(figsize=(11,3.8 if lineage=='fibro' else 4.5))
    marker_panel(ax,table,lineage,fig);save(fig,out/f'{lineage}_marker_dotplot',config['dpi'])
    fig,ax=plt.subplots(figsize=(6,4.3));signature_panel(ax,means,lineage,fig)
    save(fig,out/f'{lineage}_signature_heatmap',config['dpi'])
    # Separate panels for each signature using sample-level state means, without cell-level p-values.
    sample=pd.read_csv(out/'signature_sample_means.csv')
    groups=ORDER[lineage]
    fig,axes=plt.subplots(2,2 if len(groups)==4 else 3,figsize=(11,6.2),layout='constrained')
    rng=np.random.default_rng(config['seed'])
    for ax,signature in zip(axes.ravel(),groups):
        t=sample[sample.signature==signature]
        data=[t.loc[t.subtype==k,'mean'].to_numpy() for k in groups]
        data=[a if len(a) else np.array([np.nan]) for a in data]
        ax.boxplot(data,showfliers=False,widths=.55,medianprops=dict(color='#222222'))
        for i,(k,a) in enumerate(zip(groups,data)):
            display=rng.choice(a,min(200,len(a)),replace=False)
            ax.scatter(i+1+rng.uniform(-.2,.2,len(display)),display,s=2,c=COLORS[k],alpha=.3,rasterized=True)
        ax.set_xticks(range(1,len(groups)+1),groups,rotation=30,ha='right');ax.axhline(0,color='.7',lw=.5)
        ax.set_title(f'{signature} signature');ax.set_ylabel('Sample × state mean score')
    save(fig,out/f'{lineage}_signature_sample_boxplots',config['dpi'])
    obs=pd.read_pickle(out/'metadata.pkl');rows=np.load(out/'plot_rows.npy') if (out/'plot_rows.npy').exists() else stratified_rows(obs,config['plot_cells'],config['seed'],('subtype','tissue'))
    xy=np.load(out/'umap.npy',mmap_mode='r')[rows];scores=np.load(out/'signature_scores.npy',mmap_mode='r')
    fig,axes=plt.subplots(2,2 if len(groups)==4 else 3,figsize=(11,7),layout='constrained')
    for j,(ax,k) in enumerate(zip(axes.ravel(),groups)):
        val=scores[rows,j];vmin,vmax=np.quantile(val,[.02,.98])
        im=ax.scatter(xy[:,0],xy[:,1],c=val,s=.7,lw=0,cmap='viridis',vmin=vmin,vmax=vmax,rasterized=True)
        ax.set(xticks=[],yticks=[],title=f'{k} signature');ax.set_axis_off()
        fig.colorbar(im,ax=ax,fraction=.03,pad=.01)
    save(fig,out/f'{lineage}_signature_umap',config['dpi'])


def pseudotime_panel(ax,frame,lineage,fig):
    finite=np.isfinite(frame.pseudotime)
    ax.scatter(frame.umap1,frame.umap2,c='#D9DDE0',s=1,lw=0,rasterized=True)
    if finite.any():
        t=frame.loc[finite]
        im=ax.scatter(t.umap1,t.umap2,c=t.pseudotime,s=2,lw=0,cmap='viridis',vmin=0,vmax=1,rasterized=True)
        fig.colorbar(im,ax=ax,fraction=.03,pad=.02,label='Root-relative DPT')
    ax.set(xticks=[],yticks=[]);ax.set_axis_off()
    ax.set_title(f'{TITLES[lineage]} | state exploration',loc='left')
    ax.text(.01,-.06,f'Finite DPT: {finite.sum():,}/{len(frame):,} sampled cells; grey = outside root component',transform=ax.transAxes,fontsize=6)


def paga_panel(ax, frame, matrix, lineage, threshold=.15):
    values=frame['general'].to_numpy()
    ax.scatter(frame.umap1,frame.umap2,c=[COLORS[k] for k in values],s=1.4,lw=0,alpha=.35,rasterized=True)
    categories=list(matrix.index)
    centers={k:frame.loc[frame['general']==k,['umap1','umap2']].median().to_numpy() for k in categories}
    for i,a in enumerate(categories):
        for j,b in enumerate(categories):
            weight=float(matrix.loc[a,b])
            if i<j and weight>=threshold:
                x,y=centers[a],centers[b]
                ax.plot([x[0],y[0]],[x[1],y[1]],color='#5A6570',lw=.4+2*weight,alpha=.7,zorder=3)
    for k,v in centers.items():
        ax.scatter(*v,s=100,color=COLORS[k],edgecolors='white',lw=.8,zorder=4)
        ax.text(*v,k,fontsize=7,ha='center',va='bottom',zorder=5)
    ax.set_axis_off();ax.set_title(f'{TITLES[lineage]} | undirected PAGA',loc='left')
    ax.text(.01,-.06,f'Exploratory state relationships; edges ≥ {threshold}; DPT provided separately',transform=ax.transAxes,fontsize=6)


def _integrated_figure_v1(config,root):
    out=root/'integrated_figure';out.mkdir(exist_ok=True)
    # Draw original artists into a common layout; no bitmap panel stitching.
    fig=plt.figure(figsize=(22,17.8),layout=None)
    outer=fig.add_gridspec(5,2,height_ratios=[1.8,2.1,1.1,1.05,1.6],hspace=.60,wspace=.43,
                          left=.07,right=.93,top=.965,bottom=.07)
    xmax=np.zeros((2,2))
    for lin in ['fibro','endo']:
        for i,file in enumerate(['summary_tissue_top.csv','summary_age_group.csv']):
            t=pd.read_csv(root/lin/file,index_col=0)
            xmax[i]=np.maximum(xmax[i],t[['n_cells','n_samples']].max().to_numpy())
    conclusions=[];selected_panels={}
    for col,lin in enumerate(['fibro','endo']):
        p=root/lin;obs=pd.read_pickle(p/'metadata.pkl');rows=np.load(p/'plot_rows.npy')
        xy=np.load(p/'umap.npy',mmap_mode='r')[rows]
        info=json.loads((root/'shared'/f'data_audit_{lin}.json').read_text())
        ax=fig.add_subplot(outer[0,col]);atlas_panel(ax,obs.iloc[rows],xy,lin,info,legend=False)
        ax.text(-.10,1.08,string.ascii_uppercase[col],transform=ax.transAxes,fontsize=15,fontweight='bold')
        tissue=pd.read_csv(p/'summary_tissue_top.csv',index_col=0);age=pd.read_csv(p/'summary_age_group.csv',index_col=0)
        axes=summary_panel(fig,outer[1,col],tissue,age,lin,xmax)
        axes[0].text(-.30,1.14,string.ascii_uppercase[2+col],transform=axes[0].transAxes,fontsize=15,fontweight='bold')
        ax=fig.add_subplot(outer[2,col]);marker_panel(ax,pd.read_csv(p/'marker_summary.csv'),lin,fig)
        ax.text(-.10,1.21,string.ascii_uppercase[4+col],transform=ax.transAxes,fontsize=15,fontweight='bold')
        ax=fig.add_subplot(outer[3,col]);signature_panel(ax,pd.read_csv(p/'signature_mean.csv',index_col=0),lin,fig)
        ax.text(-.10,1.13,string.ascii_uppercase[6+col],transform=ax.transAxes,fontsize=15,fontweight='bold')
        decision_path=p/'trajectory'/'trajectory_selection.json'
        decision=json.loads(decision_path.read_text()) if decision_path.exists() else {'main_panel':'dpt'}
        selected_panels[lin]=decision.get('model',decision['main_panel'])
        ax=fig.add_subplot(outer[4,col])
        if config.get('trajectory_engine')=='metacell_three_method':
            from metacell_trajectory import draw_selected
            draw_selected(ax,lin,config,root)
        elif decision['main_panel']=='paga':
            frame=pd.read_csv(p/'trajectory'/'pseudotime.csv.gz')
            matrix=pd.read_csv(p/'trajectory'/f'paga_{decision["reference"]}.csv',index_col=0)
            paga_panel(ax,frame,matrix,lin,config['paga_threshold'])
        else:
            frame=pd.read_csv(p/'trajectory'/'pseudotime.csv.gz');pseudotime_panel(ax,frame,lin,fig)
        ax.text(-.10,1.10,string.ascii_uppercase[8+col],transform=ax.transAxes,fontsize=15,fontweight='bold')
        top=tissue.drop(index='Others',errors='ignore').n_cells.idxmax()
        unknown=info['missing']['age_group']/info['n_cells']
        concord=json.loads((p/'signature_concordance.json').read_text())['top_score_agreement']
        availability=json.loads((p/'marker_availability.json').read_text())
        missing=[g for genes in availability['missing'].values() for g in genes]
        missing_note='; missing markers: '+', '.join(missing) if missing else '; all requested markers present'
        conclusions.append(f'{TITLES[lin]}: {top} contributes most cells; age-group Unknown {unknown:.1%}; marker top-score agreement {concord:.1%} (descriptive)'+missing_note+'.')
    trajectory_note='Trajectory: a common recorded meta-cell cohort for Monocle 2, STREAM and CellRank; each meta-cell aggregates original cells (median 30, variable size).' if config.get('trajectory_engine')=='metacell_three_method' else 'Trajectory uses stratified cells and existing PCA; roots are assumptions.'
    footer='\n'.join(conclusions+['Dashed contours: 90% density regions on plot samples. Cells/samples use full metadata; marker summaries use all cells.',
        'Age A/H: historical bins include text labels; A includes embryo/infant and H includes “Aged”. Sample counts are not donor counts.',
        trajectory_note,'Roots/direction are assumptions; visual branches, pseudotime and transition arrows do not establish a uniquely correct developmental lineage.'])
    fig.text(.07,.011,footer,fontsize=7,linespacing=1.5)
    save(fig,out/'Fig_StromalAtlas_Overview',config['dpi'])
    trajectory_caption='Selected native Monocle/STREAM/CellRank fits on user-provided meta-cells; full three-method panels and all Monocle feature settings are separate. No PAGA.' if config.get('trajectory_engine')=='metacell_three_method' else 'Selected exploratory state relationship panels.'
    (out/'figure_caption.txt').write_text('A–B: broad states in existing UMAP; fine Celltype_sub and Celltype UMAPs are provided separately.\nC–D: full tissue/age cell and unique-sample coverage, shared horizontal scales across lineages.\nE–F: grouped, alphabetically sorted marker dotplots; size = detected fraction, color = per-gene min–max state mean log expression.\nG–H: matched-control marker signatures; numerals = raw state means, color = column z-score. Same marker sets are descriptive support, not independent annotation validation.\nI–J: '+trajectory_caption+' '+json.dumps(selected_panels)+'\n'+footer+'\n')
    (root/'shared'/'core_findings.txt').write_text('\n'.join(conclusions)+'\nRead trajectory_report.md for graph/root robustness before interpreting pseudotime.\n')
    if config.get('trajectory_engine')=='metacell_three_method':
        from metacell_trajectory import draw_model
        fig,axes=plt.subplots(2,3,figsize=(16,10))
        for row,lin in enumerate(['fibro','endo']):
            p=root/lin/'trajectory';ob=pd.read_pickle(p/'prepared/metadata.pkl')
            choice=json.loads((p/'trajectory_selection.json').read_text())
            for col,fam in enumerate(['monocle','stream','cellrank']):
                draw_model(axes[row,col],p/choice['best_by_family'][fam]['path'],ob,lin)
                axes[row,col].text(-.08,1.08,string.ascii_uppercase[3*row+col],transform=axes[row,col].transAxes,fontsize=13,fontweight='bold')
        fig.text(.04,.015,'Same recorded meta-cell cohort per lineage. Black edges = fitted principal tree; CellRank arrows = transition-matrix displacement. Stars = assumed roots. No developmental ground truth is available.',fontsize=8)
        save(fig,out/'Fig_StromalAtlas_Trajectory_Comparison',config['dpi'])


def integrated_figure(config, root):
    """Current overview: reference Leiden atlas, merged coverage, verified meta-cell inputs."""
    from overview_v2 import render_overview
    render_overview(config, root, options=config.get('overview_plot'))
    # Keep the independent 2x3 three-method comparison in full pipeline deliveries.
    from metacell_trajectory import draw_model
    out = root / 'integrated_figure'
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for row, lineage in enumerate(['fibro', 'endo']):
        path = root / lineage / 'trajectory'
        obs = pd.read_pickle(path / 'prepared/metadata.pkl')
        choice = json.loads((path / 'trajectory_selection.json').read_text())
        for col, family in enumerate(['monocle', 'stream', 'cellrank']):
            draw_model(axes[row, col], path / choice['best_by_family'][family]['path'], obs, lineage)
            stamp = string.ascii_uppercase[3 * row + col]
            axes[row, col].text(-.08, 1.08, stamp, transform=axes[row, col].transAxes,
                               fontsize=13, fontweight='bold')
    fig.text(.04, .015, 'Common original meta-cells per lineage. Stars = assumed roots. CellRank arrows use transition-matrix displacement; direction is a prior.', fontsize=8)
    save(fig, out / 'Fig_StromalAtlas_Trajectory_Comparison', config['dpi'])
