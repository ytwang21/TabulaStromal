#!/usr/bin/env python3
"""Rebuild Fig2 workflow, fold-4 heatmaps and final stable-feature training fits; no training."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fig2-mpl')
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / 'AgePrediction/GeneModel/TabulaStromal_Rebenchmark_2026/Stromal_Combined_Figure/outputs'
ORDER = ['XGBoost', 'ExtraTrees', 'ElasticNet', 'MLP', 'Ridge', 'Lasso',
         'RandomForest', 'GradientBoosting', 'Supervised_PCA_Ridge', 'SVR',
         'DeepSets', 'Pseudobulk', 'Ensemble']
LINEAGES = ['Fibroblast', 'Endothelial']
SPLITS = ['Train', 'Valid', 'External']
COLORS = dict(Train='#A3ADB8', Valid='#D48736', External='#00867B')
FONT = 'DejaVu Sans'
LINEAGE_COLORS = {'Fibroblast':'#DB4854', 'Endothelial':'#2474B5'}
STABLE_ROOT = ROOT/'AgePrediction/GeneModel/TabulaStromal_Rebenchmark_2026'
STABLE_PATHS = {'Fibroblast':STABLE_ROOT/'Fibro_Accelerated_20261102/results/Fibro/stable_clock_optimization',
                'Endothelial':STABLE_ROOT/'Endo/stable_clock_optimization'}
plt.rcParams.update({'font.family': FONT, 'font.size': 12, 'pdf.fonttype': 42,
                     'ps.fonttype': 42, 'axes.linewidth': .8, 'text.color': '#162A3B',
                     'axes.labelcolor': '#162A3B', 'xtick.color': '#374756',
                     'ytick.color': '#374756', 'savefig.facecolor': 'white'})


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def metrics(d):
    y, p = d.Age.to_numpy(float), d.prediction.to_numpy(float)
    assert len(y) and np.isfinite(y).all() and np.isfinite(p).all()
    return dict(n=len(y), R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum(),
                MAE=np.abs(y-p).mean(), MAD=np.abs(y-y.mean()).mean())


def prepare(out):
    data = out / 'data'
    data.mkdir(exist_ok=True)
    provenance = []
    def register(p):
        provenance.append(dict(path=str(p.resolve()), sha256=digest(p), bytes=p.stat().st_size))
    workflow = HERE/'overview/TabulaStromal_Workflow_V6.pdf'
    register(workflow)
    register(SOURCE/'Stromal_individual_aging_clock_Fibroblast_Endothelial_4.pdf')
    register(SOURCE/'completion_4.json')
    tables, points, scores, checks = [], [], [], []
    for lineage in LINEAGES:
        frames = {}
        for kind in ['benchmark_summary', 'benchmark_per_fold', 'scatter_predictions',
                     'scatter_metric_diagnostics', 'external_metric_scope_audit', 'internal_K_selection']:
            p = SOURCE/'plot_data'/f'{lineage}_{kind}_4.csv'
            register(p)
            shutil.copy2(p, data/p.name)
            frames[kind] = pd.read_csv(p)
        summary = frames['benchmark_summary']
        t = frames['benchmark_per_fold'].copy()
        assert len(t) == 39 and set(t.model) == set(ORDER)
        assert t.fold.eq(4).all()
        t['lineage'] = lineage
        t['original_split'] = t.split
        t['split'] = t.split.map(lambda x: 'External' if x.startswith('External') else x)
        assert not t.duplicated(['model','split']).any()
        assert set(t.split) == set(SPLITS)
        for _, r in t.iterrows():
            s = summary[(summary.model==r.model)&(summary.K==r.K)&(summary.split==r.original_split)]
            assert len(s)==1
            for m in ['R2','MAE','MAD']:
                np.testing.assert_allclose(r[m], s.iloc[0][m+'_mean'], rtol=1e-11, atol=1e-11)
            if r.split == 'External':
                a = frames['external_metric_scope_audit']
                a = a[(a.model==r.model)&(a.K==r.K)&(a.fold==4)]
                assert len(a)==1 and bool(a.iloc[0].same_pairs_for_all_metrics)
                for m in ['n','R2','MAE','MAD']:
                    np.testing.assert_allclose(r[m], a.iloc[0][m], rtol=1e-11, atol=1e-11)
            checks.append(dict(lineage=lineage, model=r.model, split=r.split,
                               summary_match=True, external_audit_match=r.split=='External'))
        tables.append(t[['lineage','model','K','fold','split','n','R2','MAE','MAD','evaluation_scope']])
        p = frames['scatter_predictions'].copy()
        p['lineage'] = lineage
        p['split'] = p.display_group.map(lambda x: 'External' if x.startswith('External') else x)
        assert not p.sample_id.duplicated().any()
        assert p.model.eq('ElasticNet').all() and p.K.eq(20000).all() and p.fold.eq(4).all()
        for split in SPLITS:
            sub = p[p.split==split]
            s = metrics(sub)
            r = t[(t.model=='ElasticNet')&(t.split==split)].iloc[0]
            for m in ['n','R2','MAE','MAD']:
                # Archived training metrics and CSV predictions differ at floating-point precision.
                np.testing.assert_allclose(s[m],r[m],rtol=1e-8,atol=1e-8)
            scores.append(dict(lineage=lineage,model='ElasticNet',K=20000,fold=4,split=split,**s))
        points.append(p[['lineage','sample_id','model','K','fold','split','Age','prediction',
                         'is_adaptation','is_heldout','evaluation_scope']])
    t, p, s = pd.concat(tables,ignore_index=True), pd.concat(points,ignore_index=True), pd.DataFrame(scores)
    assert t[['R2','MAE','MAD']].notna().all().all()
    t.to_csv(data/'Fig2_performance_long.csv',index=False)
    p.to_csv(data/'Fig2_fold4_scatter_predictions.csv',index=False)
    s.to_csv(data/'Fig2_fold4_scatter_metrics.csv',index=False)
    pd.DataFrame(checks).to_csv(data/'Fig2_metric_validation.csv',index=False)
    for m in ['R2','MAE','MAD']:
        t.pivot(index='model',columns=['lineage','split'],values=m).reindex(
            index=ORDER,columns=pd.MultiIndex.from_product([LINEAGES,SPLITS])).to_csv(data/f'Fig2_{m}_matrix.csv')
    (data/'source_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # Retain the original 1536x1024 baseline after the public V6 files are
    # synchronized. Re-extracting their new composite PDF would lose overlays
    # and change the pixel coordinate system used by workflow().
    baseline = HERE/'data/workflow_original_extracted.png'
    if baseline.exists():
        original=Image.open(baseline).convert('RGB')
        register(baseline)
        (data/'source_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    else:
        original=PdfReader(workflow).pages[0].images[0].image.convert('RGB')
    assert original.size==(1536,1024),'Workflow drawing requires the archived original 1536x1024 baseline.'
    if data/'workflow_original_extracted.png' != baseline or not baseline.exists():
        original.save(data/'workflow_original_extracted.png')
    return t,p,s,original


def stable_metrics(d):
    y=d.Age.to_numpy(float);p=d.prediction.to_numpy(float);err=np.abs(y-p)
    assert np.isfinite(y).all() and np.isfinite(p).all()
    return dict(n=len(d),R2=1-np.square(y-p).sum()/np.square(y-y.mean()).sum(),
                MAE=err.mean(),MAD=np.median(err),PCC=np.corrcoef(y,p)[0,1])


def load_stable_results(out):
    """Use the exact cited training-fit and selected ALL100 prediction tables."""
    folder=HERE/'benchmark/final_model'
    inputs=[folder/'ElasticNet_final_fit_predictions.csv',folder/'ElasticNet_final_fit_metrics.csv',
            folder/'ElasticNet_final_fit.png',folder/'ElasticNet_final_fit.pdf',
            STABLE_ROOT/'StableClock_Optimization/balanced_ALL100_selected_summary.csv',
            STABLE_ROOT/'StableClock_Optimization/Stromal_balanced_ALL100_selected_clock.png']
    train=pd.read_csv(inputs[0]);ts=pd.read_csv(inputs[1])
    assert train.prediction_scope.eq('development_in_sample').all()
    train['Age']=train.actual_age;train['prediction']=train.predicted_age
    assert not train.duplicated(['lineage','sample_id']).any()
    selected=pd.read_csv(inputs[4]);frames=[];scores=[]
    for lineage,key in [('Fibroblast','Fibro'),('Endothelial','Endo')]:
        d=train[train.lineage.eq(lineage)];r=ts[ts.lineage.eq(lineage)].iloc[0]
        calc=stable_metrics(d)
        for m in ['n','R2','MAE','MAD','PCC']:
            np.testing.assert_allclose(calc[m],r[m],rtol=1e-9,atol=1e-9)
        path=STABLE_PATHS[lineage]/'performance_first_all100_predictions.csv';inputs.append(path)
        q=pd.read_csv(path);q['lineage']=lineage
        assert not q.sample_id.duplicated().any()
        r=selected[selected.lineage.eq(key)].iloc[0];calc=stable_metrics(q)
        assert calc['n']==int(r.n)
        # The saved summary can predate serialization at ~2.3e-6 years.
        for m in ['R2','MAE','MAD']:np.testing.assert_allclose(calc[m],r['ALL100_'+m],rtol=1e-6,atol=3e-6)
        scores.append(dict(lineage=lineage,n=calc['n'],**{m:float(r['ALL100_'+m]) for m in ['R2','MAE','MAD']},
                           prediction_scope='selected_external_ALL100',MAD_definition='median_absolute_prediction_error'))
        frames.append(q)
    external=pd.concat(frames,ignore_index=True);es=pd.DataFrame(scores)
    train.to_csv(out/'data/Fig2_scatter_predictions.csv',index=False)
    ts.to_csv(out/'data/Fig2_scatter_metrics.csv',index=False)
    external.to_csv(out/'data/Fig2_workflow_external_predictions.csv',index=False)
    es.to_csv(out/'data/Fig2_workflow_external_metrics.csv',index=False)
    manifest=out/'data/source_manifest.json';items=json.loads(manifest.read_text())
    items.extend(dict(path=str(p.resolve()),sha256=digest(p),bytes=p.stat().st_size) for p in inputs)
    manifest.write_text(json.dumps(items,indent=2)+'\n')
    return train,ts,external,es


def axes(fig, x,y,w,h):
    fw,fh=fig.get_size_inches()
    return fig.add_axes([x/fw,y/fh,w/fw,h/fh])


def text(fig,x,y,s,**kw):
    fw,fh=fig.get_size_inches()
    return fig.text(x/fw,y/fh,s,**kw)


def heading(fig,letter,title,x,y):
    text(fig,x,y,letter,fontsize=22,weight='bold',va='baseline')
    text(fig,x+.39,y+.01,title,fontsize=16,weight='bold',va='baseline')


def workflow(fig,im,x,y,w,h,train,ts,external,es,include_heading=True):
    """Keep the workflow framework/order, updating only two requested regions."""
    if include_heading:
        heading(fig,'A','StromalClock workflow',x,y+h-.3)
    ax=axes(fig,x,y,w,h-.75 if include_heading else h)
    ax.imshow(im,interpolation='lanczos',aspect='equal')
    ax.axis('off')
    factor=w/15.4
    # Original-pixel coordinates; keep headings, separators and all other art.
    ax.add_patch(Rectangle((848,532),337,180,facecolor='white',edgecolor='none',zorder=3))
    ax.add_patch(Rectangle((1203,518),330,197,facecolor='white',edgecolor='none',zorder=3))
    def inset(left,top,width,height):
        a=ax.inset_axes([left/im.width,(im.height-top-height)/im.height,width/im.width,height/im.height],zorder=5)
        a.tick_params(labelsize=4.6*factor,length=1.4,pad=1)
        a.spines[['top','right']].set_visible(False)
        a.spines[['left','bottom']].set_linewidth(.45)
        return a
    for left,m,title in [(885,'R2','R²'),(987,'MAE','MAE'),(1089,'MAD','MAD')]:
        a=inset(left,592,83,84)
        a.plot([-5,120],[-5,120],ls='--',color='#8D99A1',lw=.5,zorder=0)
        handles=[]
        for lineage,short in [('Fibroblast','Fibro'),('Endothelial','Endo')]:
            q=train[train.lineage.eq(lineage)];r=ts[ts.lineage.eq(lineage)].iloc[0]
            a.scatter(q.Age,q.prediction,s=.65*factor**2,color=LINEAGE_COLORS[lineage],alpha=.42,linewidths=0)
            label=f'{short}: {r[m]:.3f}' if m=='R2' else f'{short}: {r[m]:.2f}'
            handles.append(Line2D([],[],marker='o',linestyle='',color=LINEAGE_COLORS[lineage],markersize=2.6*factor,label=label))
        a.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.14,1.035),frameon=False,
                 fontsize=5.3*factor,handletextpad=.25,handlelength=.7,labelspacing=.18,borderaxespad=0)
        ax.text(left+41.5,542,title,ha='center',va='center',fontsize=7*factor,weight='bold',color='#102332',zorder=6)
        a.set_xlim(-20,120)
        a.set_ylim(-20,120)
        a.set_xticks([0,60,120])
        a.set_yticks([0,60,120])
        a.set_xlabel('Actual age',fontsize=5.2*factor,labelpad=1)
        if m=='R2':a.set_ylabel('Predicted age',fontsize=5.2*factor,labelpad=1)
    for left,lineage in [(1245,'Fibroblast'),(1410,'Endothelial')]:
        a=inset(left,548,108,123);q=external[external.lineage.eq(lineage)];r=es[es.lineage.eq(lineage)].iloc[0]
        a.scatter(q.Age,q.prediction,s=2.1*factor**2,color=LINEAGE_COLORS[lineage],alpha=.8,linewidths=0)
        a.plot([-5,110],[-5,110],ls='--',color='#8D99A1',lw=.55,zorder=0)
        a.set(xlim=(-5,110),ylim=(-5,110));a.set_xticks([0,50,100]);a.set_yticks([0,50,100])
        a.set_title(lineage,fontsize=6.6*factor,color=LINEAGE_COLORS[lineage],weight='bold',pad=3)
        a.text(.035,.97,f'n={int(r.n)}\nR²={r.R2:.3f}\nMAE={r.MAE:.2f}\nMAD={r.MAD:.2f}',
               transform=a.transAxes,ha='left',va='top',fontsize=5.2*factor,linespacing=1.05,color='#172B39')
        a.set_xlabel('Actual age',fontsize=5.4*factor,labelpad=1)
        a.set_ylabel('Predicted age',fontsize=5.4*factor,labelpad=1)


def final_fit_scatter(fig,p,s,lineage,letter,x,y,size):
    """Replot every point from ElasticNet_final_fit, retaining its fit statistics."""
    heading(fig,letter,lineage,x-.15,y+size+.53)
    ax=axes(fig,x,y,size,size);q=p[p.lineage.eq(lineage)];r=s[s.lineage.eq(lineage)].iloc[0]
    a=q.Age.to_numpy();b=q.prediction.to_numpy()
    lo=min(0,int(np.floor(min(a.min(),b.min())/20))*20)
    hi=max(120,int(np.ceil(max(a.max(),b.max())/20))*20)
    dots=ax.scatter(a,b,s=5,alpha=.55,color='#2474B5',edgecolors='none',linewidths=0,zorder=3)
    assert len(dots.get_offsets())==int(r.n)
    xx=np.linspace(lo,hi,300);slope,intercept=np.polyfit(a,b,1)
    ax.plot(xx,xx,'--',color='#6E6E6E',lw=1.2,dashes=(4,3),label='Ideal line, $y=x$')
    ax.plot(xx,slope*xx+intercept,color='#D55E00',lw=1.8,label='EN regression')
    ax.text(.04,.96,f'R² = {r.R2:.3f}\nMAE = {r.MAE:.2f}\nMAD = {r.MAD:.2f}\nPCC = {r.PCC:.3f}\nn = {int(r.n)}',
            transform=ax.transAxes,ha='left',va='top',fontsize=12,linespacing=1.25)
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Chronological age (years)',ylabel='Predicted age (years)')
    ax.set_aspect('equal',adjustable='box');ax.set_xticks(np.arange(lo,hi+1,20));ax.set_yticks(np.arange(lo,hi+1,20))
    ax.grid(True,lw=.4,alpha=.12);ax.set_axisbelow(True)
    ax.legend(loc='lower right',frameon=False,fontsize=10.5,handlelength=2.3)
    return ax


def heatmap(fig,t,metric,letter,x,y,w,h):
    title={'R2':'R²','MAE':'MAE (years)','MAD':'MAD (years)'}[metric]
    heading(fig,letter,title,x-.85,y+h+.81)
    mat=t.pivot(index='model',columns=['lineage','split'],values=metric).reindex(
        index=ORDER,columns=pd.MultiIndex.from_product([LINEAGES,SPLITS])).to_numpy()
    if metric=='R2':
        cmap=plt.get_cmap('RdBu_r')  # negative blue, zero white, positive red
        norm=TwoSlopeNorm(vmin=float(mat.min()),vcenter=0.,vmax=float(mat.max()))
    else:
        cmap=plt.get_cmap('Blues')  # larger values are darker
        norm=Normalize(vmin=float(mat.min()),vmax=float(mat.max()))
    # Two separate axes create a genuine white gutter between the lineages.
    # Both axes retain the same model rows and metric color normalization.
    gap=.26
    block_width=(w-gap)/2
    blocks=[]
    for group,(lineage,col) in enumerate(zip(LINEAGES,['#B84C55','#337CB5'])):
        ax=axes(fig,x+group*(block_width+gap),y,block_width,h)
        ax.set_gid(f'{letter}_{lineage}')
        blocks.append(ax)
        values=mat[:,group*3:group*3+3]
        mesh=ax.imshow(values,aspect='auto',cmap=cmap,norm=norm,interpolation='none')
        for i in range(13):
            for j in range(3):
                rgba=cmap(norm(values[i,j]));lum=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                ax.text(j,i,f'{values[i,j]:.3f}' if metric=='R2' else f'{values[i,j]:.2f}',
                        ha='center',va='center',fontsize=10.5,color='white' if lum<.51 else '#193B46')
        ax.set_xticks(range(3),SPLITS,rotation=45,ha='right',fontsize=11)
        ax.set_yticks(range(13))
        if group==0:
            ax.set_yticklabels([m.replace('Supervised_PCA_Ridge','SPCA–Ridge').replace('GradientBoosting','Grad. boost') for m in ORDER],fontsize=9.5)
            for label in ax.get_yticklabels():
                if label.get_text()=='ElasticNet':
                    label.set_weight('bold')
                    label.set_color('#9B2C83')
                    label.set_fontsize(10)
        else:
            ax.tick_params(labelleft=False)
        ax.tick_params(length=0,pad=5)
        ax.set_xticks(np.arange(-.5,3,1),minor=True)
        ax.set_yticks(np.arange(-.5,13,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=1)
        ax.tick_params(which='minor',length=0)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.add_patch(Rectangle((-.5,-1.6),3,.71,transform=ax.transData,
                              facecolor=col,clip_on=False,linewidth=0))
        ax.text(1,-1.24,lineage,ha='center',va='center',color='white',fontsize=11,weight='bold',clip_on=False)
    cax=axes(fig,x,y-.83,w,.12)
    ticks=[float(mat.min()),-1.,0.,float(mat.max())] if metric=='R2' else np.linspace(mat.min(),mat.max(),4)
    cb=fig.colorbar(mesh,cax=cax,orientation='horizontal',ticks=ticks)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9,length=2,pad=2)
    cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v,pos:f'{v:.2f}'))
    direction='Negative (blue) · 0 · Positive (red)' if metric=='R2' else ('Darker = higher MAE' if metric=='MAE' else 'Darker = higher MAD')
    text(fig,x+w/2,y-1.18,direction,ha='center',fontsize=10,color='#49626C')
    return blocks[0]


def scatter(fig,p,s,lineage,letter,x,y,size):
    heading(fig,letter,lineage,x-.15,y+size+.53)
    text(fig,x,y+size+.22,'ElasticNet · K = 20,000 · source fold 4',fontsize=11,color='#49626C')
    ax=axes(fig,x,y,size,size)
    ax.plot([-25,125],[-25,125],color='#7D8C98',lw=1,ls=(0,(4,3)),zorder=0)
    for split,marker,ms,alpha in [('Train','o',16,.50),('Valid','^',29,.78),('External','s',35,.9)]:
        d=p[(p.lineage==lineage)&(p.split==split)]
        ax.scatter(d.Age,d.prediction,s=ms,marker=marker,color=COLORS[split],alpha=alpha,
                   edgecolor='white' if split=='External' else 'none',linewidth=.35,zorder=SPLITS.index(split)+1)
    ax.set(xlim=(-25,125),ylim=(-25,125),xlabel='Chronological age (years)',ylabel='Predicted biological age (years)')
    ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([0,25,50,75,100,125]);ax.set_yticks([0,25,50,75,100,125])
    ax.spines[['top','right']].set_visible(False)
    ax.grid(color='#E7ECEF',lw=.6,zorder=0)
    ax.set_axisbelow(True)
    handles,labels=[],[]
    for split,marker in zip(SPLITS,['o','^','s']):
        r=s[(s.lineage==lineage)&(s.split==split)].iloc[0]
        labels.append(f'{split}: n={int(r.n)}, R²={r.R2:.3f}, MAE={r.MAE:.2f}, MAD={r.MAD:.2f}')
        handles.append(Line2D([],[],marker=marker,ls='',color=COLORS[split],markersize=6))
    leg=ax.legend(handles,labels,loc='upper left',bbox_to_anchor=(-.025,-.135),frameon=False,
                  fontsize=10.3,handletextpad=.45,labelspacing=.65,borderaxespad=0)
    leg.get_texts()[-1].set_weight('bold')
    leg.get_texts()[-1].set_color(COLORS['External'])
    return ax


def save(fig,out,name,dpi):
    fig.savefig(out/f'{name}.pdf')
    fig.savefig(out/f'{name}.png',dpi=dpi)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    parser.add_argument('--dpi',type=int,default=400)
    args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(exist_ok=True,parents=True)
    t,p,s,im=prepare(out)
    train,ts,external,es=load_stable_results(out)
    # Portrait composition: intact workflow across row 1, B–D across row 2,
    # and two equal square prediction axes across row 3.
    fig=plt.figure(figsize=(16.2,27.2),facecolor='white')
    workflow_width=15.4
    workflow_height=workflow_width*im.height/im.width+.75
    workflow(fig,im,.4,16.0,workflow_width,workflow_height,train,ts,external,es)
    for x,m,l in zip([1.25,6.45,11.65],['R2','MAE','MAD'],['B','C','D']):
        heatmap(fig,t,m,l,x,10.05,3.95,4.45)
    a=final_fit_scatter(fig,train,ts,'Fibroblast','E',1.15,1.8,6.05)
    b=final_fit_scatter(fig,train,ts,'Endothelial','F',9.45,1.8,6.05)
    fig.canvas.draw()
    geometry=[]
    for letter,ax in [('E',a),('F',b)]:
        box=ax.get_window_extent()
        assert abs(box.width-box.height)<1e-6
        geometry.append(dict(panel=letter,width_inches=box.width/fig.dpi,height_inches=box.height/fig.dpi))
    save(fig,out,'Fig2_StromalClock_ModelPerformance',args.dpi)
    panels=out/'panels';panels.mkdir(exist_ok=True)
    f=plt.figure(figsize=(14,10.2));workflow(f,im,.3,.3,13.25,9.6,train,ts,external,es);save(f,panels,'Fig2A_workflow',args.dpi)
    for m,l in zip(['R2','MAE','MAD'],['B','C','D']):
        f=plt.figure(figsize=(5.5,7));heatmap(f,t,m,l,1.25,1.5,3.75,4.45);save(f,panels,f'Fig2{l}_{m}_heatmap',args.dpi)
    for lineage,l,short in [('Fibroblast','E','Fibro'),('Endothelial','F','Endo')]:
        f=plt.figure(figsize=(7.3,7.8));final_fit_scatter(f,train,ts,lineage,l,1.0,1.5,5.55);save(f,panels,f'Fig2{l}_{short}_squareplot',args.dpi)
    pd.DataFrame(geometry).to_csv(out/'data/Fig2_square_geometry.csv',index=False)
    validation=dict(benchmark_rows=len(t),annotated_heatmap_cells=len(t)*3,
                    scatter_samples=len(train),legend_rows=len(ts),metrics_verified=['n','R2','MAE','MAD','PCC'],
                    heatmap_fold=4,model='ElasticNet',stable_features={'Fibroblast':750,'Endothelial':350},
                    final_scatter_scope='development_in_sample',workflow_external_samples=len(external),
                    MAD_definition_stable='median absolute prediction error',MAD_definition_heatmap='mean absolute age deviation',
                    square_panels=geometry,training_performed=False,
                    layout='portrait: A full-width top; B–D middle; E–F bottom',
                    figure_inches=[16.2,27.2],png_pixels=[16.2*args.dpi,27.2*args.dpi])
    (out/'data/validation_report.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(ts.to_string(index=False));print(es.to_string(index=False));print(json.dumps(validation,indent=2))


if __name__=='__main__':main()
