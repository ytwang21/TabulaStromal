"""Monocle2, STREAM and CellRank comparison on the user's exact meta-cell atlases.

No PAGA. No invented cell-level velocity, fate probabilities or developmental truth.
All methods use the same meta-cell IDs; historical screenshots are reference only.
"""
import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time

# Initialize thread/cache locations before sklearn/anndata can import Numba.
import stromal_atlas_pipeline as _runtime

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from utils_data import read_frame, read_elem, dump_json, MARKERS, ORDER, COLORS, stratified_rows

HERE=Path(__file__).resolve().parent


def meta_config(config):
    return config['metacell_trajectory']


def corr(a,b):
    a=np.asarray(a);b=np.asarray(b);ok=np.isfinite(a)&np.isfinite(b)
    return float(spearmanr(a[ok],b[ok]).statistic) if ok.sum()>10 and np.ptp(a[ok])>0 and np.ptp(b[ok])>0 else np.nan


def unit(x):
    x=np.asarray(x,dtype=float);good=np.isfinite(x)
    if good.any():x[good]=(x[good]-x[good].min())/max(np.ptp(x[good]),1e-12)
    return x


def prepare(path,lineage,config,root):
    out=root/lineage/'trajectory';out.mkdir(parents=True,exist_ok=True)
    p=out/'prepared';p.mkdir(exist_ok=True)
    if (p/'complete.json').exists():
        spec=json.loads((p/'prepared.json').read_text())
        expected={f'HVG{k}' for k in meta_config(config)['monocle_hvg']}|{f'DEG{k}' for k in meta_config(config)['monocle_deg']}
        if set(spec['variants'])!=expected or spec['root_state']!=config['trajectory_root'][lineage]:
            raise ValueError('Prepared feature/root configuration changed; use a new run')
        audit=json.loads((out/'metacell_audit.json').read_text());stat=Path(path).stat()
        if stat.st_size!=audit['source_size'] or stat.st_mtime_ns!=audit['source_mtime_ns']:
            raise ValueError('Meta-cell input changed since preparation')
        apply_common_cohort(p,config,out)
        export_audit(out,root,lineage)
        return p
    cfg=meta_config(config)
    with h5py.File(path,'r') as f:
        cols=list(f['obs'].attrs['column-order'])
        keep=[k for k in ['Celltype_general','Celltype_sub','Celltype','tissue','sample','GSE','n_cells_in_meta','meta_cell_id','res2','leiden'] if k in cols]
        obs=read_frame(f['obs'],keep)
        rename={'Celltype_general':'general','Celltype_sub':'subtype','GSE':'dataset','Celltype':'celltype'}
        obs=obs.rename(columns=rename)
        for k in ['general','subtype','celltype','tissue','sample','dataset']:
            if k not in obs:raise ValueError(f'{lineage} meta-cell field missing: {k}')
            obs[k]=obs[k].astype(str)
        if set(obs.general)-set(ORDER[lineage]):raise ValueError('Unexpected broad annotation')
        var=read_frame(f['var']);genes=var.index.astype(str)
        if not obs.index.is_unique or not var.index.is_unique:raise ValueError('Nonunique meta-cell/gene IDs')
        if not isinstance(f['X'],h5py.Dataset):raise ValueError('Expected audited dense meta-cell mean expression')
        n,g=f['X'].shape
        means=np.zeros((len(ORDER[lineage]),g));sumsq=np.zeros_like(means);counts=np.zeros(len(means),int)
        codes=pd.Categorical(obs.general,categories=ORDER[lineage]).codes
        for start in range(0,n,512):
            stop=min(n,start+512);x=f['X'][start:stop].astype(np.float32)
            if not np.isfinite(x).all() or np.min(x)<0:raise ValueError('Invalid mean expression')
            for j in range(len(counts)):
                a=x[codes[start:stop]==j];counts[j]+=len(a)
                means[j]+=a.sum(axis=0,dtype=np.float64);sumsq[j]+=np.square(a,dtype=np.float64).sum(axis=0)
        total_mean=means.sum(axis=0)/n
        state_mean=means/np.maximum(counts[:,None],1)
        between=np.sum(counts[:,None]*(state_mean-total_mean)**2,axis=0)/max(len(counts)-1,1)
        within=np.maximum((sumsq-means**2/np.maximum(counts[:,None],1)).sum(axis=0)/max(n-len(counts),1),1e-10)
        fscore=between/within
        variance=np.maximum(sumsq.sum(axis=0)/n-total_mean**2,0)
        # Label-informed rankings are explicitly distinguished from unsupervised HVGs.
        hvgrank=np.argsort(-variance)
        degrank=np.argsort(-fscore)
        eligible=(total_mean>1e-5)&(~genes.str.match(r'^(MT-|RPL|RPS)'))
        hvgrank=[i for i in hvgrank if eligible[i]];degrank=[i for i in degrank if eligible[i]]
        variants={f'HVG{k}':genes[hvgrank[:k]].tolist() for k in cfg['monocle_hvg']}
        variants.update({f'DEG{k}':genes[degrank[:k]].tolist() for k in cfg['monocle_deg']})
        selected=sorted({g for gs in variants.values() for g in gs})
        indices=genes.get_indexer(selected)
        dest=np.memmap(p/'expression.f32',mode='w+',dtype='<f4',shape=(n,len(selected)))
        for start in range(0,n,512):dest[start:min(n,start+512)]=f['X'][start:min(n,start+512)][:,indices]
        dest.flush();del dest
        pd.DataFrame({'gene':selected}).to_csv(p/'genes.csv',index=False)
        pd.DataFrame({'gene':genes,'variance':variance,'subtype_association_F':fscore,'mean':total_mean,'eligible':eligible}).to_csv(p/'feature_ranking.csv',index=False)
        marker_genes=[gene for state in ORDER[lineage] for gene in sorted(MARKERS[lineage][state]) if gene in set(genes)]
        mcols=genes.get_indexer(marker_genes)
        mx=np.empty((n,len(mcols)),np.float32)
        for start in range(0,n,512):mx[start:min(n,start+512)]=f['X'][start:min(n,start+512)][:,mcols]
        np.save(p/'marker_expression.npy',mx);dump_json(p/'marker_genes.json',marker_genes)
        pca=f['obsm/X_pca'][:,:30].astype(np.float32);umap=f['obsm/X_umap'][:,:2].astype(np.float32)
        np.save(p/'pca.npy',pca);np.save(p/'umap.npy',umap)
        for k in ['connectivities','distances']:
            if f'obsp/{k}' not in f:raise ValueError('Audited meta-cell neighbor graph missing')
            sparse.save_npz(p/f'{k}.npz',read_elem(f[f'obsp/{k}']))
        obs.index.name='cell_id';obs.to_csv(p/'metadata.csv');obs.to_pickle(p/'metadata.pkl')
        size=obs.n_cells_in_meta.to_numpy(dtype=float)
        audit=dict(input=str(Path(path).resolve()),n_meta_cells=n,n_genes=g,
                   represented_cell_sum=int(size.sum()),aggregation='Mean of existing expression, per historical create_meta_cells',
                   n_cells_in_meta_quantiles=dict(zip(['min','q25','median','q75','q95','q99','max'],np.quantile(size,[0,.25,.5,.75,.95,.99,1]).tolist())),
                   n_large_meta_cells=int((size>60).sum()),
                   limitation='Historical leftover groups may contain far more than 30 cells; each meta-cell receives equal modeling weight; size-stratified diagnostics supplied.',
                   fields=cols,states=obs.general.value_counts().to_dict(),
                   root_state=config['trajectory_root'][lineage],source_size=Path(path).stat().st_size,
                   source_mtime_ns=Path(path).stat().st_mtime_ns)
        dump_json(out/'metacell_audit.json',audit)
        obs.loc[size>60].to_csv(out/'large_metacells_over60.csv')
        spec=dict(n_cells=n,n_genes=len(selected),variants=variants,root_state=config['trajectory_root'][lineage],
                  expression_scale='averaged log expression, not counts',feature_selection='HVG: empirical meta-cell variance; DEG: descriptive between-state F ranking (label-informed; not donor-independent DE testing)')
        dump_json(p/'prepared.json',spec)
    dump_json(p/'complete.json',dict(status='complete',n_meta_cells=len(obs)))
    apply_common_cohort(p,config,out)
    export_audit(out,root,lineage)
    logging.info('%s prepared %s meta-cells; size median=%s max=%s',lineage,n,np.median(size),max(size))
    return p


def export_audit(out,root,lineage):
    audit=json.loads((out/'metacell_audit.json').read_text())
    with h5py.File(audit['input'],'r') as f:
        audit.update(has_raw='raw' in f,layers=list(f['layers']),obsm=list(f['obsm']),
                     var_columns=list(f['var'].attrs['column-order']),expression_dtype=str(f['X'].dtype),
                     expression_scale='Average of preexisting expression, per source meta-cell code; raw is not assumed to contain counts')
    dump_json(out/'metacell_audit.json',audit)
    shared=root/'shared';shared.mkdir(exist_ok=True,parents=True)
    dump_json(shared/f'data_audit_metacell_{lineage}.json',audit)
    pd.DataFrame([{'item':k,'value':json.dumps(v,ensure_ascii=False)} for k,v in audit.items()]).to_csv(shared/f'data_audit_metacell_{lineage}.csv',index=False)


def apply_common_cohort(p,config,out):
    """Same recorded meta-cell cohort for all three tools; no projection/imputation."""
    cap=meta_config(config).get('max_fit_metacells',6000)
    if (p/'cohort.json').exists():
        previous=json.loads((p/'cohort.json').read_text())
        if previous['configured_cap']!=cap or previous['seed']!=config['seed']:raise ValueError('Meta-cell cohort changed; start a new run')
        return
    obs=pd.read_pickle(p/'metadata.pkl');n=len(obs)
    rows=stratified_rows(obs,min(cap,n) if cap else n,config['seed'],('subtype','tissue'))
    spec=json.loads((p/'prepared.json').read_text())
    # Preserve original full-source identity/embedding audits independently of fit data.
    for name in ['metadata.pkl','metadata.csv','pca.npy','umap.npy']:
        import shutil
        shutil.copy2(p/name,p/('source_'+name))
    pd.DataFrame({'source_row':rows,'cell_id':obs.index[rows]}).to_csv(p/'cohort_cells.csv',index=False)
    if len(rows)!=n:
        old=np.memmap(p/'expression.f32',mode='r',dtype='<f4',shape=(n,spec['n_genes']))
        dest=np.memmap(p/'expression_subset.f32',mode='w+',dtype='<f4',shape=(len(rows),spec['n_genes']))
        dest[:]=old[rows];dest.flush();del old,dest
        (p/'expression_subset.f32').replace(p/'expression.f32')
        for name in ['pca.npy','umap.npy','marker_expression.npy']:
            a=np.load(p/name);np.save(p/name,a[rows])
        sub=obs.iloc[rows].copy();sub.to_pickle(p/'metadata.pkl');sub.to_csv(p/'metadata.csv')
        import scanpy as sc
        from anndata import AnnData
        ad=AnnData(sparse.csr_matrix((len(rows),1),dtype=np.float32),dtype=np.float32)
        ad.obsm['X_pca']=np.load(p/'pca.npy')
        sc.pp.neighbors(ad,n_neighbors=30,use_rep='X_pca',random_state=config['seed'])
        for key in ['connectivities','distances']:sparse.save_npz(p/f'{key}.npz',ad.obsp[key])
        spec['n_cells']=len(rows);dump_json(p/'prepared.json',spec)
    cohort=dict(configured_cap=cap,n_source_metacells=n,n_fit_metacells=len(rows),seed=config['seed'],
                sampling='Reserve one per subtype×tissue stratum, fill proportionally from remaining meta-cells; same IDs for Monocle, STREAM, CellRank',
                graph='kNN 30 in existing meta-cell PCA for common sampled cohort; no induced sparse BBKNN subgraph' if len(rows)!=n else 'Full original meta-cell neighbor graph retained without inducing a sampled graph',
                original_metacell_size='Retained, not reaggregated; each fit point is still an original user-provided meta-cell')
    dump_json(p/'cohort.json',cohort)
    audit=json.loads((out/'metacell_audit.json').read_text());audit['modeling_cohort']=cohort;dump_json(out/'metacell_audit.json',audit)


def run_monocle(prepared,lineage,config,out):
    cfg=meta_config(config);spec=json.loads((prepared/'prepared.json').read_text())
    def fit(variant):
        p=out/'monocle'/variant;p.mkdir(parents=True,exist_ok=True)
        if (p/'method.json').exists():return
        import shutil
        runner=p/'monocle_runner_used.R';shutil.copy2(HERE/'monocle_runner.R',runner)
        command=[cfg['rscript'],str(runner),str(prepared),variant,str(p)]
        env=os.environ.copy()
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
            env[key]=str(cfg.get('r_threads',2))
        logging.info('%s Monocle %s',lineage,variant)
        with (p/'run.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=cfg.get('monocle_parallel_jobs',3)) as pool:
        futures=[pool.submit(fit,v) for v in spec['variants']]
        errors=[]
        for future in futures:
            try:future.result()
            except Exception as exc:errors.append(exc)
        if errors:raise errors[0]


def graph_adata(prepared):
    from anndata import AnnData
    obs=pd.read_pickle(prepared/'metadata.pkl')
    ad=AnnData(X=sparse.csr_matrix((len(obs),1),dtype=np.float32),obs=obs.copy(),dtype=np.float32)
    ad.obsm['X_pca']=np.load(prepared/'pca.npy');ad.obsm['X_umap']=np.load(prepared/'umap.npy')
    ad.obsp['connectivities']=sparse.load_npz(prepared/'connectivities.npz').tocsr()
    ad.obsp['distances']=sparse.load_npz(prepared/'distances.npz').tocsr()
    ad.uns['neighbors']=dict(connectivities_key='connectivities',distances_key='distances',params=dict(use_rep='X_pca',n_neighbors=30,method='umap'))
    return ad


def run_cellrank(prepared,lineage,config,out):
    import scanpy as sc
    import cellrank as cr
    from scipy.sparse.csgraph import connected_components
    p=out/'cellrank';p.mkdir(exist_ok=True)
    if (p/'method.json').exists():return
    ad=graph_adata(prepared);root_type=config['trajectory_root'][lineage]
    nc,component=connected_components(ad.obsp['connectivities'],directed=False)
    cand=np.flatnonzero(ad.obs.general.to_numpy()==root_type)
    root_comp=np.bincount(component[cand]).argmax()
    ids=np.flatnonzero(component==root_comp);ad=ad[ids].copy()
    rep=ad.obsm['X_pca'];cand=np.flatnonzero(ad.obs.general.to_numpy()==root_type)
    distances=np.linalg.norm(rep[cand]-np.median(rep[cand],axis=0),axis=1)
    roots=cand[np.argsort(distances)]
    root=int(roots[0]);alternative=int(roots[int(.75*(len(roots)-1))])
    sc.tl.diffmap(ad,n_comps=16,random_state=config['seed'])
    ad.uns['iroot']=root;sc.tl.dpt(ad,n_dcs=10,n_branchings=0)
    pt=ad.obs.dpt_pseudotime.to_numpy().copy()
    ad.uns['iroot']=alternative;sc.tl.dpt(ad,n_dcs=10,n_branchings=0)
    alternate=ad.obs.dpt_pseudotime.to_numpy().copy()
    ad.obs['oriented_pseudotime']=pt
    pk=cr.kernels.PseudotimeKernel(ad,time_key='oriented_pseudotime')
    pk.compute_transition_matrix(threshold_scheme='soft',n_jobs=1,show_progress_bar=False)
    ck=cr.kernels.ConnectivityKernel(ad).compute_transition_matrix()
    kernel=.8*pk+.2*ck;kernel.compute_transition_matrix()
    t=kernel.transition_matrix.tocsr();sparse.save_npz(p/'transition_matrix.npz',t)
    xy=ad.obsm['X_umap'];flow=t@xy-xy
    frame=pd.DataFrame({'cell_id':ad.obs_names,'x':xy[:,0],'y':xy[:,1],'pseudotime':unit(pt),'general':ad.obs.general.to_numpy(),
                        'component':root_comp,'flow_x':flow[:,0],'flow_y':flow[:,1]})
    # Without PETSc/SLEPc, CellRank's full Schur backend densifies n×n.
    # Fit GPCCA to an explicitly recorded Markov lumping, then solve fate
    # probabilities on the original sparse meta-cell transition matrix.
    from sklearn.cluster import MiniBatchKMeans
    nland=min(meta_config(config).get('cellrank_landmarks',200),len(ad))
    clustering=MiniBatchKMeans(n_clusters=nland,random_state=config['seed'],n_init=3,batch_size=1024).fit(rep)
    assignment=clustering.labels_
    h=sparse.csr_matrix((np.ones(len(ad)),(np.arange(len(ad)),assignment)),shape=(len(ad),nland))
    lumped=(sparse.diags(1/np.maximum(np.asarray(h.sum(axis=0)).ravel(),1))@h.T@t@h).tocsr()
    from anndata import AnnData
    landmark_obs=pd.DataFrame({'general':[ad.obs.general.iloc[np.flatnonzero(assignment==i)].mode().iloc[0] for i in range(nland)]},index=[f'landmark_{i}' for i in range(nland)])
    landmark_obs['general']=landmark_obs.general.astype('category')
    landmark_ad=AnnData(sparse.csr_matrix((nland,1),dtype=np.float32),obs=landmark_obs,dtype=np.float32)
    landmark_kernel=cr.kernels.PrecomputedKernel(lumped,adata=landmark_ad)
    gpcca=cr.estimators.GPCCA(landmark_kernel)
    nstates=len(ORDER[lineage])
    gpcca.compute_schur(n_components=max(12,nstates+3))
    gpcca.compute_macrostates(n_states=nstates,cluster_key='general',n_cells=1)
    macro_names=list(gpcca.macrostates_memberships.names)
    memberships=np.asarray(gpcca.macrostates_memberships)[assignment]
    frame['branch']=np.array(macro_names)[memberships.argmax(axis=1)]
    gpcca.predict_terminal_states(method='top_n',n_states=min(3,nstates-1),n_cells=1)
    terminal={}
    for state in gpcca.terminal_states.cat.categories:
        landmarks=np.flatnonzero(gpcca.terminal_states.astype(str).to_numpy()==state)
        cells=np.flatnonzero(np.isin(assignment,landmarks))
        col=macro_names.index(state)
        chosen=cells[np.argsort(-memberships[cells,col])[:min(30,len(cells))]]
        terminal[str(state)]=ad.obs_names[chosen].tolist()
    full_gpcca=cr.estimators.GPCCA(kernel)
    full_gpcca.set_terminal_states(terminal)
    # CellRank 2.0.7 uses SciPy's newer rtol spelling; SciPy 1.11 uses tol.
    import inspect
    import cellrank._utils._linear_solver as solver_module
    old_solver=solver_module._AVAIL_ITER_SOLVERS['gmres']
    if 'rtol' not in inspect.signature(old_solver).parameters:
        def compatible_gmres(a,b,rtol=1e-6,**kwargs):
            return old_solver(a,b,tol=rtol,**kwargs)
        solver_module._AVAIL_ITER_SOLVERS['gmres']=compatible_gmres
    full_gpcca.compute_fate_probabilities(solver='gmres',use_petsc=False,n_jobs=1,show_progress_bar=False,tol=1e-6)
    fate=np.asarray(full_gpcca.fate_probabilities)
    fate_names=list(full_gpcca.fate_probabilities.names)
    pd.DataFrame(fate,index=ad.obs_names,columns=fate_names).to_csv(p/'fate_probabilities.csv.gz',index_label='cell_id')
    coarse=gpcca.coarse_T
    pd.DataFrame(np.asarray(coarse),index=gpcca.macrostates_memberships.names,columns=gpcca.macrostates_memberships.names).to_csv(p/'macrostate_transitions.csv')
    sparse.save_npz(p/'landmark_transition_matrix.npz',lumped)
    pd.DataFrame({'cell_id':ad.obs_names,'landmark':assignment,'macrostate':frame['branch']}).to_csv(p/'landmark_mapping.csv',index=False)
    pd.DataFrame(memberships,index=ad.obs_names,columns=macro_names).to_csv(p/'macrostate_memberships.csv.gz',index_label='cell_id')
    dump_json(p/'terminal_core_metacells.json',terminal)
    frame.to_csv(p/'cells.csv',index=False)
    # Retain all excluded components with NA in a separate full-ID table.
    all_obs=pd.read_pickle(prepared/'metadata.pkl')
    full=all_obs[['general']].join(frame.set_index('cell_id').drop(columns='general'),how='left')
    full.to_csv(p/'all_metacells.csv.gz',index_label='cell_id')
    dump_json(p/'method.json',dict(method='CellRank PseudotimeKernel(0.8) + ConnectivityKernel(0.2) + GPCCA',version=cr.__version__,
              root_cell=str(ad.obs_names[root]),alternative_root_cell=str(ad.obs_names[alternative]),root_annotation=root_type,
              root_spearman=corr(pt,alternate),n_cells=len(ad),total_metacells=len(all_obs),n_components=int(nc),
              terminal_state_rule='top 3 stable coarse states, not manually forced to expected subtype',fate_names=fate_names,
              gpcca_approximation=f'{nland}-landmark Markov lumping H^T T H weighted by cluster sizes; GPCCA on lumped chain, fate probabilities solved on full sparse common-cohort T',
              directional_prior='DPT orientation is an assumed root; this is not RNA velocity or independent evidence of direction',
              flow_definition='Expected displacement T @ X_umap - X_umap from actual CellRank transition matrix'))
    logging.info('%s CellRank finished',lineage)


def run_stream(prepared,lineage,config,out):
    # Separate process isolates the legacy Python/R bridge from CellRank.
    p=out/'stream';p.mkdir(exist_ok=True)
    if (p/'complete.json').exists():return
    command=[sys.executable,str(HERE/'stream_runner.py'),str(prepared),str(p),config['trajectory_root'][lineage],str(HERE.parent/'tmp/Rlib')]
    with (p/'run.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)


def trajectory_stage(path,lineage,config,root):
    prepared=prepare(config['metacell_inputs'][lineage],lineage,config,root)
    out=root/lineage/'trajectory'
    for method in ['monocle','stream','cellrank']:
        {'monocle':run_monocle,'stream':run_stream,'cellrank':run_cellrank}[method](prepared,lineage,config,out)
    compare_and_plot(lineage,config,root)


def eta2(x,labels):
    frame=pd.DataFrame({'x':x,'g':labels}).dropna()
    if frame.x.var()==0:return 0.
    t=frame.groupby('g').x.agg(['mean','size'])
    return float(np.sum(t['size']*(t['mean']-frame.x.mean())**2)/np.sum((frame.x-frame.x.mean())**2))


def model_paths(out):
    return sorted([p.parent for p in out.glob('monocle/*/method.json')])+sorted([p.parent for p in out.glob('stream/*/method.json')])+([out/'cellrank'] if (out/'cellrank/method.json').exists() else [])


def load_model(path,obs):
    frame=pd.read_csv(path/'cells.csv',dtype={'cell_id':str,'branch':str}).set_index('cell_id')
    if not frame.index.is_unique:raise ValueError('Duplicate fitted meta-cell IDs')
    frame['pseudotime']=unit(frame.pseudotime)
    return frame.reindex(obs.index)


def draw_model(ax,path,obs,lineage,color='general',show_title=True):
    f=load_model(path,obs);info=json.loads((path/'method.json').read_text())
    good=np.isfinite(f.x)&np.isfinite(f.y)
    f=f.loc[good]
    if color=='general':
        sc=ax.scatter(f.x,f.y,c=[COLORS[k] for k in obs.loc[f.index,'general']],s=1.7,lw=0,alpha=.65,rasterized=True)
        from matplotlib.lines import Line2D
        handles=[Line2D([],[],ls='',marker='o',markersize=3,color=COLORS[k],label=k) for k in ORDER[lineage]]
        ax.legend(handles=handles,loc='upper left',frameon=False,fontsize=5.5,handletextpad=.25,labelspacing=.25)
    else:
        sc=ax.scatter(f.x,f.y,c=f.pseudotime,cmap='viridis',vmin=0,vmax=1,s=1.7,lw=0,rasterized=True)
    if (path/'nodes.csv').exists():
        nodes=pd.read_csv(path/'nodes.csv',dtype={'node':str}).set_index('node')
        edges=pd.read_csv(path/'edges.csv',dtype=str)
        for _,e in edges.iterrows():
            a,b=e['source'],e['target']
            if a in nodes.index and b in nodes.index:
                x,y=nodes.loc[a],nodes.loc[b]
                ax.plot([x.x,y.x],[x.y,y.y],c='#333333',lw=.75,alpha=.8)
    elif 'flow_x' in f:
        # Display T-derived expected displacements, not a fitted pseudotime gradient.
        x=f[['x','y']].to_numpy();velocity=f[['flow_x','flow_y']].to_numpy()
        from scipy.ndimage import gaussian_filter
        lo=x.min(axis=0);hi=x.max(axis=0)
        gx=np.linspace(lo[0],hi[0],35);gy=np.linspace(lo[1],hi[1],35)
        grid=np.array(np.meshgrid(gx,gy)).reshape(2,-1).T
        nn=NearestNeighbors(n_neighbors=min(30,len(x))).fit(x)
        d,idx=nn.kneighbors(grid);scale=max(np.median(d[:,-1]),1e-9)
        w=np.exp(-np.square(d/scale));vf=np.sum(velocity[idx]*w[:,:,None],axis=1)/np.maximum(w.sum(axis=1)[:,None],1e-9)
        allowed=d[:,0]<scale
        # Actual transition displacements are kernel-smoothed on the grid;
        # mask unsupported regions so streamlines do not cross empty islands.
        vx=np.ma.masked_where(~allowed.reshape(len(gy),len(gx)),vf[:,0].reshape(len(gy),len(gx)))
        vy=np.ma.masked_where(~allowed.reshape(len(gy),len(gx)),vf[:,1].reshape(len(gy),len(gx)))
        ax.streamplot(gx,gy,vx,vy,color='#303030',linewidth=.7,density=.85,arrowsize=.8,minlength=.15,maxlength=1.2)
    if len(f):
        root=f.loc[f.pseudotime.idxmin()]
        ax.scatter(root.x,root.y,marker='*',s=75,c='white',edgecolor='#222222',lw=.7,zorder=10)
    ax.set_aspect('equal',adjustable='datalim');ax.set(xticks=[],yticks=[])
    for sp in ax.spines.values():sp.set_visible(False)
    family=path.parent.name if path.name!='cellrank' else 'cellrank'
    name='Monocle 2 · '+path.name if family=='monocle' else 'STREAM · '+path.name if family=='stream' else 'CellRank · transition flow'
    if show_title:ax.set_title(f'{lineage.capitalize()} | {name}',loc='left',fontsize=9)
    ax.text(.01,-.05,f'{len(f):,} meta-cells; star = assumed root',transform=ax.transAxes,fontsize=6)
    return sc


def compare_and_plot(lineage,config,root):
    from utils_plot import save
    import matplotlib.pyplot as plt
    out=root/lineage/'trajectory';prepared=out/'prepared';obs=pd.read_pickle(prepared/'metadata.pkl')
    paths=model_paths(out)
    families={'monocle':[], 'stream':[], 'cellrank':[]}
    frames={};records=[]
    rep=np.load(prepared/'pca.npy')
    nn=NearestNeighbors(n_neighbors=16).fit(rep)
    neigh=nn.kneighbors(rep,return_distance=False)[:,1:]
    for p in paths:
        fam='cellrank' if p.name=='cellrank' else p.parent.name
        key=fam if fam=='cellrank' else f'{fam}_{p.name}'
        families[fam].append(key);frame=load_model(p,obs);frames[key]=frame
        pt=frame.pseudotime.to_numpy();finite=np.isfinite(pt)
        local=np.nanmean(np.square(pt[:,None]-pt[neigh]))
        variance=np.nanvar(pt)
        branch=pd.crosstab(frame.branch,obs.general)
        purity=float(branch.max(axis=1).sum()/max(branch.to_numpy().sum(),1))
        info=json.loads((p/'method.json').read_text())
        nbranch=info.get('n_branchpoints')
        if (p/'edges.csv').exists():
            ed=pd.read_csv(p/'edges.csv',dtype=str)
            degree=pd.concat([ed.source,ed.target]).value_counts();nbranch=int((degree>2).sum())
        record=dict(model=key,family=fam,path=str(p.relative_to(out)),n_metacells=len(obs),finite_fraction=float(finite.mean()),
                    local_roughness=float(local/max(variance,1e-12)),branch_state_purity=purity,
                    n_branchpoints=nbranch,tissue_eta2=eta2(pt,obs.tissue),study_eta2=eta2(pt,obs.dataset),
                    size_spearman=corr(pt,np.log1p(obs.n_cells_in_meta)),
                    label_informed_features=fam=='monocle' and p.name.startswith('DEG'))
        records.append(record)
        fig,axes=plt.subplots(1,2,figsize=(11,4.2))
        draw_model(axes[0],p,obs,lineage,'general')
        im=draw_model(axes[1],p,obs,lineage,'pseudotime')
        fig.colorbar(im,ax=axes[1],fraction=.025,pad=.02,label='Root-relative pseudotime (0–1)')
        save(fig,p/'trajectory',config['dpi'])
    if any(not v for v in families.values()):raise ValueError(f'Three-method comparison incomplete: {families}')
    # Equal method-family weighting avoids six Monocle fits dominating consensus.
    concord=pd.DataFrame(index=frames,columns=frames,dtype=float)
    for a in frames:
        for b in frames:concord.loc[a,b]=corr(frames[a].pseudotime,frames[b].pseudotime)
    concord.to_csv(out/'pseudotime_concordance.csv')
    for rec in records:
        others=[np.nanmedian([concord.loc[rec['model'],k] for k in keys]) for fam,keys in families.items() if fam!=rec['family']]
        rec['cross_method_concordance']=float(np.nanmean(others))
        # Presentation ranking is explicitly heuristic, not a correctness metric.
        rec['descriptive_score']=float(rec['finite_fraction']*(.5/(1+rec['local_roughness'])+.5*(rec['cross_method_concordance']+1)/2))
    metrics=pd.DataFrame(records).sort_values('descriptive_score',ascending=False)
    metrics.to_csv(out/'three_method_comparison.csv',index=False)
    best_by_family={fam:metrics[metrics.family==fam].iloc[0].to_dict() for fam in families}
    winner=metrics.iloc[0]
    selected=dict(reference=winner.model,main_panel='metacell_three_method',method=winner.family,
                  model=winner.model,path=winner.path,n_meta_cells=len(obs),best_by_family=best_by_family,
                  rationale='Highest descriptive score combining local smoothness and equal-family root-relative pseudotime concordance. Not a measure of true developmental correctness. DEG-informed feature selection is explicitly flagged; no forced endpoint or subtype arrows.')
    dump_json(out/'trajectory_selection.json',selected)
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    settings=[f'HVG{k}' for k in meta_config(config)['monocle_hvg']]+[f'DEG{k}' for k in meta_config(config)['monocle_deg']]
    mono=sorted([p for p in paths if p.parent.name=='monocle'],key=lambda p:settings.index(p.name))
    for ax,p in zip(axes.ravel(),mono):draw_model(ax,p,obs,lineage)
    for ax in axes.ravel()[len(mono):]:ax.set_axis_off()
    save(fig,out/f'{lineage}_monocle_feature_comparison',config['dpi'])
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for j,fam in enumerate(['monocle','stream','cellrank']):
        p=out/best_by_family[fam]['path']
        draw_model(axes[0,j],p,obs,lineage)
        im=draw_model(axes[1,j],p,obs,lineage,'pseudotime')
        fig.colorbar(im,ax=axes[1,j],fraction=.025,pad=.02)
    save(fig,out/f'{lineage}_three_method_comparison',config['dpi'])
    mx=np.load(prepared/'marker_expression.npy');genes=json.loads((prepared/'marker_genes.json').read_text())
    trends=[];associations=[]
    for rec in records:
        pt=frames[rec['model']].pseudotime.to_numpy();valid=np.flatnonzero(np.isfinite(pt))
        bins=pd.qcut(pt[valid],20,duplicates='drop',labels=False)
        for j,gene in enumerate(genes):
            associations.append(dict(model=rec['model'],gene=gene,spearman=corr(pt,mx[:,j])))
        for b in np.unique(bins):
            ids=valid[bins==b]
            for j,gene in enumerate(genes):trends.append(dict(model=rec['model'],bin=int(b),gene=gene,pseudotime=float(np.median(pt[ids])),mean_expression=float(mx[ids,j].mean()),n_metacells=len(ids)))
    trend=pd.DataFrame(trends);trend.to_csv(out/'marker_pseudotime_trends.csv',index=False)
    pd.DataFrame(associations).to_csv(out/'marker_pseudotime_associations.csv',index=False)
    fig,axes=plt.subplots(1,len(ORDER[lineage]),figsize=(14,3),layout='constrained')
    for ax,state in zip(np.atleast_1d(axes),ORDER[lineage]):
        for gene in sorted(MARKERS[lineage][state]):
            t=trend[(trend.model==winner.model)&(trend.gene==gene)]
            if len(t):ax.plot(t.pseudotime,t.mean_expression,label=gene,lw=1)
        ax.set(title=state,xlabel='Pseudotime',ylabel='Mean log expression')
        if ax.lines:ax.legend(frameon=False,fontsize=5)
    save(fig,out/f'{lineage}_marker_trends',config['dpi'])
    sizes=obs.n_cells_in_meta.to_numpy()
    diagnostic=[]
    for rec in records:
        pt=frames[rec['model']].pseudotime
        for label,mask in [('15–60 cells',sizes<=60),('>60 cells',sizes>60)]:
            diagnostic.append(dict(model=rec['model'],size_group=label,n_metacells=int(mask.sum()),median_pseudotime=float(np.nanmedian(pt[mask])),pseudotime_mean=float(np.nanmean(pt[mask]))))
    pd.DataFrame(diagnostic).to_csv(out/'metacell_size_diagnostic.csv',index=False)
    full=frames[winner.model].join(obs[['sample','dataset','tissue','n_cells_in_meta']])
    full.to_csv(out/'selected_pseudotime.csv.gz',index_label='cell_id')
    report=f'''# {lineage}: Monocle 2 / STREAM / CellRank on meta-cells

Input: {config['metacell_inputs'][lineage]}. The original source has {json.loads((out/'metacell_audit.json').read_text())['n_meta_cells']:,} meta-cells; the same recorded representative {len(obs):,} meta-cells are used for all three tools. No single cells are reaggregated or substituted. Set max_fit_metacells=0 in a NEW run to fit the entire original meta-cell object; Monocle2's dense ordering step is substantially more expensive. No PAGA is run or shown.

Selected display: {winner.model}; score {winner.descriptive_score:.3f}. This is a descriptive presentation ranking, not evidence of the true differentiation tree. See three_method_comparison.csv and pseudotime_concordance.csv for all actual metrics.

Monocle: six DDRTree feature settings (HVG 750/1000/1500; subtype-associated DEG-informed 50/150/350). These inputs are averaged log expression: uninormal family, no count likelihood and no second log normalization. DEG-informed rankings use annotation labels and therefore cannot independently validate those labels or a desired Y shape.

STREAM: actual seed AND refined elastic-principal-tree fits in PCA10 and existing UMAP. Historical Notebook seeded seven centers and some later refinement calls were commented out; a seed MST alone is not an optimized trajectory. PCA10 is compared against the historical 2D approach rather than silently assuming either is correct.

CellRank: genuine PseudotimeKernel/ConnectivityKernel transition matrix and solved absorption/fate probabilities. GPCCA is applied to a recorded {meta_config(config).get('cellrank_landmarks',200)}-landmark Markov lumping to avoid the unavailable PETSc backend densifying the full matrix; sparse fate probabilities are then solved on the original common-cohort transition matrix. This is an explicitly recorded approximation to macrostate discovery, not the full-matrix Schur decomposition. Expected transition displacements T@UMAP−UMAP are plotted. Root-oriented DPT supplies a directional prior; no spliced/unspliced data exist, so arrows are not RNA velocity. Initial/terminal states and root assumptions are recorded. Historical one-hot macrostate labels are not used as probabilities, and no gradient-based or UMAP-distance fallback is presented as a CellRank inference.

Root prior: {config['trajectory_root'][lineage]}. Snapshot atlases without lineage tracing, actual time or velocity do not identify a uniquely correct developmental direction. Root orientation, subtype purity, visible branches and marker trends remain descriptive. Tissue, study, disease and sample composition may explain structure; unadjusted eta² is biased upward when many groups exist. MCs and LEC are retained in Endo and may be distinct populations rather than vascular developmental intermediates.

Meta-cell sizes: median {np.median(sizes):g}, maximum {sizes.max():g}, >60 cells: {(sizes>60).sum()}. The historical aggregation function can combine all leftover cells of a group, creating a long tail. Each meta-cell is modeled with equal weight; high-size associations are exported. Original_cells records local group indices in the historical code, so it must not be interpreted as globally unique lineage tracing. Inputs are unchanged.

Outputs: native topology/pseudotime figures for every fit; Monocle six-panel grid; three-method comparison; real CellRank transition/fate tables; all-ID pseudotime, marker trends, size/tissue/study diagnostics and exact feature lists.

Sources: https://cole-trapnell-lab.github.io/monocle-release/docs/ ; https://github.com/pinellolab/STREAM ; https://cellrank.readthedocs.io/en/latest/api/kernels.html
'''
    (out/'trajectory_report.md').write_text(report)
    dump_json(out/'comparison_complete.json',dict(status='complete',n_models=len(records),selected=winner.model))


def draw_selected(ax,lineage,config,root):
    out=root/lineage/'trajectory';obs=pd.read_pickle(out/'prepared/metadata.pkl')
    choice=json.loads((out/'trajectory_selection.json').read_text())
    return draw_model(ax,out/choice['path'],obs,lineage)


def main():
    import stromal_atlas_pipeline
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=HERE/'config.json')
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--lineage',choices=['fibro','endo'],required=True)
    parser.add_argument('--method',choices=['prepare','monocle','stream','cellrank','compare'],required=True)
    args=parser.parse_args()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    cfg=json.loads(args.config.read_text())
    cfg['metacell_inputs']={k:str((args.config.resolve().parent/v).resolve()) for k,v in cfg['metacell_inputs'].items()}
    p=prepare(cfg['metacell_inputs'][args.lineage],args.lineage,cfg,args.run_dir)
    if args.method in ['monocle','stream','cellrank']:
        {'monocle':run_monocle,'stream':run_stream,'cellrank':run_cellrank}[args.method](p,args.lineage,cfg,args.run_dir/args.lineage/'trajectory')
    elif args.method=='compare':compare_and_plot(args.lineage,cfg,args.run_dir)


if __name__=='__main__':main()
