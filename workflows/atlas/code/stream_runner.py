#!/usr/bin/env python
"""Run actual STREAM seed AND elastic refinement; seeded MST alone is not a fit."""
import os
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ.setdefault(key,'4')
os.environ.setdefault('NUMBA_CACHE_DIR',str(HERE.parent/'tmp/numba_cache'))
os.environ.setdefault('MPLCONFIGDIR',str(HERE.parent/'tmp/matplotlib_cache'))

prepared=Path(sys.argv[1]);output=Path(sys.argv[2]);root_type=sys.argv[3]
rlib=Path(sys.argv[4]).resolve()
# Isolated project library first, existing global installation left untouched.
os.environ['R_LIBS_USER']=str(rlib)

import numpy as np
import pandas as pd
import networkx as nx
if not hasattr(nx,'from_numpy_matrix'):nx.from_numpy_matrix=nx.from_numpy_array
from anndata import AnnData
from scipy import sparse
import stream as st
from rpy2 import robjects
from rpy2.robjects import numpy2ri,pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.rinterface import ListSexpVector
from utils_data import dump_json

# Current rpy2 removed global activate(); use a scoped converter in this worker.
pandas2ri.activate=lambda:None
converter=robjects.default_converter+numpy2ri.converter+pandas2ri.converter
# STREAM expects nested R ListVectors with rx2(), not new rpy2 NamedLists.
converter.rpy2py.register(ListSexpVector,lambda obj:robjects.vectors.ListVector(obj))

robjects.r(f'.libPaths(c({str(rlib)!r}, .libPaths()))')
obs=pd.read_pickle(prepared/'metadata.pkl')
rep=np.load(prepared/'pca.npy')[:,:10]
xy=np.load(prepared/'umap.npy')
for representation,data in [('pca10',rep),('historical_umap',xy)]:
    dest=output/representation;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'method.json').exists():continue
    ad=AnnData(X=sparse.csr_matrix((len(obs),1),dtype=np.float32),obs=obs.copy(),dtype=np.float32)
    ad.obs['label']=ad.obs.general.astype(str)
    ad.obsm['X_dr']=data.astype(float)
    st.set_workdir(ad,str(dest))
    st.seed_elastic_principal_graph(ad,n_clusters=7,use_vis=False)
    with localconverter(converter):
        st.elastic_principal_graph(ad,epg_n_nodes=50,epg_alpha=.02,epg_mu=.1,epg_lambda=.02,
                                  epg_n_processes=1,save_fig=False,PlotGraph=False,ShowTimer=False)
    graph=ad.uns['epg'];flat=ad.uns['flat_tree']
    nodes=nx.get_node_attributes(flat,'pos');labels=nx.get_node_attributes(flat,'label')
    candidates=[n for n in flat if flat.degree(n)==1]
    root_center=np.median(data[obs.general.to_numpy()==root_type],axis=0)
    root_node=min(candidates,key=lambda n:np.linalg.norm(nodes[n]-root_center))
    alternative=sorted(candidates,key=lambda n:np.linalg.norm(nodes[n]-root_center))[min(1,len(candidates)-1)]
    root_label=labels[root_node];alt_label=labels[alternative]
    pt=ad.obs[root_label+'_pseudotime'].to_numpy(dtype=float)
    frame=pd.DataFrame({'cell_id':obs.index,'x':data[:,0],'y':data[:,1],'pseudotime':pt,
                        'general':obs.general.to_numpy(),'branch':ad.obs.branch_id_alias.astype(str).to_numpy(),
                        'alternative_root_pseudotime':ad.obs[alt_label+'_pseudotime'].to_numpy(dtype=float)})
    frame.to_csv(dest/'cells.csv',index=False)
    nodepos=nx.get_node_attributes(graph,'pos')
    pd.DataFrame([{'node':str(n),'x':v[0],'y':v[1],**{f'dim{i+1}':v[i] for i in range(len(v))}} for n,v in nodepos.items()]).to_csv(dest/'nodes.csv',index=False)
    pd.DataFrame([(str(a),str(b)) for a,b in graph.edges()],columns=['source','target']).to_csv(dest/'edges.csv',index=False)
    dump_json(dest/'flat_tree.json',nx.node_link_data(flat))
    dump_json(dest/'method.json',dict(method='STREAM seed + fitted ElPiGraph.R principal tree',stream_version=st.__version__,
              elpigraph_version=str(robjects.r('as.character(packageVersion("ElPiGraph.R"))')[0]),
              representation=representation,dimensions=data.shape[1],seed_clusters=7,elastic_nodes=50,
              root_annotation=root_type,root_graph_node=str(root_node),root_label=root_label,
              n_branchpoints=sum(flat.degree(n)>2 for n in flat),n_leaves=len(candidates),
              root_choice='Nearest principal-tree leaf to root-state median; direction is an assumption.',
              interpretive_limit='historical_umap reproduces legacy low-dimensional approach; it can inherit UMAP distortions and is compared against PCA10.'))
dump_json(output/'complete.json',dict(status='complete',representations=['pca10','historical_umap']))
print('STREAM_SUCCESS',len(obs),'meta-cells')
