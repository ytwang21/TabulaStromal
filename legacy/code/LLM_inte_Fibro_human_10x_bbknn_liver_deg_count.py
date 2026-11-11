################################## Test the bbknn result based on lung Fibroblasts lineage
import pickle
import os
import re
import pandas as pd
import numpy as np
import sys
import scanpy as sc
# import scatlasvae  ##
import scipy.sparse as sp
import scipy
pickle_file = './data/scalex/fibro/human/Liver/Fibro_human_liver_FibDEG_counts.pickle'
with open(pickle_file,'rb') as f:
    adata_list_merge = pickle.load(f)

adata_summary_df = pd.read_csv('./data/scalex/fibro/human/Liver/Fibro_human_liver_FibDEG_counts.tsv',sep='\t')
adata_summary_df = adata_summary_df[(adata_summary_df['cell_num']>= 50) & (adata_summary_df['cell_num']<adata_summary_df['qc_cn']) ]
dataset_index = list(adata_summary_df.index)
adata_list_merge_new = [adata_list_merge[i] for i in dataset_index]
del adata_list_merge
adata = sc.concat(adata_list_merge_new,join='outer')   ### 972891 cells (后来经过筛选后，会少部分 cells)
import harmonypy as hm
del adata_list_merge_new
######## First do the scanpy pipeline
cells = adata.shape[0]
if cells <= 5000:
    n_pcs = 15
    resolution = 0.6
elif cells <= 10000:
    n_pcs = 20
    resolution = 0.8
elif cells <= 40000:
    n_pcs = 30
    resolution = 1
elif cells <= 80000:
    n_pcs = 40
    resolution =1
elif cells <= 150000:
    n_pcs = 50
    resolution = 1
else:
    n_pcs = 75
    resolution = 1
adata.var_names_make_unique()   #### gene name unique
adata.obs_names_make_unique()
sc.pl.highest_expr_genes(adata, n_top=20)  ### 先不找 HVGs
#### basic filtering, quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)  #### LLM 中是 min_cells是 500, 这步会filter掉  gene
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)  ### mito genes
sc.pp.normalize_total(adata, target_sum=1e4)   ### do normalization
sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)  ## 跑完这步会丢失很多gene
adata.raw = adata
# adata = adata[:, adata.var.highly_variable]  ## 只挑选 highly_variable也不对,其实跑完sc.pp.highly_variable_genes, adata的维度就变了
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")   ### PCA
sc.external.pp.bbknn(adata, batch_key="GSM")    ###  GSM 是batch
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs) #### Find neighborhoods, 用BBKNN做integration，不需要Find Neighbors
sc.tl.umap(adata)
sc.tl.leiden( adata, resolution=resolution,random_state=0,  # flavor="igraph",  # 加上这个参数会报错,
    n_iterations=2, directed=False)
adata.write('./data/scalex/fibro/human/Liver/Fibro_human_liver_FibDEG_counts.h5ad')


# adata.write('./small_file/Fibro_human_10x_TissueNoSplit_scatlasvae.h5ad')
