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
pickle_file_1 = './data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge.pickle'
with open(pickle_file_1,'rb') as f:
    adata_list_merge1 = pickle.load(f)

adata_summary_df1 = pd.read_csv('./data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge.tsv',sep='\t')
adata_summary_df1 = adata_summary_df1.reset_index()
adata_summary_df1 = adata_summary_df1[(adata_summary_df1['cell_num']>= 50) & (adata_summary_df1['cell_num']<adata_summary_df1['qc_cn']) ]
dataset_index1 = list(adata_summary_df1.index)
adata_list_merge_new1 = [adata_list_merge1[i] for i in dataset_index1]
del adata_list_merge1

pickle_file_2 = './data/scalex/fibro/human/Third_Fibro_human_FibDEG_count_notAge.pickle'
with open(pickle_file_2,'rb') as f:
    adata_list_merge2 = pickle.load(f)

adata_summary_df2 = pd.read_csv('./data/scalex/fibro/human/Third_Fibro_human_FibDEG_count_notAge.tsv',sep='\t')
adata_summary_df2 = adata_summary_df2[(adata_summary_df2['cell_num']>= 50) & (adata_summary_df2['cell_num']<adata_summary_df2['qc_cn']) ]
dataset_index2 = list(adata_summary_df2.index)
adata_list_merge_new2 = [adata_list_merge2[i] for i in dataset_index2]
del adata_list_merge2
adata_list_merge_new = adata_list_merge_new1 + adata_list_merge_new2
del adata_list_merge_new1
del adata_list_merge_new2

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
    resolution = 2
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
# sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")   ### PCA
sc.external.pp.bbknn(adata, batch_key="GSM")    ###  GSM 是batch
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs) #### Find neighborhoods, 用BBKNN做integration，不需要Find Neighbors
sc.tl.umap(adata)
sc.tl.leiden( adata, resolution=resolution,random_state=0,  # flavor="igraph",  # 加上这个参数会报错,
    n_iterations=2, directed=False)
adata.write('./data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge_GSMbatch.h5ad')


# adata.write('./small_file/Fibro_human_10x_TissueNoSplit_scatlasvae.h5ad')


########################## 换一种去 batch 的方法,以 tissue 为 batch   ###########
#################### 换一种batch removal的方法，用文文那个试一下
import bbknn
# adata = adata.raw.to_adata()
adata.var_names_make_unique()   #### gene name unique
adata.obs_names_make_unique()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)   # 不直接挑选有多少个基因
# sc.pp.highly_variable_genes(adata, n_top_genes=8000)   # 细胞数比较多，多选一些 top_genes
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
# 需要重新矫正批次
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
print("Start BBKNN")
bbknn.bbknn(adata, batch_key='tissue', n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8,key_added='leiden0.8')
sc.tl.leiden(adata, resolution=0.9,key_added='leiden0.9')
sc.tl.leiden(adata, resolution=1 ,key_added='leiden1')
sc.tl.leiden(adata, resolution=1.1 ,key_added='leiden1.1')
sc.tl.leiden(adata, resolution=1.2 ,key_added='leiden1.2')
sc.tl.leiden(adata, resolution=1.3 ,key_added='leiden1.3')
sc.tl.leiden(adata, resolution=1.4 ,key_added='leiden1.4')
sc.tl.leiden(adata, resolution=1.5 ,key_added='leiden1.5')
sc.tl.leiden(adata, resolution=1.6 ,key_added='leiden1.6')
sc.tl.leiden(adata, resolution=1.7 ,key_added='leiden1.7')
sc.tl.leiden(adata, resolution=1.8 ,key_added='leiden1.8')
sc.tl.leiden(adata, resolution=1.9 ,key_added='leiden1.9')
sc.tl.leiden(adata, resolution=2 ,key_added='leiden2')





############### 之前那种去 batch 的方法有好些 tissue没 mix到一起
############### Solution: merge the small cell number tissue into One tissue
def quick_merge_small_tissues(adata, threshold=5000, mixed_label='mix'):
    if 'tissue_original' not in adata.obs.columns:
        adata.obs['tissue_original'] = adata.obs['tissue'].copy()
    tissue_counts = adata.obs['tissue'].value_counts()
    small_tissues = tissue_counts[tissue_counts < threshold].index
    mapping_dict = {tissue: (mixed_label if tissue in small_tissues else tissue)
                   for tissue in adata.obs['tissue'].unique()}
    adata.obs['tissue'] = adata.obs['tissue'].map(mapping_dict)
    print(f"合并了 {len(small_tissues)} 个细胞数少于 {threshold} 的tissue")
    print(f"新的tissue分布:")
    print(adata.obs['tissue'].value_counts())
    return adata

adata = quick_merge_small_tissues(adata, threshold=5000, mixed_label='mix')  ## 然后再跑上面的 bbknn 过程
import bbknn
# adata = adata.raw.to_adata()
adata.var_names_make_unique()   #### gene name unique
adata.obs_names_make_unique()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)   # 不直接挑选有多少个基因
# sc.pp.highly_variable_genes(adata, n_top_genes=8000)   # 细胞数比较多，多选一些 top_genes
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
# 需要重新矫正批次
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
print("Start BBKNN")
bbknn.bbknn(adata, batch_key='tissue', n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=resolution,random_state=0, n_iterations=2, directed=False)
