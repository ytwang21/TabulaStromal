######################################################### 下面开始用stream来做 tragectory分析         #########
import h5py
import numpy as np
import scanpy as sc
import scipy.sparse as sp_sparse
import os
import pandas as pd
import gzip
import time
import re
#from pyensembl import EnsemblRelease
import pwd
import itertools
#### scanpy parameters
import omicverse as ov
import sys
import pandas as pd
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
import scvelo as scv
import cellrank as cr
from anndata import AnnData
import matplotlib.pyplot as plt
import dynamo as dyn
from scipy.stats import entropy
from sklearn.neighbors import NearestNeighbors
import stream as st
############################################################# 一、Determine the threshhold
import h5py
import numpy as np
import scanpy as sc
import scipy.sparse as sp_sparse
import os
import pandas as pd
import gzip
import time
import re
from pyensembl import EnsemblRelease
import pwd
import itertools
import glob
import json
import matplotlib.pyplot as plt
#### scanpy parameters
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

####### list the TISCH2 dataset with the expression.h5 & meta file
os.chdir('/mnt/Mercury/home/shixiaoying/Project/HUSCH_table/HUSCH_data')
husch_dataset = os.listdir('/mnt/Mercury/home/shixiaoying/Project/HUSCH_table/HUSCH_data')
husch_dataset = sorted([item for item in husch_dataset if item.startswith('HU_')])
directory ='/mnt/Mercury/home/shixiaoying/Project/HUSCH_table/HUSCH_data/'
Fib_marker_general = sorted(list(set(["COL1A1","COL1A2","COL3A1","COL6A1","IL6","MCAM","MMP2","MYL9", "MYLK","PDGFA"])))
Fib_marker_spec =sorted(list(set(['COL1A1','CFD','CTHRC1','POSTN','FAP','LRRC15','HAND2','SFRP1','IL6','CXCL3','APOD','PI16',
  'GPX1','IGF2','CCL19','SAA1','VWA1','CCL5','CCL4','CILP','THBS4','MYH11','MYLK','GJA4','ACTA2','RGS5'])))  # 'PTPRC'
Fib_marker_yilv = ['IGF2','ENG','FAP','PDGFRA','PECAM1','ITGB1','SPON1','MMP2','MYO1B','ARHGDIB','NDUFA4L2','HIGD1B','COX4I2','RGS5',
    'CNN1','DES','FGFR4','VEGFD','MME','LEFR','PI16','C3','CFD','C7','ADH1B','GREM1','CTHRC1','LRRC15','POSTN','MYH11','ACTA2','PDPN',
    'CXCL14','CXCL12','CCL2','COL12A1','COL11A1','COL10A1','COL8A1','COL5A1','COL4A2','COL4A1','COL3A1','COL1A2','COL1A1'] # ,'CD34'
Fib_marker_union = list(set(Fib_marker_general) | set(Fib_marker_spec) | set(Fib_marker_yilv))
Endo_marker_general =["PECAM1","VWF" ,"CHGB",'CD144','CDH5','KDR','FLK1','CD54','CD105','CD62E','SELE']      # ["PECAM1","VWF" ,"CHGB"]
Endo_marker_spec = ['NR2F2','TCIM','VCAM1','RGCC','APLN','ESM1','BMX','CXCR4','IGKC','SFN','SPARC','PROX1','PDPN','CHGB']  # 'PTPRC','ACTA2','CCL4','GJA4'
Endo_marker_union = list(set(Endo_marker_general).union(set(Endo_marker_spec)))
MoMac_marker_general = ["CD14","CD68" ,"CSF1R","CST3","FABP4","MMP9","S100A9"] #
MoMac_marker_spec = ["FCN1","VCAN","S100A8","S100A9","CD14","FCGR3A","CX3CR1","CD68","MERTK","PPARG","MME","FABP4","MCEMP1","MARCO",
            "MRC1","CD163","CSF1R","GPNMB","APOE","SPP1","SELENOP","TREM2","MSR1","VEGFA","HIF1A","SERPINB9","TNF","CD86","CD40","HLA-DRA","HLA-DQA1","HLA-DPA1",
             "IL1B","CXCL3","CXCL8","CXCL9","CXCL10","CXCL11","CCL18","IDO1","CD274","SIRPA","ENTPD1","VSIR",
            "C1QC","C1QA","C3","CD55","CFP","SLC40A1","IL2RA","CCR2","FOLR2",'CSF3']
MoMac_marker_union = list(set(MoMac_marker_general).union(set(MoMac_marker_spec)))


def find_indices_of_substring(lst, substring):
    indices = []
    for index, element in enumerate(lst):
        if substring in element:
            indices.append(index)
    return indices

mouse_index = find_indices_of_substring(husch_dataset, 'mouse')
def list_files(directory):
    return [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(directory)) for f in fn]

# 筛选出所有.h5文件
def list_h5_files(directory):
    return [f for f in list_files(directory) if f.endswith('_gene_count.h5')]

husch_h5_files = list_h5_files(directory)  ## 只保留 HU_开头的
husch_h5_files = sorted([item for item in husch_h5_files if 'HU_' in item])

################ get rank genes ###############
adata = sc.read_h5ad('./data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge_GSMbatch.h5ad')
############## 查看每个cluster的 top5 DEGs
sc.tl.rank_genes_groups(adata, groupby="leiden",  method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden", standard_scale="var", n_genes=5
)
plt.savefig('./data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge_GSMbatch_degs.png', dpi=300, bbox_inches='tight')
plt.close()

groups = adata.obs['leiden'].cat.categories
df = sc.get.rank_genes_groups_df(adata, group=groups)
df.to_csv('./data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge_GSMbatch_degs.tsv',sep='\t',index=False)
