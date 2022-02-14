# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: read_data.py
# @time: 2022/2/9 21:34

import os
import scanpy as sc
import pandas as pd

path_save = '/root/scATAC/ATAC_data/Alzheimer_Morabito/raw_data/'
adata = sc.read_h5ad(os.path.join(path_save, 'AD_Control.h5ad'))

df_meta = pd.read_csv(os.path.join(path_save, 'GSE174367_snATAC-seq_cell_meta.csv'))
df_meta.index = df_meta['Barcode']
df_meta['celltype'] = df_meta.apply(lambda x: f"{x['Diagnosis']}_{x['Cell.Type']}", axis=1)
adata = adata[df_meta.index, :]
adata.obs = df_meta.loc[adata.obs.index, :]

adata.write_h5ad(filename=os.path.join(path_save, 'AD_Control_celltype.h5ad'))
