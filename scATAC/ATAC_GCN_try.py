#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: ATAC_GCN.py
# @time: 2021/8/24 14:39

from time import time
import os
import scanpy as sc
import episcanpy.api as epi
import pandas as pd
import numpy as np
import anndata as ad
import torch
from torch_geometric.datasets import TUDataset
from torch_geometric.data import DataLoader
import torch.nn as nn
import torch.nn.functional as nn_func
import torch_geometric.nn as geo_nn


path_human_brain = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_brain'
path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
file_csv = os.path.join(path_human_brain, 'GSM5289636_s3atac.hg38.counts.csv')
df_csv = pd.read_csv(file_csv, index_col=0)
df_csv = df_csv.T
colnames = [x.replace('-', ':', 1) for x in df_csv.columns]
df_csv.columns = colnames
file_count_tsv = os.path.join(path_data_root, 'counts.tsv')
df_csv.to_csv(file_count_tsv, sep='\t')

file_meta_csv = os.path.join(path_human_brain, 'GSM5289636_s3atac.hg38.metadata.csv')
df_meta_csv = pd.read_csv(file_meta_csv, index_col=0)
df_meta_csv = df_meta_csv.loc[:, ['cellID', 'celltype']]
file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
df_meta_csv.to_csv(file_meta_tsv, sep='\t')

adata = ad.read_text(file_count_tsv, delimiter='\t', first_column_names=True, dtype='int')
df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
adata.obs['celltype'] = df_meta.loc[adata.obs.index, 'celltype']

print(np.max(adata.X))
if np.max(adata.X) > 1:
    epi.pp.binarize(adata)
    print(np.max(adata.X))
epi.pp.filter_cells(adata, min_features=1)
epi.pp.filter_features(adata, min_cells=1)
# QC
adata.obs['log_nb_features'] = [np.log10(x) for x in adata.obs['nb_features']]
epi.pl.violin(adata, ['nb_features'])
epi.pl.violin(adata, ['log_nb_features'])
epi.pp.coverage_cells(adata, binary=False, log=False, bins=50, threshold=50000, save=None)
epi.pp.coverage_cells(adata, binary=False, log=10, bins=50, threshold=1000, save=None)

epi.pp.coverage_features(adata, binary=False, log=False, threshold=5)
epi.pp.coverage_features(adata, binary=False, log=True, threshold=5)

min_features = 1000
max_features = 50000
epi.pp.filter_cells(adata, min_features=min_features)
epi.pp.filter_cells(adata, max_features=max_features)
min_cells = 5
epi.pp.filter_features(adata, min_cells=min_cells)

# save the raw matrix
adata_raw = adata.copy()

epi.pp.cal_var(adata)
min_score_value = 0.515
nb_feature_selected = 100000
epi.pl.variability_features(adata,log=None,
                            min_score=min_score_value, nb_features=nb_feature_selected,
                            save=None)
# create a new AnnData containing only the most variable features
adata = epi.pp.select_var_feature(adata,
                                  nb_features=nb_feature_selected,
                                  show=False,
                                  copy=True)

# save the current version of the matrix (binary, not normalised) in a layer of the Anndata.
adata.layers['binary'] = adata.X.copy()
# normalization
sc.pp.normalize_total(adata)
# save the current version of the matrix (normalised) in a layer of the Anndata.
adata.layers['normalised'] = adata.X.copy()

epi.pp.lazy(adata)
epi.pl.pca_overview(adata, color=['nb_features', 'celltype'])
epi.pl.umap(adata, color=['nb_features', 'celltype'], wspace=0.3)


adata = adata_raw.copy()
adata.layers['raw'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=100000, flavor='seurat')
adata = adata[:, adata.var.highly_variable]
# sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['nb_features', 'celltype'])

# bed
file_peaks = os.path.join(path_data_root, 'peaks.bed')
fmt_peak = "{chrom}\t{start}\t{end}\t{peak_id}\n"
with open(file_peaks, 'w') as w_peak:
    for one_peak in adata_raw.var.index:
        chrom = one_peak.strip().split(':')[0]
        locs = one_peak.strip().split(':')[1]
        start = locs.strip().split('-')[0]
        end = locs.strip().split('-')[1]
        peak_id = one_peak
        w_peak.write(fmt_peak.format(**locals()))

file_peaks_sort = os.path.join(path_data_root, 'peaks.sort.bed')
os.system(f"bedtools sort -i {file_peaks} > {file_peaks_sort}")

path_process = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain/processed_data'
file_gene_hg38 = '/root/scATAC/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
file_peaks_promoter = os.path.join(path_process, 'peaks_promoter.txt')
os.system(f"bedtools intersect -a {file_peaks_sort} -b {file_gene_hg38} -wao "
          f"> {file_peaks_promoter}")


