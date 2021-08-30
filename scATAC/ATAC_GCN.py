#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: ATAC_GCN.py
# @time: 2021/8/24 14:39

from time import time
import os
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


class ATACDataset(object):
    def __init__(self, data_root, raw_filename, meta_filename, n_pc=50,
                 n_neigh=30):
        self.data_root = data_root
        self.raw_filename = raw_filename
        self.meta_filename = meta_filename

        self.adata = self.load_matrix()
        self.n_pc = n_pc
        self.n_neigh = n_neigh

        self.adata = self.preprocess()

    def load_matrix(self):
        adata = ad.read_text(
            os.path.join(self.data_root, self.raw_filename),
            delimiter='\t', first_column_names=True, dtype='int')
        adata = epi.pp.load_metadata(adata, os.path.join(self.data_root, self.meta_filename))
        return adata

    def preprocess(self):
        adata = epi.pp.normalize_total(self.adata)
        adata = epi.pp.pca(adata, self.n_pc)
        adata = epi.pp.neighbors(
            adata, n_neighbors=self.n_neigh, n_pcs=self.n_pc, metric='cosine')
        return adata


if __name__ == '__main__':
    time_start = time()
    path_human_brain = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_brain'
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    file_csv = os.path.join(path_human_brain,
                            'GSM5289636_s3atac.hg38.counts.csv')
    df_csv = pd.read_csv(file_csv, index_col=0)
    df_csv = df_csv.T
    file_count_tsv = os.path.join(path_data_root, 'counts.tsv')
    df_csv.to_csv(file_count_tsv, sep='\t')

    file_meta_csv = os.path.join(path_human_brain,
                                 'GSM5289636_s3atac.hg38.metadata.csv')
    df_meta_csv = pd.read_csv(file_meta_csv, index_col=0)
    df_meta_csv = df_meta_csv.loc[:, ['cellID', 'celltype']]
    file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
    df_meta_csv.to_csv(file_meta_tsv, sep='\t')

    adata_ATAC = ATACDataset(
        data_root=path_data_root, raw_filename='counts.tsv',
        meta_filename='metadata.tsv')

    adata = ad.read_text(file_count_tsv,
                         delimiter='\t', first_column_names=True, dtype='int')
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
    epi.pp.coverage_cells(adata, binary=False, log=False, bins=50,
                          threshold=50000, save=None)
    epi.pp.coverage_cells(adata, binary=False, log=10, bins=50,
                          threshold=1000, save=None)

    epi.pp.coverage_features(adata, binary=False, log=False, threshold=5)
    epi.pp.coverage_features(adata, binary=False, log=True, threshold=5)

    min_features = 1000
    max_features = 50000
    epi.pp.filter_cells(adata, min_features=min_features)
    epi.pp.filter_cells(adata, max_features=max_features)
    min_cells = 5
    epi.pp.filter_features(adata, min_cells=min_cells)

    epi.pp.cal_var()

    time_end = time()
    print(time_end - time_start)
