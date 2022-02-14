#!/usr/bin/env python
# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: ATACProcess.py
# @time: 2022/1/6 14:54

from time import time
import os
import episcanpy.api as epi
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
import pickle
import torch
from torch_geometric.datasets import TUDataset
from torch_geometric.data import InMemoryDataset, Data, DataLoader
import torch.nn as nn
import torch.nn.functional as nn_func
import torch_geometric.nn as geo_nn
from captum.attr import Saliency, IntegratedGradients
import random


class ATACDataset(object):
    def __init__(self, data_root, raw_filename):
        self.data_root = data_root
        self.raw_filename = raw_filename
        self.adata = self.load_matrix()
        self.path_process = os.path.join(data_root, 'processed_files')
        if not os.path.exists(self.path_process):
            os.mkdir(self.path_process)
        self.file_peaks_sort = os.path.join(self.path_process, 'peaks.sort.bed')
        self.all_genes = None
        self.adata_merge = None
        self.other_peaks = None
        self.df_graph = None
        self.list_graph = None
        self.array_peak = None
        self.array_celltype = None

    def load_matrix(self):
        adata = ad.read_text(
            os.path.join(self.data_root, self.raw_filename),
            delimiter='\t', first_column_names=True, dtype='int')
        return adata

    def generate_peaks_file(self):
        file_peaks = os.path.join(self.path_process, 'peaks.bed')
        fmt_peak = "{chrom}\t{start}\t{end}\t{peak_id}\n"
        with open(file_peaks, 'w') as w_peak:
            for one_peak in self.adata.var.index:
                chrom = one_peak.strip().split(':')[0]
                locs = one_peak.strip().split(':')[1]
                start = locs.strip().split('-')[0]
                end = locs.strip().split('-')[1]
                peak_id = one_peak
                w_peak.write(fmt_peak.format(**locals()))

        os.system(f"bedtools sort -i {file_peaks} > {self.file_peaks_sort}")

    def quality_control(self, min_features=1000, max_features=50000, min_cells=5):
        epi.pp.filter_cells(self.adata, min_features=min_features)
        epi.pp.filter_cells(self.adata, max_features=max_features)
        epi.pp.filter_features(self.adata, min_cells=min_cells)

    def find_neighbors(self, num_peak=120000, num_pc=50, num_neighbor=30):
        adata = self.adata
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=num_peak, flavor='seurat')
        adata = adata[:, adata.var.highly_variable]
        # sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack', n_comps=num_pc)
        sc.pp.neighbors(adata, n_neighbors=num_neighbor, n_pcs=num_pc, metric='cosine', knn=True)
        self.adata = adata

    def plot_umap(self):
        adata = self.adata
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
        sc.tl.umap(adata)
        out_plot = sc.pl.umap(self.adata, color=['nb_features', 'celltype'])
        return out_plot

    def add_promoter(self, file_promoter):
        if not os.path.exists(self.file_peaks_sort):
            self.generate_peaks_file()
        file_peaks_promoter = os.path.join(self.path_process, 'peaks_promoter.txt')
        os.system(f"bedtools intersect -a {self.file_peaks_sort} -b {file_promoter} -wao "
                  f"> {file_peaks_promoter}")
        dict_promoter = defaultdict(list)
        with open(file_peaks_promoter, 'r') as w_pro:
            for line in w_pro:
                list_line = line.strip().split('\t')
                if list_line[4] == '.':
                    continue
                gene = list_line[7].strip().split('<-')[0]
                peak = list_line[3]
                dict_promoter[gene].append(peak)

        all_genes = dict_promoter.keys()
        list_peaks_1 = []
        list_genes_1 = []
        list_peaks_2 = []
        list_genes_2 = []
        for gene in all_genes:
            sub_peaks = dict_promoter[gene]
            if len(sub_peaks) == 1:
                list_peaks_1.extend(sub_peaks)
                list_genes_1.append(gene)
            else:
                list_genes_2.extend([gene for _ in range(len(sub_peaks))])
                list_peaks_2.extend(sub_peaks)
        adata_gene_1 = self.adata[:, list_peaks_1]
        df_gene_peak_1 = pd.DataFrame(adata_gene_1.X, index=adata_gene_1.obs.index,
                                      columns=list_genes_1)
        adata_gene_2 = self.adata[:, list_peaks_2]
        df_gene_peak_2 = pd.DataFrame(
            adata_gene_2.X,
            index=adata_gene_2.obs.index,
            columns=pd.MultiIndex.from_arrays([list_genes_2, list_peaks_2], names=['gene', 'peak']))
        df_gene_peak_2_t = df_gene_peak_2.T
        df_gene_peak_2_t_gene = df_gene_peak_2_t.groupby('gene').apply(lambda x: x.sum())
        df_gene_peak_2 = df_gene_peak_2_t_gene.T
        all_cols = set(list_peaks_1 + list_peaks_2)
        other_cols = set(self.adata.var.index).difference(all_cols)
        self.other_peaks = other_cols
        adata_other = self.adata[:, [one_peak for one_peak in self.adata.var.index
                                     if one_peak in other_cols]]
        adata_other.var['cRE_type'] = np.full(adata_other.n_vars, 'Other')

        df_gene = pd.concat([df_gene_peak_1, df_gene_peak_2], axis=1)
        adata_promoter = \
            ad.AnnData(X=df_gene,
                       var=pd.DataFrame(data={'cRE_type': np.full(df_gene.shape[1], 'Promoter')},
                                        index=df_gene.columns),
                       obs=pd.DataFrame(index=df_gene.index))
        self.all_genes = set(df_gene.columns)
        adata_merge = ad.concat([adata_promoter, adata_other], axis=1)
        self.adata_merge = adata_merge

        return

    def build_graph(self, path_interaction):
        file_pp = os.path.join(path_interaction, 'PP.txt')
        file_po = os.path.join(path_interaction, 'PO.txt')
        df_pp = pd.read_csv(file_pp, sep='\t', header=None)
        df_pp = df_pp.loc[
                df_pp.apply(lambda x:
                            x.iloc[0] in self.all_genes and x.iloc[1] in self.all_genes, axis=1), :]
        df_pp.columns = ['region1', 'region2']
        file_po_peaks = os.path.join(self.path_process, 'peaks_PO.bed')
        os.system(f"bedtools intersect -a {self.file_peaks_sort} -b {file_po} -wao "
                  f"> {file_po_peaks}")
        list_dict = []
        with open(file_po_peaks, 'r') as r_po:
            for line in r_po:
                list_line = line.strip().split('\t')
                peak = list_line[3]
                gene = list_line[8]
                if peak in self.other_peaks and gene in self.all_genes:
                    list_dict.append({"region1": gene, "region2": peak})
        df_po = pd.DataFrame(list_dict)
        df_interaction = pd.concat([df_pp, df_po])
        self.df_graph = df_interaction.drop_duplicates()

        return

    def generate_data_list(self):
        graph_data = self.df_graph
        adata_atac = self.adata
        adata_merge = self.adata_merge
        all_peaks = set(graph_data['region1']).union(set(graph_data['region2']))
        adata_merge_peak = adata_merge[:, [one_peak for one_peak in adata_merge.var.index
                                           if one_peak in all_peaks]]
        array_peak = np.array(adata_merge_peak.var.index)
        array_celltype = np.unique(np.array(adata_atac.obs['celltype']))
        array_region1 = graph_data['region1'].apply(lambda x: np.argwhere(array_peak == x)[0, 0])
        array_region2 = graph_data['region2'].apply(lambda x: np.argwhere(array_peak == x)[0, 0])
        df_graph_index = torch.tensor([np.array(array_region1), np.array(array_region2)],
                                      dtype=torch.int64)
        df_merge_peak = adata_merge_peak.to_df()
        list_graph = []
        for i in range(0, adata_atac.n_obs):
            cell = adata_atac.obs.index[i]
            label = adata_atac.obs.loc[cell, 'celltype']
            label_idx = torch.tensor(np.argwhere(array_celltype == label)[0], dtype=torch.int64)
            cell_data = Data(x=torch.reshape(torch.Tensor(df_merge_peak.loc[cell, :]),
                                             (adata_merge_peak.shape[1], 1)),
                             edge_index=df_graph_index, y=label_idx)
            list_graph.append(cell_data)

        self.list_graph = list_graph
        self.array_peak = array_peak
        self.array_celltype = array_celltype

        return


if __name__ == '__main__':
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC'
    dataset_ATAC = ATACDataset(data_root=path_data_root, raw_filename='counts.tsv')
    file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
    dataset_ATAC.adata.obs['celltype'] = df_meta.loc[dataset_ATAC.adata.obs.index, 'celltype']

    dataset_ATAC.quality_control(min_features=3000)
    file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_ATAC.add_promoter(file_gene_hg38)
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    dataset_ATAC.build_graph(path_hic)
    dataset_ATAC.generate_data_list()

    # dataset_ATAC.find_neighbors()
    # dataset_ATAC.plot_umap()
    time_end = time()
    print(time_end - time_start)
