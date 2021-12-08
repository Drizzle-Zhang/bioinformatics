#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: ATAC_GCN.py
# @time: 2021/8/24 14:39

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
import torch
from torch_geometric.datasets import TUDataset
from torch_geometric.data import InMemoryDataset, Data, DataLoader
import torch.nn as nn
import torch.nn.functional as nn_func
import torch_geometric.nn as geo_nn
torch.multiprocessing.set_sharing_strategy('file_system')


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

    def sum_peaks(self, df_gene_peak, dict_promoter, gene):
        gene_peaks = dict_promoter[gene]
        if len(gene_peaks) == 1:
            df_sum = df_gene_peak.loc[:, gene_peaks].squeeze('columns')
        elif len(gene_peaks) > 1:
            df_sum = np.sum(df_gene_peak.loc[:, gene_peaks], axis=1)
        else:
            return
        df_sum.name = gene
        return df_sum

    def add_promoter(self, file_promoter, num_threads=20):
        if not os.path.exists(self.file_peaks_sort):
            self.generate_peaks_file()
        file_peaks_promoter = os.path.join(self.path_process, 'peaks_promoter.txt')
        os.system(f"bedtools intersect -a {self.file_peaks_sort} -b {file_promoter} -wao "
                  f"> {file_peaks_promoter}")
        dict_promoter = defaultdict(list)
        all_peaks = set()
        with open(file_peaks_promoter, 'r') as w_pro:
            for line in w_pro:
                list_line = line.strip().split('\t')
                if list_line[4] == '.':
                    continue
                gene = list_line[7].strip().split('<-')[0]
                peak = list_line[3]
                dict_promoter[gene].append(peak)
                all_peaks.add(peak)

        all_genes = dict_promoter.keys()
        adata_gene = self.adata[:, [one_peak for one_peak in self.adata.var.index
                                    if one_peak in all_peaks]]
        df_gene_peak = pd.DataFrame(adata_gene.X, index=adata_gene.obs.index,
                                    columns=adata_gene.var.index)
        all_cols = df_gene_peak.columns
        other_cols = set(self.adata.var.index).difference(all_cols)
        self.other_peaks = other_cols
        adata_other = self.adata[:, [one_peak for one_peak in self.adata.var.index
                                     if one_peak in other_cols]]
        adata_other.var['cRE_type'] = np.full(adata_other.n_vars, 'Other')

        pool = Pool(num_threads)
        func_sum = partial(self.sum_peaks, df_gene_peak, dict_promoter)
        result = pool.map(func_sum, all_genes)
        pool.close()
        # result = [one_df for one_df in result if one_df is not None]
        df_gene = pd.concat(result, axis=1)
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
        return df_interaction


def one_test(all_genes, adata_ATAC, adata_merge_peak, array_peak, array_celltype, i):
    # i = 0
    cell = adata_ATAC.obs.index[i]
    label = adata_ATAC.obs.loc[cell, 'celltype']
    label_idx = torch.tensor(np.argwhere(array_celltype == label)[0], dtype=torch.int64)
    cell_cluster = np.nonzero(adata_ATAC.obsp['connectivities'].getcol(i).toarray())[0]
    cell_cluster = np.append(cell_cluster, i)
    df_cell = adata_merge_peak[cell_cluster, :].to_df()
    list_df = []
    for gene in all_genes:
        sub_distal = df_graph.loc[df_graph['region1'] == gene, 'region2']
        df_distal = df_cell.loc[:, sub_distal]
        if type(df_distal) == 'pandas.core.series.Series':
            df_distal = df_cell.loc[:, sub_distal].to_frame()
        vec_corr = df_distal.corrwith(df_cell.loc[:, gene], method='spearman')
        df_corr = pd.DataFrame({'gene': np.full(len(sub_distal), gene), 'distal': vec_corr.index,
                                'corr': vec_corr.tolist()})
        list_df.append(df_corr)

    df_corr_cell = pd.concat(list_df).fillna(0)
    array_region1 = df_corr_cell['gene'].apply(lambda x: np.argwhere(array_peak == x)[0, 0])
    array_region2 = df_corr_cell['distal'].apply(lambda x: np.argwhere(array_peak == x)[0, 0])
    df_graph_index = torch.tensor([np.array(array_region1), np.array(array_region2)],
                                  dtype=torch.int64)
    cell_data = Data(x=torch.reshape(torch.Tensor(df_cell.loc[cell, :].T), (df_cell.shape[1], 1)),
                     edge_index=df_graph_index,
                     edge_attr=torch.reshape(torch.Tensor(df_corr_cell['corr'].tolist()),
                                             (df_corr_cell.shape[0], 1)),
                     y=label_idx)
    print(cell)

    return cell_data


def generate_data_list(graph_data, dataset_atac, num_threads):
    all_genes = set(graph_data['region1'])
    all_peaks = set(graph_data['region1']).union(set(graph_data['region2']))

    adata_ATAC = dataset_atac.adata
    adata_merge = dataset_atac.adata_merge
    adata_merge_peak = adata_merge[:, [one_peak for one_peak in adata_merge.var.index
                                       if one_peak in all_peaks]]
    array_peak = np.array(adata_merge_peak.var.index)
    array_celltype = np.unique(np.array(adata_ATAC.obs['celltype']))

    list_idx = list(range(0, adata_ATAC.n_obs))
    pool = Pool(num_threads)
    func_test = partial(one_test, all_genes, adata_ATAC, adata_merge_peak,
                        array_peak, array_celltype)
    result = pool.map(func_test, list_idx)
    pool.close()

    return result, array_peak, array_celltype


class ATACGraphDataset(InMemoryDataset):
    def __init__(self, root, data_list):
        self.data_list = data_list
        super(ATACGraphDataset, self).__init__(root)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ['some_file_1']

    @property
    def processed_file_names(self):
        return ['data.pt']

    def download(self):
        pass

    def process(self):
        # Read data into huge `Data` list.
        data_list = self.data_list
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])


class GCN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super(GCN, self).__init__()
        torch.manual_seed(12345)
        self.conv1 = geo_nn.GCNConv(dataset.num_node_features, hidden_channels)
        self.conv2 = geo_nn.GCNConv(hidden_channels, hidden_channels)
        self.conv3 = geo_nn.GCNConv(hidden_channels, hidden_channels)
        self.lin = nn.Linear(hidden_channels, dataset.num_classes)

    def forward(self, x, edge_index, batch):
        # 1. Obtain node embeddings
        x = self.conv1(x, edge_index)
        x = x.relu()
        x = self.conv2(x, edge_index)
        x = x.relu()
        x = self.conv3(x, edge_index)

        # 2. Readout layer
        x = geo_nn.global_mean_pool(x, batch)  # [batch_size, hidden_channels]

        # 3. Apply a final classifier
        x = nn_func.dropout(x, p=0.5, training=self.training)
        x = self.lin(x)

        return x


def train(loader):
    model.train()
    for data in loader:  # Iterate in batches over the training dataset.
        data = data.to(device)
        out = model(data.x, data.edge_index, data.batch)  # Perform a single forward pass.
        loss = criterion(out, data.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.


def test(loader):
    model.eval()
    correct = 0
    for data in loader:  # Iterate in batches over the training/test dataset.c
        data = data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        correct += int((pred == data.y).sum())  # Check against ground-truth labels.
    return correct / len(loader.dataset)  # Derive ratio of correct predictions.


if __name__ == '__main__':
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    dataset_ATAC = ATACDataset(data_root=path_data_root, raw_filename='counts.tsv')
    file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
    dataset_ATAC.adata.obs['celltype'] = df_meta.loc[dataset_ATAC.adata.obs.index, 'celltype']

    dataset_ATAC.quality_control()
    dataset_ATAC.find_neighbors()
    # adata_ATAC.plot_umap()
    file_gene_hg38 = '/root/scATAC/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_ATAC.add_promoter(file_gene_hg38)

    path_hic = '/root/scATAC/pcHi-C/Interactions_by_tissue/Dorsolateral_Prefrontal_Cortex'
    df_graph = dataset_ATAC.build_graph(path_hic)
    time_end = time()
    print(time_end - time_start)

    time_start = time()
    list_graph_data, peaks, celltypes = generate_data_list(df_graph, dataset_ATAC, 30)
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input, list_graph_data)
    time_end = time()
    print(time_end - time_start)

    # train model
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:800]
    test_dataset = dataset[800:]

    model = GCN(hidden_channels=64).to(device)
    print(model)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=800, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=200, shuffle=False)

    for epoch in range(1, 201):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')



