#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: GraphConv.py
# @time: 2021/9/9 11:58

from time import time
import os
from torch_geometric.data import InMemoryDataset, DataLoader
import networkx as nx
import numpy as np
from torch_geometric.utils import to_networkx
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, global_add_pool, global_mean_pool, TopKPooling
from captum.attr import Saliency, IntegratedGradients
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from scipy.stats import kstest
import episcanpy.api as epi
import scanpy as sc
import pandas as pd


class ATACGraphDataset(InMemoryDataset):
    def __init__(self, root, data_list=None):
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
    def __init__(self, input_channels, output_channels, hidden_channels, num_nodes):
        super(GCN, self).__init__()
        torch.manual_seed(12345)

        self.conv1 = GraphConv(input_channels, hidden_channels//4)
        # self.conv2 = GraphConv(hidden_channels//4, hidden_channels)
        # self.conv3 = GraphConv(hidden_channels, hidden_channels//4)
        # self.conv4 = GraphConv(hidden_channels, hidden_channels)
        self.conv5 = GraphConv(hidden_channels//4, input_channels)

        self.lin1 = nn.Linear(num_nodes, num_nodes//5)
        self.lin2 = nn.Linear(num_nodes//5, num_nodes//25)
        self.lin3 = nn.Linear(num_nodes//25, num_nodes//125)
        self.lin4 = nn.Linear(num_nodes//125, output_channels)

    def forward(self, x, edge_index, batch, edge_weight=None):
        x = self.conv1(x, edge_index, edge_weight).relu()
        # x = self.conv2(x, edge_index, edge_weight).relu()
        # x = self.conv3(x, edge_index, edge_weight).relu()
        # x = self.conv4(x, edge_index, edge_weight).relu()
        x = self.conv5(x, edge_index, edge_weight)
        batch_size = len(torch.unique(batch))
        x = x.view(batch_size, x.shape[0]//batch_size)
        x = self.lin1(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin2(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin3(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin4(x)
        return F.log_softmax(x, dim=1)


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


def model_forward(edge_mask, data):
    batch = torch.zeros(data.x.shape[0], dtype=int).to(device)
    out = model(data.x, data.edge_index, batch, edge_mask)
    return out


def explain(method, data, target=0):
    input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
    if method == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=target,
                            additional_forward_args=(data,),
                            internal_batch_size=data.edge_index.shape[1])
    elif method == 'saliency':
        saliency = Saliency(model_forward)
        mask = saliency.attribute(input_mask, target=target,
                                  additional_forward_args=(data,))
    else:
        raise Exception('Unknown explanation method')

    edge_mask = np.abs(mask.cpu().detach().numpy())
    if edge_mask.max() > 0:  # avoid division by zero
        edge_mask = edge_mask / edge_mask.max()
    return edge_mask


def aggregate_edge_directions(edge_mask, data):
    edge_mask_dict = defaultdict(float)
    for val, u, v in list(zip(edge_mask, *data.edge_index)):
        u, v = u.item(), v.item()
        if u > v:
            u, v = v, u
        edge_mask_dict[(u, v)] += val
    return edge_mask_dict


if __name__ == '__main__':
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:1700]
    test_dataset = dataset[1700:]

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = GCN(input_channels=dataset.num_node_features,
                output_channels=dataset.num_classes, hidden_channels=32,
                num_nodes=dataset_atac_graph[0].num_nodes).to(device)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    # train model
    time_start = time()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, '
              f'Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')
        if test_acc > 0.99:
            break
    time_end = time()
    print(time_end - time_start)

    # explain model
    def one_cell(method, device, data):
        data = data.to(device)
        edge_mask = explain(method, data, target=data.y.to(device))
        edge_mask_dict = aggregate_edge_directions(edge_mask, data)
        return edge_mask_dict

    method = 'saliency'
    pool = Pool(20)
    func_explain = partial(one_cell, method, device)
    list_dict = pool.map(func_explain, dataset_atac_graph)
    pool.close()

    adata_edge = ad.AnnData(X=df_weight, obs=dataset_ATAC.adata.obs)
    adata_edge.raw = adata_edge.copy()
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata_edge, n_top_genes=10000, flavor='seurat')
    adata = adata_edge[:, adata_edge.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['celltype'])

    # rank edge
    sc.tl.rank_genes_groups(adata_edge, groupby='celltype', method='wilcoxon', use_raw=True,
                            pts=True, tie_correct=True)
    array_names = adata_edge.uns['rank_genes_groups']['names']
    array_scores = adata_edge.uns['rank_genes_groups']['scores']
    array_pval = adata_edge.uns['rank_genes_groups']['pvals_adj']

    file_peaks_celltypes = os.path.join(path_graph_input, 'peaks_celltypes.npz')
    filenpz = np.load(file_peaks_celltypes, allow_pickle=True)
    peaks = filenpz['peaks']
    celltypes = filenpz['celltypes']
    dict_cell_scores = {}
    for idx_cell in range(len(celltypes)):
        cell_name = celltypes[idx_cell]
        list_edge_names = []
        list_edge_scores = []
        for idx_edge in range(array_names.shape[0]):
            edge_name = array_names[idx_edge][idx_cell]
            list_edge_names.append((peaks[edge_name[0]], peaks[edge_name[1]]))
            list_edge_scores.append(array_scores[idx_edge][idx_cell])
        dict_cell_scores[cell_name] = pd.Series(list_edge_scores, index=list_edge_names)

    # cell type specific interactome
    df_graph_pp = df_graph.loc[df_graph['region2'].apply(lambda x: x[:3] != 'chr'), :]
    list_pp = df_graph_pp.apply(lambda x: (x['region1'], x['region2']), axis=1).tolist()
    df_graph_po = df_graph.loc[df_graph['region2'].apply(lambda x: x[:3] == 'chr'), :]
    list_po = df_graph_po.apply(lambda x: (x['region1'], x['region2']), axis=1).tolist()
    path_cell_interatome = os.path.join(path_data_root, 'cell_interatome')
    file_po = os.path.join(path_cell_interatome, 'PO.bed')
    with open(file_po, 'w') as w_po:
        for pair in list_po:
            pair_1 = pair[0]
            pair_2 = pair[1]
            chrom = pair_2.split(':')[0]
            start = pair_2.split(':')[1].split('-')[0]
            end = pair_2.split(':')[1].split('-')[1]
            w_po.write(f"{chrom}\t{start}\t{end}\t{pair_1}\t{pair_2}\n")

    list_cell = ['Microglia', 'Neuron', 'Oligo']
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    list_df_pp = []
    list_df_po = []
    for cell in list_cell:
        path_cell_hg38 = os.path.join(path_hic, cell)
        file_pp_cell = os.path.join(path_cell_hg38, 'PP.txt')
        file_po_cell = os.path.join(path_cell_hg38, 'PO.txt')
        path_cell = os.path.join(path_cell_interatome, cell)
        if not os.path.exists(path_cell):
            os.mkdir(path_cell)
        file_intersect = os.path.join(path_cell, 'PO_intersect.txt')
        os.system(f"bedtools intersect -a {file_po} -b {file_po_cell} -wao > {file_intersect}")
        df_intersect = pd.read_csv(file_intersect, sep='\t', header=None)
        df_po_cell = df_intersect.loc[df_intersect.iloc[:, 5] != '.', [3, 4]]
        df_po_cell = df_po_cell.drop_duplicates()
        df_po_cell.columns = ['region1', 'region2']
        df_pp_cell = pd.read_csv(file_pp_cell, sep='\t', header=None)
        df_pp_cell.columns = ['region1', 'region2']
        set_pp_cell = set(df_pp_cell.apply(lambda x: (x['region1'], x['region2']), axis=1).tolist())
        list_df_pp.append(pd.Series([pair in set_pp_cell for pair in list_pp], index=list_pp))
        set_po_cell = set(df_po_cell.apply(lambda x: (x['region1'], x['region2']), axis=1).tolist())
        list_df_po.append(pd.Series([pair in set_po_cell for pair in list_po], index=list_po))

    df_pp = pd.concat(list_df_pp, axis=1)
    df_po = pd.concat(list_df_po, axis=1)
    df_cell = pd.concat([df_pp, df_po], axis=0)
    df_cell.columns = list_cell

    # ks test
    df_pval = pd.DataFrame(np.full(shape=(len(list_cell), len(celltypes)), fill_value=1),
                           index=list_cell, columns=celltypes)
    for cell in list_cell:
        cell_pair = df_cell.index[df_cell[cell]]
        for celltype in celltypes:
            all_score = dict_cell_scores[celltype]
            cell_score = all_score[cell_pair]
            df_pval.loc[cell, celltype] = \
                kstest(np.array(cell_score), np.array(all_score), alternative='less')[1]
    file_pvals = os.path.join(path_cell_interatome, 'interactome_pvals.txt')
    df_pval.to_csv(file_pvals, sep='\t')
