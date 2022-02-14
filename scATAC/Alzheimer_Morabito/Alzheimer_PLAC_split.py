# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: ATAC_GCN.py
# @time: 2021/8/24 14:39


from time import time
import os
from torch_geometric.data import InMemoryDataset, Data, DataLoader
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, TopKPooling
from captum.attr import Saliency, IntegratedGradients
from collections import defaultdict
from scipy.stats import kstest
import episcanpy.api as epi
import scanpy as sc
import pandas as pd
import anndata as ad
import pickle
import captum.attr as attr
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
        if self.raw_filename[-5:] == '.h5ad':
            adata = sc.read_h5ad(os.path.join(self.data_root, self.raw_filename))
        elif self.raw_filename[-4:] == '.tsv':
            adata = ad.read_text(
                os.path.join(self.data_root, self.raw_filename),
                delimiter='\t', first_column_names=True, dtype='int')
        else:
            raise ImportError("Input format error!")
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

    def quality_control(self, min_features=1000, max_features=50000,
                        min_percent=None, min_cells=None):
        self.adata.raw = self.adata.copy()
        adata_atac = self.adata
        epi.pp.filter_cells(adata_atac, min_features=min_features)
        epi.pp.filter_cells(adata_atac, max_features=max_features)
        if min_percent is not None:
            df_count = pd.DataFrame(adata_atac.X, index=adata_atac.obs.index,
                                    columns=adata_atac.var.index)
            array_celltype = np.array(adata_atac.obs['celltype'])
            celltypes = np.unique(array_celltype)
            df_percent = pd.DataFrame(
                np.full(shape=(df_count.shape[1], len(celltypes)), fill_value=0),
                index=adata_atac.var.index, columns=celltypes)
            for cell in celltypes:
                sub_count = df_count.loc[array_celltype == cell, :]
                sub_percent = np.sum(sub_count != 0, axis=0) / sub_count.shape[0]
                df_percent[cell] = sub_percent
            df_percent_max = np.max(df_percent, axis=1)
            sel_peaks = df_percent_max.index[df_percent_max > min_percent]
            self.adata = self.adata[:, sel_peaks]
        elif min_cells is not None:
            epi.pp.filter_features(adata_atac, min_cells=min_cells)

    def select_genes(self, num_peak=120000):
        adata = self.adata
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=num_peak, flavor='seurat')
        self.adata = self.adata[:, adata.var.highly_variable]

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
        if os.path.exists(self.file_peaks_sort):
            os.remove(self.file_peaks_sort)
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
        df_gene_peak_1 = pd.DataFrame(adata_gene_1.X.toarray(), index=adata_gene_1.obs.index,
                                      columns=list_genes_1)
        adata_gene_2 = self.adata[:, list_peaks_2]
        df_gene_peak_2 = pd.DataFrame(
            adata_gene_2.X.toarray(),
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
                       obs=self.adata.obs.loc[df_gene.index, :])
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
                             edge_index=df_graph_index, y=label_idx, cell=cell)
            list_graph.append(cell_data)

        self.list_graph = list_graph
        self.array_peak = array_peak
        self.array_celltype = array_celltype

        return


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

        self.conv1 = GraphConv(input_channels, hidden_channels)
        # self.pool1 = TopKPooling(hidden_channels, ratio=0.5)
        self.conv2 = GraphConv(hidden_channels, hidden_channels)

        # num_nodes = num_nodes//int(1/0.5)
        self.lin1 = nn.Linear(num_nodes, num_nodes//5)
        self.lin2 = nn.Linear(num_nodes//5, num_nodes//25)
        self.lin3 = nn.Linear(num_nodes//25, num_nodes//125)
        self.lin4 = nn.Linear(num_nodes//125, num_nodes//625)
        self.lin5 = nn.Linear(num_nodes//625, output_channels)

    def forward(self, x, edge_index, batch, edge_weight=None):
        x = self.conv1(x, edge_index, edge_weight).relu()
        # x, edge_index, edge_weight, batch, _, _ = self.pool1(x, edge_index, edge_weight, batch)
        x = self.conv2(x, edge_index, edge_weight)
        batch_size = len(torch.unique(batch))
        x = torch.mean(x, dim=1, keepdim=True)
        x = x.view(batch_size, x.shape[0]//batch_size)
        x = self.lin1(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin2(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin3(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin4(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin5(x)
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


def model_forward(edge_mask, data, model):
    out = model(data.x, data.edge_index, data.batch, edge_mask)
    return out


if __name__ == '__main__':
    # AD
    time_start = time()
    path_AD_Control = '/root/scATAC/ATAC_data/Alzheimer_Morabito/AD_Control'
    dataset_AD_Control = ATACDataset(data_root=path_AD_Control, raw_filename='AD_Control.h5ad')
    file_meta_tsv = os.path.join(path_AD_Control, 'metadata.csv')
    df_meta = pd.read_csv(file_meta_tsv)
    df_meta.index = df_meta['Barcode']
    df_meta['celltype'] = df_meta.apply(lambda x: f"{x['Diagnosis']}_{x['Cell.Type']}", axis=1)
    dataset_AD_Control.adata = dataset_AD_Control.adata[df_meta.index, :]
    dataset_AD_Control.adata.obs = df_meta.loc[dataset_AD_Control.adata.obs.index, :]

    dataset_AD_Control.quality_control(min_features=300, min_cells=30)
    dataset_AD_Control.select_genes(num_peak=120000)
    file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_AD_Control.add_promoter(file_gene_hg38)

    # PLAC
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    dataset_AD_Control.build_graph(path_hic)
    df_graph_PLAC = dataset_AD_Control.df_graph

    dataset_AD_Control.generate_data_list()
    list_graph_data = dataset_AD_Control.list_graph
    path_graph_input = os.path.join(path_AD_Control, 'input_graph')
    os.system(f"rm -rf {path_graph_input}")
    os.mkdir(path_graph_input)
    dataset_atac_graph = ATACGraphDataset(path_graph_input, list_graph_data)

    dataset_AD_Control.find_neighbors()
    dataset_AD_Control.plot_umap()
    time_end = time()
    print(time_end - time_start)

    # save data
    file_atac_AD_Control = os.path.join(path_AD_Control, 'dataset_atac.pkl')
    with open(file_atac_AD_Control, 'wb') as w_pkl:
        str_pkl = pickle.dumps(dataset_AD_Control)
        w_pkl.write(str_pkl)

    # read data
    path_AD_Control = '/root/scATAC/ATAC_data/Alzheimer_Morabito/AD_Control'
    file_atac_AD_Control = os.path.join(path_AD_Control, 'dataset_atac.pkl')
    with open(file_atac_AD_Control, 'rb') as r_pkl:
        dataset_AD_Control = pickle.loads(r_pkl.read())

    path_AD_Control = '/root/scATAC/ATAC_data/Alzheimer_Morabito/AD_Control'
    path_graph_input = os.path.join(path_AD_Control, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:100000]
    test_dataset = dataset[100000:]

    device = torch.device("cuda:3" if torch.cuda.is_available() else "cpu")
    model = GCN(input_channels=dataset.num_node_features,
                output_channels=dataset.num_classes, hidden_channels=8,
                num_nodes=dataset_atac_graph[0].num_nodes).to(device)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    # train model
    time_start = time()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
    for epoch in range(1, 50):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, '
              f'Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')
        if test_acc > 0.97:
            break
    time_end = time()
    print(time_end - time_start)

    # explain model
    all_loader = DataLoader(dataset_atac_graph, batch_size=32, shuffle=True)
    list_dict = []
    method = 'ig'
    # i = 0
    for data in all_loader:
        data = data.to(device)
        target = data.y
        input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
        # dl = attr.DeepLift(model_forward)
        # mask = dl.attribute(input_mask, target=target,
        #                     additional_forward_args=(data, model))
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(
            input_mask, target=target, n_steps=50,
            additional_forward_args=(data, model),
            internal_batch_size=data.edge_index.shape[1])
        batch_size = len(torch.unique(data.batch))
        num_col = mask.shape[0]//batch_size
        mask = mask.view(batch_size, num_col)
        edge_mask = np.abs(mask.cpu().detach().numpy())
        edge_mask = edge_mask / np.max(edge_mask, axis=1)[:, np.newaxis]
        sub_edge_index = data.edge_index.cpu().numpy()
        col_edge = [(sub_edge_index[0, i], sub_edge_index[1, i]) for i in range(num_col)]
        list_dict.append(pd.DataFrame(edge_mask, columns=col_edge, index=data.cell))
        # i = i + 1
        # if i >= 10:
        #     break
    df_weight = pd.concat(list_dict)
    # df_weight = df_weight.dropna()
    # df_weight.index = dataset_ATAC.adata.obs.index

    # # save weight
    # file_weight = os.path.join(path_data_root, 'weight_atac.pkl')
    # with open(file_weight, 'wb') as w_pkl:
    #     str_pkl = pickle.dumps(df_weight)
    #     w_pkl.write(str_pkl)
    #
    # # read weight
    # path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_merge_HiC'
    # file_weight = os.path.join(path_data_root, 'weight_atac.pkl')
    # with open(file_weight, 'rb') as r_pkl:
    #     df_weight = pickle.loads(r_pkl.read())

    adata_edge = ad.AnnData(X=df_weight, obs=dataset_ATAC.adata.obs.loc[df_weight.index, :])
    adata_edge.raw = adata_edge.copy()
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata_edge, n_top_genes=10000, flavor='seurat')
    adata = adata_edge[:, adata_edge.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.2)
    sc.pl.umap(adata, color=['nb_features', 'celltype'])

    # sample
    time_start = time()
    path_sample_110 = '/root/scATAC/ATAC_data/Alzheimer_Morabito/AD_Control_110'
    dataset_AD_Control = ATACDataset(data_root=path_sample_110, raw_filename='AD_Control.h5ad')
    file_meta_tsv = os.path.join(path_sample_110, 'metadata.csv')
    df_meta = pd.read_csv(file_meta_tsv)
    df_meta.index = df_meta['Barcode']
    df_meta['celltype'] = df_meta.apply(lambda x: f"{x['Diagnosis']}_{x['Cell.Type']}", axis=1)
    sel_cells = random.sample(list(df_meta.index), len(df_meta.index)//10)
    dataset_AD_Control.adata = dataset_AD_Control.adata[sel_cells, :]
    dataset_AD_Control.adata.obs = df_meta.loc[sel_cells, :]

    dataset_AD_Control.quality_control(min_features=300, min_cells=10)
    dataset_AD_Control.select_genes(num_peak=120000)
    file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_AD_Control.add_promoter(file_gene_hg38)

    # PLAC
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    dataset_AD_Control.build_graph(path_hic)
    df_graph_PLAC = dataset_AD_Control.df_graph

    dataset_AD_Control.generate_data_list()
    list_graph_data = dataset_AD_Control.list_graph
    path_graph_input = os.path.join(path_sample_110, 'input_graph')
    os.system(f"rm -rf {path_graph_input}")
    os.mkdir(path_graph_input)
    dataset_atac_graph = ATACGraphDataset(path_graph_input, list_graph_data)

    dataset_AD_Control.find_neighbors()
    dataset_AD_Control.plot_umap()
    time_end = time()
    print(time_end - time_start)

    path_sample_110 = '/root/scATAC/ATAC_data/Alzheimer_Morabito/AD_Control_110'
    path_graph_input = os.path.join(path_sample_110, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:10000]
    test_dataset = dataset[10000:]

    device = torch.device("cuda:3" if torch.cuda.is_available() else "cpu")
    model = GCN(input_channels=dataset.num_node_features,
                output_channels=dataset.num_classes, hidden_channels=8,
                num_nodes=dataset_atac_graph[0].num_nodes).to(device)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    # train model
    time_start = time()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
    for epoch in range(1, 50):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, '
              f'Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')
        if test_acc > 0.97:
            break
    time_end = time()
    print(time_end - time_start)

    # explain model
    all_loader = DataLoader(dataset_atac_graph, batch_size=32, shuffle=True)
    list_dict = []
    method = 'ig'
    # i = 0
    for data in all_loader:
        data = data.to(device)
        target = data.y
        input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
        # dl = attr.DeepLift(model_forward)
        # mask = dl.attribute(input_mask, target=target,
        #                     additional_forward_args=(data, model))
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(
            input_mask, target=target, n_steps=50,
            additional_forward_args=(data, model),
            internal_batch_size=data.edge_index.shape[1])
        batch_size = len(torch.unique(data.batch))
        num_col = mask.shape[0]//batch_size
        mask = mask.view(batch_size, num_col)
        edge_mask = np.abs(mask.cpu().detach().numpy())
        edge_mask = edge_mask / np.max(edge_mask, axis=1)[:, np.newaxis]
        sub_edge_index = data.edge_index.cpu().numpy()
        col_edge = [(sub_edge_index[0, i], sub_edge_index[1, i]) for i in range(num_col)]
        list_dict.append(pd.DataFrame(edge_mask, columns=col_edge, index=data.cell))
        # i = i + 1
        # if i >= 10:
        #     break
    df_weight = pd.concat(list_dict)
    # df_weight = df_weight.dropna()
    # df_weight.index = dataset_ATAC.adata.obs.index

    # # save weight
    # file_weight = os.path.join(path_data_root, 'weight_atac.pkl')
    # with open(file_weight, 'wb') as w_pkl:
    #     str_pkl = pickle.dumps(df_weight)
    #     w_pkl.write(str_pkl)
    #
    # # read weight
    # path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_merge_HiC'
    # file_weight = os.path.join(path_data_root, 'weight_atac.pkl')
    # with open(file_weight, 'rb') as r_pkl:
    #     df_weight = pickle.loads(r_pkl.read())

    adata_edge = ad.AnnData(X=df_weight, obs=dataset_AD_Control.adata.obs.loc[df_weight.index, :])
    adata_edge.raw = adata_edge.copy()
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata_edge, n_top_genes=10000, flavor='seurat')
    adata = adata_edge[:, adata_edge.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.2)
    sc.pl.umap(adata, color=['nb_features', 'celltype'])
    sc.pl.umap(adata, color=['Diagnosis', 'Cell.Type'])

    # rank edge
    sc.tl.rank_genes_groups(adata_edge, groupby='celltype', method='wilcoxon', use_raw=True,
                            pts=True, tie_correct=True)
    array_names = adata_edge.uns['rank_genes_groups']['names']
    array_scores = adata_edge.uns['rank_genes_groups']['scores']
    array_pval = adata_edge.uns['rank_genes_groups']['pvals_adj']

    df_graph = dataset_AD_Control.df_graph
    peaks = dataset_AD_Control.array_peak
    celltypes = dataset_AD_Control.array_celltype
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

    # cell type
    sc.tl.rank_genes_groups(adata_edge, groupby='Cell.Type', method='wilcoxon', use_raw=True,
                            pts=True, tie_correct=True)
    array_names = adata_edge.uns['rank_genes_groups']['names']
    array_scores = adata_edge.uns['rank_genes_groups']['scores']
    array_pval = adata_edge.uns['rank_genes_groups']['pvals_adj']
    array_logfoldchanges = adata_edge.uns['rank_genes_groups']['logfoldchanges']
    peaks = dataset_AD_Control.array_peak
    celltypes = ['ASC', 'EX', 'INH', 'MG', 'ODC', 'OPC', 'PER.END']
    dict_cell_genes = {}
    dict_cell_genes_idx = {}
    for idx_cell in range(len(celltypes)):
        cell_name = celltypes[idx_cell]
        list_edge_idx = []
        list_edge_names = []
        list_edge_logfoldchanges = []
        list_edge_pval = []
        for idx_edge in range(array_names.shape[0]):
            edge_name = array_names[idx_edge][idx_cell]
            list_edge_idx.append(edge_name)
            list_edge_names.append((peaks[edge_name[0]], peaks[edge_name[1]]))
            list_edge_logfoldchanges.append(array_logfoldchanges[idx_edge][idx_cell])
            list_edge_pval.append(array_pval[idx_edge][idx_cell])
        dict_cell_genes[cell_name] = \
            [list_edge_names[i] for i in range(len(list_edge_names))
             if list_edge_logfoldchanges[i] > 1 and list_edge_pval[i] < 0.01]
        dict_cell_genes_idx[cell_name] = \
            [list_edge_idx[i] for i in range(len(list_edge_idx))
             if list_edge_logfoldchanges[i] > 1 and list_edge_pval[i] < 0.01]

    cell_genes = dict_cell_genes_idx['MG']
    cell_barcodes = adata_edge.obs.index[adata_edge.obs['Cell.Type'] == 'MG']
    adata_edge_cell = adata_edge[cell_barcodes, pd.Index(cell_genes)]
    sc.tl.rank_genes_groups(adata_edge_cell, groupby='Diagnosis', method='wilcoxon', use_raw=True,
                            pts=True, tie_correct=True)
    array_names = adata_edge_cell.uns['rank_genes_groups']['names']
    array_scores = adata_edge_cell.uns['rank_genes_groups']['scores']
    peaks = dataset_AD_Control.array_peak
    status = ['AD', 'Control']
    dict_cell_scores = {}
    for idx_status in range(len(status)):
        cell_name = status[idx_status]
        list_edge_names = []
        list_edge_scores = []
        for idx_edge in range(array_names.shape[0]):
            edge_name = array_names[idx_edge][idx_status]
            list_edge_names.append((peaks[edge_name[0]], peaks[edge_name[1]]))
            list_edge_scores.append(array_scores[idx_edge][idx_status])
        dict_cell_scores[cell_name] = pd.Series(list_edge_scores, index=list_edge_names)

    time_start = time()
    path_Control_3000 = '/root/scATAC/ATAC_data/Alzheimer_Morabito/Control_3000'
    dataset_Control_3000 = ATACDataset(data_root=path_Control_3000, raw_filename='counts.tsv')
    file_meta_tsv = os.path.join(path_Control_3000, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
    df_meta.index = df_meta['Barcode']
    dataset_Control_3000.adata.obs['celltype'] = df_meta.loc[dataset_Control_3000.adata.obs.index, 'Cell.Type']

    # dataset_Control_3000.quality_control(min_features=3000)
    file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_Control_3000.add_promoter(file_gene_hg38)

    dataset_Control_3000.find_neighbors()
    # dataset_Control_3000.plot_umap()
    time_end = time()
    print(time_end - time_start)

    # save data
    file_atac_Control_3000 = os.path.join(path_Control_3000, 'dataset_atac.pkl')
    with open(file_atac_Control_3000, 'wb') as w_pkl:
        str_pkl = pickle.dumps(dataset_Control_3000)
        w_pkl.write(str_pkl)

    # read data
    path_Control_3000 = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC'
    file_atac_Control_3000 = os.path.join(path_Control_3000, 'dataset_atac.pkl')
    with open(file_atac_Control_3000, 'rb') as r_pkl:
        dataset_Control_3000 = pickle.loads(r_pkl.read())

    time_start = time()
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    df_graph_Control_3000 = dataset_Control_3000.build_graph(path_hic)
    df_graph_Control_3000 = df_graph_Control_3000.drop_duplicates()
    torch.multiprocessing.set_sharing_strategy('file_system')
    list_graph_data_Control_3000, peaks_Control_3000, celltypes_Control_3000 = \
        generate_data_list(df_graph_Control_3000, dataset_Control_3000, 20)
    path_graph_input_Control_3000 = os.path.join(path_Control_3000, 'input_graph')
    dataset_atac_graph_Control_3000 = ATACGraphDataset(path_graph_input_Control_3000, list_graph_data_Control_3000)
    # dataset_atac_graph = ATACGraphDataset(path_graph_input)
    file_peaks_celltypes_Control_3000 = os.path.join(path_graph_input_Control_3000, 'peaks_celltypes')
    np.savez(file_peaks_celltypes_Control_3000, peaks=peaks_Control_3000, celltypes=celltypes_Control_3000)
    time_end = time()
    print(time_end - time_start)

    # build edge matrix
    list_edgr_attr = [pd.Series(one_data.edge_attr.numpy().reshape(-1),
                                index=[(one_data.edge_index.numpy()[0, i],
                                        one_data.edge_index.numpy()[1, i])
                                       for i in range(one_data.edge_index.numpy().shape[1])])
                      for one_data in dataset_atac_graph_Control_3000]
    df_edge = pd.concat(list_edgr_attr, axis=1)
    # df_edge = df_edge.applymap(lambda x: np.abs(x))
    adata_edge = ad.AnnData(X=df_edge.T, obs=dataset_Control_3000.adata.obs)
    adata_edge.raw = adata_edge.copy()
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata_edge, n_top_genes=10000, flavor='seurat')
    # adata_edge = adata_edge[:, adata_edge.var.highly_variable]
    sc.pp.scale(adata_edge, max_value=10)
    sc.tl.pca(adata_edge, svd_solver='arpack')
    sc.pp.neighbors(adata_edge, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata_edge)
    sc.pl.umap(adata_edge, color=['nb_features', 'celltype'])

    # rank edge
    sc.tl.rank_genes_groups(adata_edge, groupby='celltype', method='wilcoxon', use_raw=True)
    array_names = adata_edge.uns['rank_genes_groups']['names']
    array_scores = adata_edge.uns['rank_genes_groups']['scores']

    file_peaks_celltypes = os.path.join(path_graph_input_Control_3000, 'peaks_celltypes.npz')
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
