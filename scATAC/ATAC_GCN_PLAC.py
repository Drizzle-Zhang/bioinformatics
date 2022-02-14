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
import pickle
import torch
from torch_geometric.datasets import TUDataset
from torch_geometric.data import InMemoryDataset, Data, DataLoader
import torch.nn as nn
import torch.nn.functional as nn_func
import torch_geometric.nn as geo_nn
from captum.attr import Saliency, IntegratedGradients
import random
from scipy.stats import kstest


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

    def add_promoter(self, file_promoter, num_threads=3):
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
    print(i)
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

    list_idx = range(0, adata_ATAC.n_obs)
    # list_idx = range(1550, 1967)
    pool = Pool(num_threads)
    func_test = partial(one_test, all_genes, adata_ATAC, adata_merge_peak,
                        array_peak, array_celltype)
    result = pool.map(func_test, list_idx)
    pool.close()

    return result, array_peak, array_celltype


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
    def __init__(self, input_channels, output_channels, hidden_channels):
        super(GCN, self).__init__()
        torch.manual_seed(12345)
        self.lin1 = nn.Linear(input_channels, 8)
        self.conv1 = geo_nn.GraphConv(8, hidden_channels)
        self.conv2 = geo_nn.GraphConv(hidden_channels, hidden_channels)
        self.conv3 = geo_nn.GraphConv(hidden_channels, hidden_channels)
        self.lin2 = nn.Sequential(
            nn.Linear(hidden_channels, hidden_channels//2),
            nn.ReLU(),
            nn.Linear(hidden_channels//2, output_channels)
        )
        # self.softmax = nn.Softmax(dim=0)

    def forward(self, x, edge_index, edge_weight, batch):
        # 1. Obtain node embeddings
        x = self.lin1(x)
        x = x.relu()
        x = self.conv1(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv2(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv3(x, edge_index, edge_weight)

        # 2. Readout layer
        x = geo_nn.global_mean_pool(x, batch)  # [batch_size, hidden_channels]

        # 3. Apply a final classifier
        x = nn_func.dropout(x, p=0.5, training=self.training)
        x = self.lin2(x)
        # x = self.softmax(x)

        return x


def train(loader):
    model.train()
    for data in loader:  # Iterate in batches over the training dataset.
        data = data.to(device)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        # out = model(data.x, data.edge_index, data.batch)
        loss = criterion(out, data.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.


def test(loader):
    model.eval()
    correct = 0
    for data in loader:  # Iterate in batches over the training/test dataset.c
        data = data.to(device)
        # out = model(data.x, data.edge_index, data.batch)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        correct += int((pred == data.y).sum())  # Check against ground-truth labels.
    return correct / len(loader.dataset)  # Derive ratio of correct predictions.


def model_forward(edge_mask, data):
    batch = torch.zeros(data.x.shape[0], dtype=int).to(device)
    out = model(data.x, data.edge_index, edge_mask, batch)
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
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC'
    dataset_ATAC = ATACDataset(data_root=path_data_root, raw_filename='counts.tsv')
    file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
    dataset_ATAC.adata.obs['celltype'] = df_meta.loc[dataset_ATAC.adata.obs.index, 'celltype']

    dataset_ATAC.quality_control(min_features=3000)
    file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
    dataset_ATAC.add_promoter(file_gene_hg38)

    dataset_ATAC.find_neighbors()
    # dataset_ATAC.plot_umap()
    time_end = time()
    print(time_end - time_start)

    # save data
    file_atac_test = os.path.join(path_data_root, 'dataset_atac.pkl')
    with open(file_atac_test, 'wb') as w_pkl:
        str_pkl = pickle.dumps(dataset_ATAC)
        w_pkl.write(str_pkl)

    # read data
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC'
    file_atac_test = os.path.join(path_data_root, 'dataset_atac.pkl')
    with open(file_atac_test, 'rb') as r_pkl:
        dataset_ATAC = pickle.loads(r_pkl.read())

    time_start = time()
    path_hic = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
    df_graph = dataset_ATAC.build_graph(path_hic)
    df_graph = df_graph.drop_duplicates()
    torch.multiprocessing.set_sharing_strategy('file_system')
    list_graph_data, peaks, celltypes = generate_data_list(df_graph, dataset_ATAC, 25)
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input, list_graph_data)
    # dataset_atac_graph = ATACGraphDataset(path_graph_input)
    file_peaks_celltypes = os.path.join(path_graph_input, 'peaks_celltypes')
    np.savez(file_peaks_celltypes, peaks=peaks, celltypes=celltypes)
    time_end = time()
    print(time_end - time_start)

    # train model
    time_start = time()
    device = torch.device("cuda:2" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:1600]
    test_dataset = dataset[1600:]

    model = GCN(input_channels=dataset.num_node_features,
                output_channels=dataset.num_classes, hidden_channels=32)
    print(model)
    model.to(device)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.0025)
    for epoch in range(1, 201):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(1, 201):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')
    time_end = time()
    print(time_end - time_start)

    model.eval()
    correct = 0
    list_pred = []
    list_true = []
    for data in test_loader:  # Iterate in batches over the training/test dataset.c
        data = data.to(device)
        # out = model(data.x, data.edge_index, data.batch)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        list_pred.extend(pred.cpu().numpy())
        list_true.extend(data.y.cpu().numpy())
    df_res = pd.DataFrame(dict(pred=list_pred, ture=list_true))
    pd.crosstab(index=df_res['pred'], columns=df_res['ture'])

    data = random.choice([t for t in test_dataset if not t.y.item()]).to(device)
    method = 'ig'
    # method = 'saliency'
    edge_mask = explain(method, data, target=0)
    edge_mask_dict = aggregate_edge_directions(edge_mask, data)

    train_loader = DataLoader(dataset_atac_graph, batch_size=40, shuffle=True)
    train(train_loader)
    model.eval()
    correct_data = []
    for data in train_loader:  # Iterate in batches over the training/test dataset.c
        data = data.to(device)
        # out = model(data.x, data.edge_index, data.batch)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        df_res = pd.DataFrame(dict(pred=pred.cpu().numpy(), true=data.y.cpu().numpy()))
        correct_data.extend(
            data[list(df_res.index[df_res.apply(lambda x: x['pred'] == x['true'], axis=1)])])

    list_dict = []
    list_labels = []
    method = 'saliency'
    i = 0
    for data in correct_data:
        data = data.to(device)
        edge_mask = explain(method, data, target=data.y.to(device))
        edge_mask_dict = aggregate_edge_directions(edge_mask, data)
        list_dict.append(edge_mask_dict)
        list_labels.extend(data.y.cpu().numpy())
        i = i + 1
        if i > 10:
            break
    df_weight = pd.DataFrame(list_dict)

    adata_edge = ad.AnnData(X=df_weight,
                            obs=pd.DataFrame(dict(celltype=[str(label) for label in list_labels])))
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata_edge, n_top_genes=10000, flavor='seurat')
    adata = adata_edge[:, adata_edge.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['celltype'])

    # build edge matrix
    list_edgr_attr = [pd.Series(one_data.edge_attr.numpy().reshape(-1),
                                index=[(one_data.edge_index.numpy()[0, i],
                                        one_data.edge_index.numpy()[1, i])
                                       for i in range(one_data.edge_index.numpy().shape[1])])
                      for one_data in dataset_atac_graph]
    df_edge = pd.concat(list_edgr_attr, axis=1)
    # df_edge = df_edge.applymap(lambda x: np.abs(x))
    adata_edge = ad.AnnData(X=df_edge.T, obs=dataset_ATAC.adata.obs)
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

    # significant interaction
    dict_cell_scores = {}
    for idx_cell in range(len(celltypes)):
        cell_name = celltypes[idx_cell]
        list_edge_names = []
        list_edge_scores = []
        list_edge_pval = []
        for idx_edge in range(array_names.shape[0]):
            edge_name = array_names[idx_edge][idx_cell]
            list_edge_names.append((peaks[edge_name[0]], peaks[edge_name[1]]))
            list_edge_scores.append(array_scores[idx_edge][idx_cell])
            list_edge_pval.append(array_pval[idx_edge][idx_cell])
        dict_cell_scores[cell_name] = \
            pd.DataFrame({'scores':list_edge_scores, 'pval':list_edge_pval}, index=list_edge_names)

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

    for cell in list_cell:
        cell_pair = df_cell.index[df_cell[cell]]
        df_scores = pd.DataFrame(np.full(shape=(len(cell_pair), len(celltypes)), fill_value=0.0),
                                         index=cell_pair, columns=celltypes)
        for celltype in celltypes:
            all_score = dict_cell_scores[celltype]
            for sub_inter in cell_pair:
                df_scores.loc[df_scores.index == sub_inter, celltype] = \
                    all_score.loc[all_score.index == sub_inter, 'scores']
        file_scores = os.path.join(path_cell_interatome,
                                   f'{cell}/interactome_celltypes_scores.txt')
        df_scores.to_csv(file_scores, sep='\t')

    # peaks specificity
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC'
    dataset_ATAC = ATACDataset(data_root=path_data_root, raw_filename='counts.tsv')
    file_meta_tsv = os.path.join(path_data_root, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta_tsv, sep='\t', index_col=0)
    dataset_ATAC.adata.obs['celltype'] = df_meta.loc[dataset_ATAC.adata.obs.index, 'celltype']
    sc.pp.normalize_total(dataset_ATAC.adata)
    sc.pp.log1p(dataset_ATAC.adata)
    dataset_ATAC.adata.raw = dataset_ATAC.adata.copy()
    adata_peaks = dataset_ATAC.adata
    sc.tl.rank_genes_groups(adata_peaks, groupby='celltype', method='wilcoxon', use_raw=True,
                            pts=True, tie_correct=True)
    array_names = adata_peaks.uns['rank_genes_groups']['names']
    array_scores = adata_peaks.uns['rank_genes_groups']['scores']
    file_peaks = dataset_ATAC.file_peaks_sort
    list_peak_names = []
    list_peak_scores = []
    dict_peak_cell_scores = {}
    for idx_cell in range(len(celltypes)):
        cell_name = celltypes[idx_cell]
        list_peak_names = []
        list_peak_scores = []
        for idx_peak in range(array_names.shape[0]):
            list_peak_names.append(array_names[idx_peak][idx_cell])
            list_peak_scores.append(array_scores[idx_peak][idx_cell])
        dict_peak_cell_scores[cell_name] = pd.Series(list_peak_scores, index=list_peak_names)

    # dict promoter
    dict_promoter = {}
    file_peak_promoter = \
        '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain_PLAC/processed_files/peaks_promoter.txt'
    with open(file_peak_promoter, 'r') as w_pro:
        for line in w_pro:
            list_line = line.strip().split('\t')
            peak = list_line[3]
            gene = list_line[7].split('<-')[0]
            if gene != '.':
                dict_promoter[peak] = gene

    # AD sites
    path_AD = '/root/scATAC/GWAS/Jansen_NG_2019/'
    file_AD_sites = os.path.join(path_AD, 'AD_sumstats_Jansenetal_2019sept.txt')
    df_AD_sites = pd.read_csv(file_AD_sites, sep='\t')
    df_AD_signif = df_AD_sites.loc[df_AD_sites['P'] < 1*10**-5, :]
    path_AD_process = os.path.join(path_AD, 'process')
    if not os.path.exists(path_AD_process):
        os.mkdir(path_AD_process)
    file_AD_hg19 = os.path.join(path_AD_process, 'hg19.bed')
    with open(file_AD_hg19, 'w') as w_hg19:
        for i in df_AD_signif.index:
            w_hg19.write(f"chr{df_AD_signif.loc[i, 'CHR']}\t{df_AD_signif.loc[i, 'BP'] - 1}\t"
                         f"{df_AD_signif.loc[i, 'BP'] + 1}\t{df_AD_signif.loc[i, 'uniqID.a1a2']}\n")
    file_chain = '/root/tools/files_liftOver/hg19ToHg38.over.chain.gz'
    liftover = '/root/tools/liftOver'
    file_AD_hg38 = os.path.join(path_AD_process, 'hg38.bed')
    file_ummap = os.path.join(path_AD_process, 'unmap.bed')
    os.system(f"{liftover} {file_AD_hg19} {file_chain} {file_AD_hg38} {file_ummap}")
    # df_old = pd.read_csv(file_AD_hg38, sep='\t', header=None)
    # length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
    # df_old['length'] = length_old
    path_ATAC_AD = os.path.join(path_data_root, 'AD')
    if not os.path.exists(path_ATAC_AD):
        os.mkdir(path_ATAC_AD)
    file_intersect = os.path.join(path_ATAC_AD, 'peaks_AD.txt')
    os.system(f"bedtools intersect -a {file_peaks} -b {file_AD_hg38} -wao > {file_intersect}")
    df_peaks_AD = pd.read_csv(file_intersect, sep='\t', header=None)
    peaks_AD = df_peaks_AD.loc[df_peaks_AD.iloc[:, 4] != '.', 3].tolist()
    df_AD_peak = pd.DataFrame(np.full(shape=(len(peaks_AD), len(celltypes)), fill_value=0.0),
                              index=peaks_AD, columns=celltypes)
    for celltype in celltypes:
        for sub_peak in peaks_AD:
            all_score = dict_peak_cell_scores[celltype]
            df_AD_peak.loc[sub_peak, celltype] = all_score.loc[sub_peak]
    file_AD_scores = os.path.join(path_ATAC_AD, 'peaks_celltypes_scores.txt')
    df_AD_peak.to_csv(file_AD_scores, sep='\t')
    # AD interactome
    list_AD_merge_peaks = []
    for peak in peaks_AD:
        if peak in dict_promoter.keys():
            list_AD_merge_peaks.append(dict_promoter[peak])
            continue
        list_AD_merge_peaks.append(peak)
    list_interatome = []
    for pair in df_cell.index:
        pair_1 = pair[0]
        pair_2 = pair[1]
        if pair_2[0:3] == 'chr':
            if pair_2 in list_AD_merge_peaks:
                list_interatome.append(pair)
        else:
            if (pair_1 in list_AD_merge_peaks) | (pair_2 in list_AD_merge_peaks):
                list_interatome.append(pair)
    df_AD_interactome = pd.DataFrame(np.full(shape=(len(list_interatome), len(celltypes)),
                                             fill_value=0.0),
                                     index=list_interatome, columns=celltypes)
    for celltype in celltypes:
        all_score = dict_cell_scores[celltype]
        for sub_inter in set(list_interatome):
            df_AD_interactome.loc[df_AD_interactome.index == sub_inter, celltype] = \
                all_score.loc[all_score.index == sub_inter, 'scores']
    file_interactome_AD_scores = os.path.join(path_ATAC_AD, 'interactome_celltype_scores.txt')
    df_AD_interactome.to_csv(file_interactome_AD_scores, sep='\t')

    # SCZ
    path_SCZ = '/root/scATAC/GWAS/Antonio_NG_2018/'
    file_SCZ_sites = os.path.join(path_SCZ, 'clozuk_pgc2.meta.sumstats.txt')
    df_SCZ_sites = pd.read_csv(file_SCZ_sites, sep='\t')
    df_SCZ_signif = df_SCZ_sites.loc[df_SCZ_sites['P'] < 1*10**-5, :]
    path_SCZ_process = os.path.join(path_SCZ, 'process')
    if not os.path.exists(path_SCZ_process):
        os.mkdir(path_SCZ_process)
    file_SCZ_hg19 = os.path.join(path_SCZ_process, 'hg19.bed')
    with open(file_SCZ_hg19, 'w') as w_hg19:
        for i in df_SCZ_signif.index:
            w_hg19.write(f"chr{df_SCZ_signif.loc[i, 'CHR']}\t{df_SCZ_signif.loc[i, 'BP'] - 1}\t"
                         f"{df_SCZ_signif.loc[i, 'BP'] + 1}\t{df_SCZ_signif.loc[i, 'SNP']}\n")
    file_chain = '/root/tools/files_liftOver/hg19ToHg38.over.chain.gz'
    liftover = '/root/tools/liftOver'
    file_SCZ_hg38 = os.path.join(path_SCZ_process, 'hg38.bed')
    file_ummap = os.path.join(path_SCZ_process, 'unmap.bed')
    os.system(f"{liftover} {file_SCZ_hg19} {file_chain} {file_SCZ_hg38} {file_ummap}")
    # df_old = pd.read_csv(file_SCZ_hg38, sep='\t', header=None)
    # length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
    # df_old['length'] = length_old
    path_ATAC_SCZ = os.path.join(path_data_root, 'SCZ')
    if not os.path.exists(path_ATAC_SCZ):
        os.mkdir(path_ATAC_SCZ)
    file_intersect = os.path.join(path_ATAC_SCZ, 'peaks_SCZ.txt')
    os.system(f"bedtools intersect -a {file_peaks} -b {file_SCZ_hg38} -wao > {file_intersect}")
    df_peaks_SCZ = pd.read_csv(file_intersect, sep='\t', header=None)
    peaks_SCZ = df_peaks_SCZ.loc[df_peaks_SCZ.iloc[:, 4] != '.', 3].tolist()
    df_SCZ_peak = pd.DataFrame(np.full(shape=(len(peaks_SCZ), len(celltypes)), fill_value=0.0),
                               index=peaks_SCZ, columns=celltypes)
    for celltype in celltypes:
        for sub_peak in peaks_SCZ:
            all_score = dict_peak_cell_scores[celltype]
            df_SCZ_peak.loc[sub_peak, celltype] = all_score.loc[sub_peak]
    file_SCZ_scores = os.path.join(path_ATAC_SCZ, 'peaks_celltypes_scores.txt')
    df_SCZ_peak.to_csv(file_SCZ_scores, sep='\t')
    # SCZ interactome
    list_SCZ_merge_peaks = []
    for peak in peaks_SCZ:
        if peak in dict_promoter.keys():
            list_SCZ_merge_peaks.append(dict_promoter[peak])
            continue
        list_SCZ_merge_peaks.append(peak)
    list_interatome = []
    for pair in df_cell.index:
        pair_1 = pair[0]
        pair_2 = pair[1]
        if pair_2[0:3] == 'chr':
            if pair_2 in list_SCZ_merge_peaks:
                list_interatome.append(pair)
        # else:
        #     if (pair_1 in list_SCZ_merge_peaks) | (pair_2 in list_SCZ_merge_peaks):
        #         list_interatome.append(pair)
    df_SCZ_interactome = pd.DataFrame(np.full(shape=(len(list_interatome), len(celltypes)),
                                             fill_value=0.0),
                                     index=list_interatome, columns=celltypes)
    for celltype in celltypes:
        all_score = dict_cell_scores[celltype]
        for sub_inter in set(list_interatome):
            df_SCZ_interactome.loc[df_SCZ_interactome.index == sub_inter, celltype] = \
                all_score.loc[all_score.index == sub_inter, 'scores']
    file_interactome_SCZ_scores = os.path.join(path_ATAC_SCZ, 'interactome_celltype_scores.txt')
    df_SCZ_interactome.to_csv(file_interactome_SCZ_scores, sep='\t')

    # ks test
    # dict_peaks_disease = {'AD': peaks_AD, 'SCZ': peaks_SCZ}
    # list_disease = ['AD', 'SCZ']
    # df_pval_disease = pd.DataFrame(np.full(shape=(len(list_disease), len(celltypes)), fill_value=1),
    #                                index=list_disease, columns=celltypes)
    # for disease in list_disease:
    #     peaks_disease = dict_peaks_disease[disease]
    #     for celltype in celltypes:
    #         all_score = dict_peak_cell_scores[celltype]
    #         cell_score = all_score[peaks_disease]
    #         df_pval_disease.loc[disease, celltype] = \
    #             kstest(np.array(cell_score), np.array(all_score), alternative='less')[1]

