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
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
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
adata_gene = adata[:, [one_peak for one_peak in adata.var.index if one_peak in all_peaks]]
df_gene_peak = pd.DataFrame(adata_gene.X, index=adata_gene.obs.index,
                            columns=adata_gene.var.index)


all_cols = df_gene_peak.columns


def sum_peaks(df_gene_peak, dict_promoter, all_cols, gene):
    gene_peaks = dict_promoter[gene]
    gene_peaks = set(gene_peaks).intersection(set(all_cols))
    if len(gene_peaks) == 1:
        df_sum = df_gene_peak.loc[:, gene_peaks].squeeze('columns')
    elif len(gene_peaks) > 1:
        df_sum = np.sum(df_gene_peak.loc[:, gene_peaks], axis=1)
    else:
        return
    df_sum.name = gene
    return df_sum


pool = Pool(30)
func_sum = partial(sum_peaks, df_gene_peak, dict_promoter, all_cols)
result = pool.map(func_sum, all_genes)
pool.close()
result = [one_df for one_df in result if one_df is not None]
df_gene = pd.concat(result, axis=1)



# process PO file
path_hic = '/root/scATAC/pcHi-C'
file_po = os.path.join(path_hic, 'PO.txt')
file_po_bed = os.path.join(path_hic, 'PO.bed')
fmt_po = "{chrom}\t{start}\t{end}\t{peak_id}\t{gene}\t{tissue}\n"
with open(file_po_bed, 'w') as w_bed:
    with open(file_po, 'r') as r_po:
        for line in r_po:
            list_line = line.strip().split('\t')
            chrom = list_line[1].strip().split('.')[0]
            start = list_line[1].strip().split('.')[1]
            end = list_line[1].strip().split('.')[2]
            peak_id = f"{chrom}:{start}-{end}"
            gene = list_line[0]
            tissue = list_line[2]
            w_bed.write(fmt_po.format(**locals()))

# hg19 to hg38
file_chain = '/root/tools/files_liftOver/hg19ToHg38.over.chain.gz'
liftover = '/root/tools/liftOver'

path_hic_tissue = '/root/scATAC/pcHi-C/PO_by_tissue'
df_po = pd.read_csv(file_po_bed, sep='\t', header=None)
length_o = df_po.iloc[:, 2] - df_po.iloc[:, 1]
df_po['length'] = length_o
df_po = df_po.loc[df_po['length'] < 20000, :]
all_tissue = set(df_po.iloc[:, 5].tolist())
for one_tissue in all_tissue:
    path_tissue = os.path.join(path_hic_tissue, one_tissue.replace(' ', '_'))
    if not os.path.exists(path_tissue):
        os.mkdir(path_tissue)
    file_hg19 = os.path.join(path_tissue, 'hg19.bed')
    df_tissue = df_po.loc[df_po.iloc[:, 5] == one_tissue, :]
    df_tissue.loc[:, 'interaction_id'] = df_tissue.apply(lambda x: f"{x.iloc[3]}_{x.iloc[4]}", axis=1)
    df_tissue.to_csv(file_hg19, sep='\t', index=None, header=None)
    file_hg38 = os.path.join(path_tissue, 'hg38.bed')
    file_prefix = file_hg19 + '.prefix'
    file_suffix = file_hg19 + '.suffix'
    file_hg38_prefix = file_hg38 + '.prefix'
    file_hg38_format = file_hg38 + '.format'
    file_ummap = os.path.join(path_tissue, 'unmap.bed')
    os.system(f"cut -f 1,2,3,8 {file_hg19} > {file_prefix}")
    os.system(f"cut -f 4,5,6,8 {file_hg19} > {file_suffix}")
    os.system(f"{liftover} {file_prefix} {file_chain} "
              f"{file_hg38_prefix} {file_ummap}")
    dict_peak_score = defaultdict(list)
    with open(file_suffix, 'r') as r_f:
        for line in r_f:
            list_line = line.strip().split('\t')
            dict_peak_score[list_line[3]].append(list_line[0:3])
    with open(file_hg38_format, 'w') as w_f:
        fmt = "{chrom}\t{start}\t{end}\t{interaction_id}\t{peak_id}\t{gene}\t{tissue}\n"
        with open(file_hg38_prefix, 'r') as r_hg38:
            for line in r_hg38:
                list_line = line.strip().split('\t')
                list_suffix = dict_peak_score[list_line[3]][0]
                dict_hg38 = dict(
                    chrom=list_line[0], start=list_line[1], end=list_line[2],
                    interaction_id=list_line[3],
                    peak_id=f"{list_line[0]}:{list_line[1]}-{list_line[2]}",
                    gene=list_suffix[1], tissue=list_suffix[2]
                )
                w_f.write(fmt.format(**dict_hg38))

    df_old = pd.read_csv(file_prefix, sep='\t', header=None)
    length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
    df_old['length'] = length_old
    df_bed = pd.read_csv(file_hg38_format, sep='\t', header=None)
    length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
    # df_bed['length'] = length
    df_bed = df_bed.loc[length < 20000, :]
    df_bed = df_bed.drop_duplicates()
    df_bed.to_csv(file_hg38, sep='\t', index=None, header=None)

# interactions for each tissue
path_interaction = '/root/scATAC/pcHi-C/Interactions_by_tissue'
file_pp = os.path.join(path_hic, 'PP.txt')
df_pp = pd.read_csv(file_pp, sep='\t', header=None)
all_tissue = set(df_pp.iloc[:, 2].tolist())
for one_tissue in all_tissue:
    path_po_tissue = os.path.join(path_hic_tissue, one_tissue.replace(' ', '_'))
    path_tissue = os.path.join(path_interaction, one_tissue.replace(' ', '_'))
    if not os.path.exists(path_tissue):
        os.mkdir(path_tissue)
    file_PP = os.path.join(path_tissue, 'PP.txt')
    df_pp_tissue = df_pp.loc[df_pp.iloc[:, 2] == one_tissue, [0, 1, 2]]
    with open(file_PP, 'w') as w_pp:
        for sub_dict in df_pp_tissue.to_dict('records'):
            set_gene1 = sub_dict[0].strip().split(';')
            set_gene2 = sub_dict[1].strip().split(';')
            for gene1 in set_gene1:
                for gene2 in set_gene2:
                    w_pp.write(f"{gene1}\t{gene2}\n")

    file_po_tissue = os.path.join(path_po_tissue, 'hg38.bed')
    file_PO = os.path.join(path_tissue, 'PO.txt')
    with open(file_PO, 'w') as w_po:
        with open(file_po_tissue, 'r') as r_po:
            for line in r_po:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                peak_id = list_line[4]
                set_gene = list_line[5].strip().split(';')
                for gene in set_gene:
                    w_po.write(f"{chrom}\t{start}\t{end}\t{peak_id}\t{gene}\n")


