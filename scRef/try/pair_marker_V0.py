#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: pair_marker.py
# @time: 7/14/20 9:54 PM

from time import time
import pandas as pd
import numpy as np
from multiprocessing import Pool


def mat_bool_cell(col_in):
    mat_diff = col_in[:, np.newaxis] - col_in[np.newaxis, :]
    mat_diff = (mat_diff > 0)

    return mat_diff

########## MCA
time_start = time()
file_MCA = '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse.txt'
df_MCA = pd.read_csv(file_MCA, sep='\t', index_col=0)
genes_MCA = df_MCA.index
list_vec_cell = [df_MCA[cell] for cell in df_MCA.columns]
# list_vec_cell = list_vec_cell[:50]

num_cpu = 5
num_batch = 10
num_cell = len(list_vec_cell)
for i in range(num_cell // num_batch + 1):
    if i == num_cell // num_batch:
        sub_list = list_vec_cell[10 * i:]
    else:
        sub_list = list_vec_cell[10 * i:10 * (i + 1)]
    if not sub_list:
        continue
    pool = Pool(num_cpu)
    list_mat_cell = pool.map(mat_bool_cell, sub_list)
    pool.close()
    array_3d = np.dstack(list_mat_cell)
    if i == 0:
        array_1 = np.sum(array_3d, axis=2, dtype=np.int16)
    else:
        array_1 = array_1 + np.sum(array_3d, axis=2, dtype=np.int16)
    list_mat_cell = []
    array_3d = []

array_prob = (array_1 / num_cell).astype(np.float16)
series_genes = pd.Series(genes_MCA)
# df_prob = pd.DataFrame(array_prob, index=genes_MCA, columns=genes_MCA,
#                        dtype=np.float16)
file_prob = '/home/zy/scRef/try_data/Prob_MCA_combined_by_cell.npy'
file_genes = '/home/zy/scRef/try_data/Genes_MCA_combined_by_cell.txt'
np.save(file_prob, array_prob)
series_genes.to_csv(file_genes, sep='\t', index=False, header=False)

time_end = time()
print(time_end - time_start)


time_start = time()
file_ref = '/home/zy/scRef/try_data/scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
df_ref = pd.read_csv(file_ref, sep='\t', index_col=0)
df_col_1 = df_ref.loc[:, 'Astrocytes']
col_1 = np.array(df_col_1)
mat_diff = col_1[:, np.newaxis] - col_1[np.newaxis, :]
mat_diff = (mat_diff > 0)
# mat_diff = np.where(mat_diff > 0, 1, 0).astype(np.int8)
df_diff = pd.DataFrame(mat_diff, index=df_ref.index, columns=df_ref.index)

# matrix of probability
df_prob = pd.read_csv(file_prob, sep='\t', index_col=0)
# overlap genes
overlap_genes = list(set(df_prob.index).intersection(set(df_diff.index)))
df_diff = df_diff.loc[overlap_genes, overlap_genes]
# select high expression genes
cutoff_85 = np.percentile(col_1, 85)
df_high = df_col_1.loc[df_col_1 > cutoff_85]

df_prob_high = df_prob.loc[df_high.index, overlap_genes]
df_prob_high = df_prob_high.dropna()


def sub_get_pairs(dict_in):
    sub_prob = dict_in['sub_prob']
    sub_diff = dict_in['sub_diff']
    gene = dict_in['gene']
    sub_combine = pd.concat([sub_prob, sub_diff], axis=1, sort=False)
    sub_combine.columns = ['prob', 'bool']
    genes2 = list(sub_combine.loc[
        (sub_combine['prob'] < 0.01) & sub_combine['bool'], :].index)

    return {gene: genes2}


num_batch = 400
num_cpu = 40
list_input = []
list_pairs = []
for idx, gene in enumerate(df_prob_high.index):
    sub_diff = df_diff.loc[gene, :]
    sub_prob = df_prob_high.loc[gene, :]
    sub_dict = dict(sub_prob=sub_prob, sub_diff=sub_diff, gene=gene)
    list_input.append(sub_dict)
    print(idx, gene)
    if len(list_input) % 200 == 0:
        pool = Pool(num_cpu)
        pairs_batch = pool.map(sub_get_pairs, list_input)
        pool.close()
        list_input = []
        list_pairs.extend(pairs_batch)

pool = Pool(num_cpu)
pairs_batch = pool.map(sub_get_pairs, list_input)
pool.close()
list_input = []
list_pairs.extend(pairs_batch)


def sub_test():

    return


# Zeisel
file_exp_sc = '/home/zy/scRef/try_data/summary/Zeisel_exp_sc_mat.txt'
mat_sc = pd.read_csv(file_exp_sc, sep='\t', index_col=0)
overlap_genes_sc = list(set(overlap_genes).intersection(set(mat_sc.index)))
mat_sc = mat_sc.loc[overlap_genes_sc, :]
mat_sc.columns = [f"x{col.replace('_', '.')}" for col in mat_sc.columns]


time_end = time()
print(time_end - time_start)

