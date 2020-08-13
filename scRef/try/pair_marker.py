#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: pair_marker.py
# @time: 7/14/20 9:54 PM

from time import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
from scipy.stats import mannwhitneyu, ttest_1samp
from functools import partial


def mat_bool_cell(col_in):
    mat_diff = col_in[:, np.newaxis] - col_in[np.newaxis, :]
    mat_diff = (mat_diff > 0)

    return mat_diff


def prepare_prob_mat():
    # MCA
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

    return


def sub_get_pairs(df_prob, df_ref, cell):
    df_col = df_ref.loc[:, cell]
    array_col = np.array(df_col)
    mat_diff = array_col[:, np.newaxis] - array_col[np.newaxis, :]
    mat_diff = np.triu(mat_diff)
    mat_diff = (mat_diff > 0)
    df_diff = pd.DataFrame(mat_diff, index=df_ref.index,
                           columns=df_ref.index)
    # select high expression genes
    cutoff_90 = np.percentile(array_col, 90)
    df_high = df_col.loc[df_col > cutoff_90]
    df_prob_high = df_prob.loc[df_high.index, :]

    # get pairs
    list_pairs = []
    for idx, gene in enumerate(df_prob_high.index):
        sub_diff = df_diff.loc[gene, :]
        sub_prob = df_prob_high.loc[gene, :]
        sub_combine = pd.concat([sub_prob, sub_diff], axis=1, sort=False)
        sub_combine.columns = ['prob', 'bool']
        genes2 = list(sub_combine.loc[
                      (sub_combine['prob'] < 0.01) & sub_combine['bool'],
                      :].index)
        list_pairs.append({gene: genes2})

    return list_pairs


def get_pairs_marker():
    # matrix of probability
    file_prob = '/home/zy/scRef/try_data/Prob_MCA_combined_by_cell.npy'
    file_genes = '/home/zy/scRef/try_data/Genes_MCA_combined_by_cell.txt'

    array_prob = np.load(file_prob)
    array_prob = np.triu(array_prob)
    series_genes = pd.read_csv(file_genes, sep='\t', header=None)
    genes_MCA = series_genes.iloc[:, 0].tolist()
    df_prob = pd.DataFrame(array_prob, index=genes_MCA, columns=genes_MCA,
                           dtype=np.float16)
    array_prob = []
    series_genes = []

    # reference
    file_ref = '/home/zy/scRef/try_data/scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
    df_ref = pd.read_csv(file_ref, sep='\t', index_col=0)

    # overlap genes
    overlap_genes = list(set(df_prob.index).intersection(set(df_ref.index)))
    df_ref = df_ref.loc[overlap_genes, ]
    df_prob = df_prob.loc[overlap_genes, overlap_genes]

    pool = Pool(10)
    func_get_pairs = partial(sub_get_pairs, df_prob, df_ref)
    list_out = pool.map(func_get_pairs, df_ref.columns)
    pool.close()

    dict_cell_pairs = {}
    for i, cell in enumerate(df_ref.columns):
        dict_cell_pairs[cell] = list_out[i]

    return dict_cell_pairs, overlap_genes


def unfold_dict(dict_ref, set_genes, cell_ref):
    pairs_marker = dict_ref[cell_ref]
    list_gene1 = []
    list_gene2 = []
    for i, dict_pairs in enumerate(pairs_marker):
        gene1 = list(dict_pairs.keys())[0]
        genes2 = list(dict_pairs.values())[0]
        if gene1 not in set_genes:
            continue
        if not genes2:
            continue
        sub_genes2 = list(set(genes2).intersection(set_genes))
        list_gene2.extend(sub_genes2)
        sub_gene1 = [gene1] * len(sub_genes2)
        list_gene1.extend(sub_gene1)

    dict_cell = dict(gene1=list_gene1, gene2=list_gene2)

    return dict_cell


def sub_background(dict_in):
    sub_col = dict_in['sub_col']

    cells = dict_cell_pairs.keys()
    dict_out = {}
    for cell in cells:
        pairs_marker = dict_cell_ref[cell]
        num_batch = 10000
        num_permutation_marker = 10
        list_gene1 = pairs_marker['gene1']
        list_gene2 = pairs_marker['gene2']
        array_idx_marker = np.arange(len(list_gene2))
        list_num_1_marker = []
        for i in range(num_permutation_marker):
            sub_array_idx = np.random.choice(array_idx_marker, num_batch)
            sub_gene1 = [list_gene1[j] for j in sub_array_idx]
            sub_gene2 = [list_gene2[j] for j in sub_array_idx]
            array_diff = \
                np.array(sub_col.loc[sub_gene1]) - np.array(
                    sub_col.loc[sub_gene2])
            list_num_1_marker.append(np.sum(array_diff > 0))
        dict_out[cell] = sum(list_num_1_marker) / num_permutation_marker

    return dict_out


def sub_test(dict_in):
    sub_col = dict_in['sub_col']
    sub_tag = dict_in['sub_tag']
    scref_tag = sub_tag.loc['scRef.tag']
    pairs_marker = dict_cell_ref[scref_tag]
    cell_id = sub_col.name
    # num of background
    cell_background = (df_background.sample(100)).loc[:, scref_tag]

    # num of 1 (marker)
    num_batch = 10000
    num_permutation_marker = 100
    list_gene1 = pairs_marker['gene1']
    list_gene2 = pairs_marker['gene2']
    array_idx_marker = np.arange(len(list_gene2))
    list_num_1_marker = []
    for i in range(num_permutation_marker):
        sub_array_idx = np.random.choice(array_idx_marker, num_batch)
        sub_gene1 = [list_gene1[j] for j in sub_array_idx]
        sub_gene2 = [list_gene2[j] for j in sub_array_idx]
        array_diff = \
            np.array(sub_col.loc[sub_gene1]) - np.array(sub_col.loc[sub_gene2])
        list_num_1_marker.append(np.sum(array_diff > 0))

    _, pval = mannwhitneyu(list_num_1_marker, cell_background,
                           alternative='greater')
    # ttest_1samp(np.array(list_num_1), num_1_marker)

    return {'cell_id': cell_id, 'p_value': pval,
            'fold_change': np.mean(list_num_1_marker)/np.mean(cell_background)}


if __name__ == '__main__':
    time_start = time()

    dict_cell_pairs, overlap_genes = get_pairs_marker()

    # MCA samples
    file_samples_HCA = \
        '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/sample.txt'
    df_samples = pd.read_csv(file_samples_HCA, sep='\t', index_col=0)

    # Zeisel
    file_exp_sc = '/home/zy/scRef/try_data/summary/Zeisel_exp_sc_mat.txt'
    mat_sc = pd.read_csv(file_exp_sc, sep='\t', index_col=0)
    mat_sc = mat_sc.astype(np.int32)
    overlap_genes_sc = list((set(overlap_genes).intersection(
        set(mat_sc.index))).intersection(set(df_samples.index)))
    set_overlap_genes_sc = set(overlap_genes_sc)
    # simplify dict
    cells_ref = dict_cell_pairs.keys()
    dict_cell_ref = {}
    for cell_ref in cells_ref:
        dict_cell = unfold_dict(
            dict_cell_pairs, set_overlap_genes_sc, cell_ref)
        dict_cell_ref[cell_ref] = dict_cell

    # background distribution
    df_samples_1 = (df_samples.loc[overlap_genes_sc, :]).sample(1000, axis=1)
    list_input = []
    for col in df_samples.columns:
        sub_dict = dict(sub_col=df_samples.loc[:, col])
        list_input.append(sub_dict)
    pool = Pool(40)
    list_background = pool.map(sub_background, list_input)
    pool.close()
    df_background = pd.DataFrame(list_background, index=df_samples.columns)

    mat_sc = mat_sc.loc[overlap_genes_sc, :]
    mat_sc.columns = [f"X{col.replace('_', '.')}" for col in mat_sc.columns]
    file_tags = '/home/zy/scRef/try_data/tags_Zeisel1.txt'
    df_tag = pd.read_csv(file_tags, sep='\t', index_col=0)
    col = mat_sc.columns[0]
    sub_dict = dict(sub_col=mat_sc.loc[:, col], sub_tag=df_tag.loc[col, :])
    list_input = []
    for col in mat_sc.columns:
        sub_dict = dict(sub_col=mat_sc.loc[:, col], sub_tag=df_tag.loc[col, :])
        list_input.append(sub_dict)
    pool = Pool(40)
    # func_test = partial(sub_test, dict_cell_ref)
    list_pval = pool.map(sub_test, list_input[-40:])
    pool.close()

    # test
    df_samples_2 = (df_samples.loc[overlap_genes_sc, :]).sample(100, axis=1)
    list_input = []
    for col in df_samples_2.columns:
        sub_dict = dict(sub_col=df_samples_2.loc[:, col],
                        sub_tag=df_tag.iloc[0, :])
        list_input.append(sub_dict)
    pool = Pool(40)
    # func_test = partial(sub_test, dict_cell_ref)
    list_pval = pool.map(sub_test, list_input)
    pool.close()

    df_pval = pd.DataFrame(list_pval)
    df_pval.index = df_pval['cell_id']
    df_pval = df_pval.drop('cell_id', axis=1)
    df_meta = pd.concat([df_tag, df_pval], axis=1, join='inner')
    file_res = '/home/zy/scRef/try_data/tags_Zeisel.res1.txt'
    df_meta.to_csv(file_res, sep='\t')

    time_end = time()
    print(time_end - time_start)
