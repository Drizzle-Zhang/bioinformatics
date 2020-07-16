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


def sub_get_pairs(dict_in):
    sub_prob = dict_in['sub_prob']
    sub_diff = dict_in['sub_diff']
    gene = dict_in['gene']
    sub_combine = pd.concat([sub_prob, sub_diff], axis=1, sort=False)
    sub_combine.columns = ['prob', 'bool']
    genes2 = list(sub_combine.loc[
        (sub_combine['prob'] < 0.01) & sub_combine['bool'], :].index)

    return {gene: genes2}


def get_pairs_marker():
    # matrix of probability
    file_prob = '/home/zy/scRef/try_data/Prob_MCA_combined_by_cell.npy'
    file_genes = '/home/zy/scRef/try_data/Genes_MCA_combined_by_cell.txt'

    array_prob = np.load(file_prob)
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

    dict_cell_pairs = {}
    for cell in df_ref.columns:
        df_col = df_ref.loc[:, cell]
        array_col = np.array(df_col)
        mat_diff = array_col[:, np.newaxis] - array_col[np.newaxis, :]
        mat_diff = (mat_diff > 0)
        df_diff = pd.DataFrame(mat_diff, index=df_ref.index,
                               columns=df_ref.index)
        # select high expression genes
        cutoff_85 = np.percentile(array_col, 85)
        df_high = df_col.loc[df_col > cutoff_85]
        df_prob_high = df_prob.loc[df_high.index, ]

        # get pairs
        num_batch = 400
        num_cpu = 40
        list_input = []
        list_pairs = []
        for idx, gene in enumerate(df_prob_high.index):
            sub_diff = df_diff.loc[gene, :]
            sub_prob = df_prob_high.loc[gene, :]
            sub_dict = dict(sub_prob=sub_prob, sub_diff=sub_diff, gene=gene)
            list_input.append(sub_dict)
            # print(idx, gene)
            if len(list_input) % num_batch == 0:
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
        dict_cell_pairs[cell] = list_pairs

    return dict_cell_pairs, overlap_genes


def get_pairs_marker_2():
    # matrix of probability
    file_prob = '/home/zy/scRef/try_data/Prob_MCA_combined_by_cell.npy'
    file_genes = '/home/zy/scRef/try_data/Genes_MCA_combined_by_cell.txt'

    array_prob = np.load(file_prob)
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

    dict_cell_pairs = {}
    for cell in df_ref.columns:
        df_col = df_ref.loc[:, cell]
        array_col = np.array(df_col)
        mat_diff = array_col[:, np.newaxis] - array_col[np.newaxis, :]
        mat_diff = (mat_diff > 0)
        df_diff = pd.DataFrame(mat_diff, index=df_ref.index,
                               columns=df_ref.index)
        # select high expression genes
        cutoff_85 = np.percentile(array_col, 85)
        df_high = df_col.loc[df_col > cutoff_85]
        df_prob_high = df_prob.loc[df_high.index, ]

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
            list_pairs.extend({gene: genes2})
        dict_cell_pairs[cell] = list_pairs

    return dict_cell_pairs, overlap_genes


def sub_test(dict_in):
    sub_col = dict_in['sub_col']
    sub_tag = dict_in['sub_tag']
    scref_tag = sub_tag.loc['scRef.tag']
    pairs_marker = dict_cell_pairs[scref_tag]
    # num of 1 (marker)
    num_1_marker = 0
    num_sample = 0
    for i, dict_pairs in enumerate(pairs_marker):
        gene1 = list(dict_pairs.keys())[0]
        genes2 = list(dict_pairs.values())[0]
        if gene1 not in set_overlap_genes_sc:
            continue
        if not genes2:
            continue
        col_genes2 = sub_col.loc[genes2].dropna()
        array_diff = np.array(sub_col.loc[gene1] - col_genes2)
        if i == 0:
            array_bool = (array_diff > 0)
        else:
            array_bool = np.concatenate([array_bool, (array_diff > 0)])
        # num_sample = num_sample + array_diff.shape[0]
        # num_1_marker = num_1_marker + np.sum(array_diff > 0)

    num_batch = 20000
    num_permutation_marker = 20
    list_num_1_marker = []
    for i in range(num_permutation_marker):
        sub_array_bool = np.random.choice(array_bool, num_batch)
        list_num_1_marker.append(np.sum(sub_array_bool > 0))

    # num of 1 (permutation)
    num_permutation = 100
    num_genes = len(overlap_genes_sc)
    list_num_1 = []
    for i in range(num_permutation):
        list_idx1 = np.random.randint(0, num_genes, size=num_batch)
        list_idx2 = np.random.randint(0, num_genes, size=num_batch)
        array_diff = \
            np.array(
                sub_col.loc[[overlap_genes_sc[num] for num in list_idx1]]) - \
            np.array(sub_col.loc[[overlap_genes_sc[num] for num in list_idx2]])
        list_num_1.append(np.sum(array_diff > 0))

    _, pval = mannwhitneyu(list_num_1_marker, list_num_1,
                           alternative='greater')
    # ttest_1samp(np.array(list_num_1), num_1_marker)

    return {'cell_id': sub_col.name, 'p_value': pval}


if __name__ == '__main__':
    time_start = time()

    dict_cell_pairs, overlap_genes = get_pairs_marker()

    # Zeisel
    file_exp_sc = '/home/zy/scRef/try_data/summary/Zeisel_exp_sc_mat.txt'
    mat_sc = pd.read_csv(file_exp_sc, sep='\t', index_col=0)
    overlap_genes_sc = list(set(overlap_genes).intersection(set(mat_sc.index)))
    set_overlap_genes_sc = set(overlap_genes_sc)
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
    list_pval = pool.map(sub_test, list_input)
    pool.close()

    time_end = time()
    print(time_end - time_start)
