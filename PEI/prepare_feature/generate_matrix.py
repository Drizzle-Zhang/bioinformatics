#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_matrix.py
# @time: 2020/5/3 22:55

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool
from functools import partial


def sub_dhs_term(df_all, path_tmp, dict_in):
    file_index = dict_in['file_index']
    term = dict_in['term']
    str_term = dict_in['str_term']
    df_term = pd.read_csv(file_index, sep='\t', header=None,
                          index_col=1, usecols=[1, 3], names=[term, ''])
    # print(np.min(df_term))
    # print(np.max(df_term))
    df_merge = pd.merge(df_all, df_term, how='outer',
                        left_index=True, right_index=True)
    df_merge = df_merge.fillna(np.min(df_term)[0] - 1)
    file_tmp = os.path.join(path_tmp, str_term + '.txt')
    if df_merge.shape[0] != df_all.shape[0]:
        print(term)
    df_merge.to_csv(file_tmp, sep='\t', index=None, columns=[term])

    return file_tmp


def dhs_matrix():
    df_all = pd.read_csv(
        file_all_index, sep='\t', header=None, usecols=[3], names=['dhs_id'])
    path_tmp = os.path.join(path_matrix, 'tmp_DHS')
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    file_dhs_id = os.path.join(path_tmp, 'dhs_id.txt')
    df_all.to_csv(file_dhs_id, sep='\t', index=None, columns=['dhs_id'])

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_dhs_cell, 'meta.reference.tsv'), sep='\t')
    list_tmp = [file_dhs_id]
    list_input = []
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_dhs_cell, str_term)
        file_index = os.path.join(path_term, 'index.txt')
        list_input.append({'file_index': file_index, 'term': term,
                           'str_term': str_term})

    # tissue
    df_meta_tissue = pd.read_csv(
        os.path.join(path_dhs_tissue_stan, 'meta.reference.tsv'), sep='\t')
    for i in range(df_meta_tissue.shape[0]):
        life_organ = df_meta_tissue.loc[i, 'Biosample life_organ']
        term = df_meta_tissue.loc[i, 'Biosample term name']
        str_life_organ = life_organ.replace(' ', '_')
        str_term = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        path_term = os.path.join(
            path_dhs_tissue_stan, f"{str_life_organ}/{str_term}")
        file_index = os.path.join(path_term, 'index.txt')
        list_input.append({'file_index': file_index,
                           'term': f"{life_organ} | {term}",
                           'str_term': f"{str_life_organ}+{str_term}"})

    pool = Pool(num_cpu)
    func_dhs = partial(sub_dhs_term, df_all, path_tmp)
    list_res = pool.map(func_dhs, list_input)
    pool.close()

    list_tmp.extend(list_res)
    file_matrix_dhs = os.path.join(path_matrix, 'DHS_matrix.txt')
    os.system(f"paste {' '.join(list_tmp)} > {file_matrix_dhs}")

    return


def h3k4me3_matrix():
    os.system(
        f"bedtools intersect -wa -wb -a {file_all_index} -b {file_promoter} | "
        f"cut -f 4,5,9 > {file_dhs_promoter}")
    df_dhs_promoter = pd.read_csv(file_dhs_promoter, sep='\t', header=None)
    df_dhs_promoter['gene'] = df_dhs_promoter[2].apply(
        lambda x: x.split('<-')[0])
    print(df_dhs_promoter.shape[0])
    df_dhs_promoter = df_dhs_promoter.drop_duplicates()
    print(df_dhs_promoter.shape[0])
    df_dhs_promoter.to_csv(
        file_dhs_promoter, sep='\t', index=None, header=None)

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_h3k4me3_cell, 'meta.reference.tsv'), sep='\t')


    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    file_all_index = \
        '/local/zy/PEI/mid_data/database_feature/DHS_index/all_index.txt'
    path_matrix = '/local/zy/PEI/mid_data/database_feature/matrix'
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'

    # DHS
    path_dhs_cell = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_stan = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_standard'
    dhs_matrix()

    # H3K4me3
    path_h3k4me3_cell = '/local/zy/PEI/mid_data/cell_line/DHS/reference_map'
    path_h3k4me3_tissue = '/local/zy/PEI/mid_data/tissue/DHS/reference_map'
    file_dhs_promoter = \
        '/local/zy/PEI/mid_data/database_feature/DHS_index/promoter_index.txt'

    time_end = time()
    print(time_end - time_start)
