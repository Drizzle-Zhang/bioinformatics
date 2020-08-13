#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: HCA.py
# @time: 7/3/20 11:16 AM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
import random
import codecs
import chardet

# path of input
path_HCA = '/home/disk/scRef/HumanAtlas_SingleCell_Han2020'
path_count = os.path.join(path_HCA, 'new')
path_annotation = os.path.join(path_HCA, 'annotation_rmbatch_data_revised417')
# path_annotation_utf8 = os.path.join(path_HCA, 'annotation')
path_tissue = os.path.join(path_HCA, 'HCA_by_tissue')
path_combine = os.path.join(path_HCA, 'combinedHCA')

time_start = time()
# ANSI to utf8
# files_ansi = os.listdir(path_annotation)
# for file in files_ansi:
#     file_ansi = os.path.join(path_annotation, file)
#     file_utf8 = os.path.join(path_annotation_utf8, file)
#     f_in = open(file_ansi, 'rb')
#     str_check = f_in.read()
#     f_in.close()
#     res_check = chardet.detect(str_check)
#     if res_check == 'utf-8':
#         os.system(f"cp {file_ansi} {file_utf8}")
#     else:
#         f = codecs.open(file_ansi, 'r', 'ansi')
#         ff = f.read()
#         file_object = codecs.open(file_utf8, 'w', 'utf-8')
#         file_object.write(ff)
#         f.close()
#         file_object.close()


count_files = os.listdir(path_count)
count_tissue = [file.split('.')[0] for file in count_files]
df_count = pd.DataFrame({'count_file': count_files, 'tissue': count_tissue})
annotation_files = os.listdir(path_annotation)
annotation_tissue = \
    [file.split('_')[0].replace('-', '') for file in annotation_files]
df_annotation = pd.DataFrame(
    {'annotation_file': annotation_files, 'tissue': annotation_tissue})
df_tissue = pd.merge(df_count, df_annotation, on='tissue', how='inner')
file_count_annotation = os.path.join(path_HCA, 'count_annotation.txt')
df_tissue.to_csv(file_count_annotation, sep='\t', index=None, na_rep='NA')


def sum_count_and_sample(dict_in):
    file_count = os.path.join(path_count, dict_in['count_file'])
    file_annotation = \
        os.path.join(path_annotation, dict_in['annotation_file'])
    file_sum = os.path.join(path_tissue, f"{dict_in['tissue']}.sum.txt")
    mtx_count = pd.read_csv(file_count, sep=' ', index_col=0)
    try:
        mtx_annotation = pd.read_csv(file_annotation, sep=',')
    except UnicodeDecodeError:
        mtx_annotation = pd.read_csv(file_annotation, sep=',', encoding='GBK')
    tissue = mtx_annotation.loc[1, 'Sample']
    mtx_annotation = mtx_annotation.loc[:, ['Cell_id', 'Celltype']]
    celltypes = mtx_annotation['Celltype'].unique().tolist()
    # list_series = []
    # for cell in celltypes:
    #     cell_ids = mtx_annotation.loc[
    #         mtx_annotation['Celltype'] == cell, 'Cell_id'].tolist()
    #     try:
    #         sub_mtx = np.sum(mtx_count.loc[:, cell_ids], axis=1)
    #     except KeyError:
    #         pre_cell_id = mtx_count.columns[0].split('.')[0]
    #         cell_ids = \
    #             [pre_cell_id + '.' + cell.split('.')[1] for cell in cell_ids]
    #         sub_mtx = np.sum(mtx_count.loc[:, cell_ids], axis=1)
    #     sub_mtx.name = cell + '_' + tissue
    #     list_series.append(sub_mtx)
    # mtx_out = pd.concat(list_series, sort=False, axis=1)
    # mtx_out.to_csv(file_sum, sep='\t')

    # sample
    all_cell_ids = list(mtx_count.columns)
    num_sample = int(len(all_cell_ids) / 100)
    print(dict_in['tissue'], '          ', num_sample)
    sample_cell_ids = random.sample(all_cell_ids, num_sample)
    sample_mtx_count = mtx_count.loc[:, sample_cell_ids]

    return sample_mtx_count


pool = Pool(30)
list_sample = pool.map(sum_count_and_sample, df_tissue.to_dict('records'))
pool.close()

df_sample = pd.concat(list_sample, axis=1, join='inner')
file_sample = os.path.join(path_combine, 'sample.txt')
df_sample.to_csv(file_sample, sep='\t')

sum_files = os.listdir(path_tissue)
list_df_sum = []
for file in sum_files:
    df_sum = pd.read_csv(
        os.path.join(path_tissue, file), sep='\t', index_col=0)
    list_df_sum.append(df_sum)

df_combine = pd.concat(list_df_sum, axis=1, sort=False)
list_cells = \
    [col_name.split('_')[0].strip() for col_name in df_combine.columns]
series_cells = pd.Series(list_cells, name='cell_type')
set_cells = set(list_cells)
list_combine = []
for cell_type in set_cells:
    col_names = list(series_cells[series_cells == cell_type].index)
    sub_combine = np.sum(df_combine.iloc[:, col_names], axis=1)
    sub_combine.name = cell_type
    list_combine.append(sub_combine)
df_combine_cell = pd.concat(list_combine, axis=1, sort=True)

file_combine = os.path.join(path_combine, 'HCA_combined.txt')
file_combine_cell = os.path.join(path_combine, 'HCA_combined_by_cell.txt')
df_combine.to_csv(file_combine, sep='\t')
df_combine_cell.to_csv(file_combine_cell, sep='\t')

time_end = time()
print(time_end - time_start)
