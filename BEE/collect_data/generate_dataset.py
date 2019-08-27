#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_dataset.py
# @time: 8/20/19 9:03 PM

from time import time
import pandas as pd
from math import isnan
import numpy as np
import os
from multiprocessing import Pool
from functools import partial


def stat_na(file_mtx):
    file_mtx = '/lustre/tianlab/zhangyu/BEE/collect_data/mtx_lab_cell.txt'
    df_mtx = pd.read_csv(file_mtx, sep='\t', index_col=0)
    df_na = df_mtx.applymap(lambda x: isnan(x))
    cell_no_na = df_mtx.shape[1] - np.sum(df_na, axis=1)
    cell_no_na = cell_no_na.sort_values(ascending=False)
    col_no_na = df_mtx.shape[0] - np.sum(df_na, axis=0)
    col_no_na = col_no_na.sort_values(ascending=False)

    select_cell = cell_no_na.loc[cell_no_na > 1].index
    df_select = df_mtx.loc[select_cell, :]

    df_select_na = df_select.applymap(lambda x: isnan(x))
    col_no_na = df_select.shape[0] - np.sum(df_select_na, axis=0)
    col_no_na = col_no_na.sort_values(ascending=False)

    return


def extract_data(path_data, num_cell, subinput):
    sub_path = os.path.join(
        path_data, f"{subinput['folder']}/{subinput['cell']}")
    file_mat = os.path.join(sub_path, f"{subinput['cell']}_mtx.txt")
    df_mat = pd.read_csv(file_mat, sep='\t')
    df_mat = df_mat.rename(columns={'Unnamed: 0': 'Gene'})
    df_mat = df_mat.drop_duplicates(subset='Gene')
    df_mat.index = df_mat['Gene']
    df_mat = df_mat.drop(labels=['Gene'], axis=1)

    # try:
    #     df_mat_sample = df_mat.sample(n=num_cell, random_state=1234, axis=1)
    # except ValueError:
    #     print(subinput)
    if df_mat.shape[1] < 50:
        df_mat_sample = df_mat
        print(subinput)
    else:
        df_mat_sample = df_mat.sample(n=num_cell, random_state=1234, axis=1)

    return df_mat_sample


def data_set1(path_data, path_out):
    # the same sequencing platform; same cell numbers;
    # select several batches; different cell types for each batches
    # path_data = '/lustre/tianlab/zhangyu/BEE/collect_data'
    # path_out = ''
    # path_out = ''
    file_mtx = os.path.join(path_data, 'mtx_lab_cell.txt')
    file_info = os.path.join(path_data, 'lab_info.txt')

    df_mtx = pd.read_csv(file_mtx, sep='\t', index_col=0)
    df_info = pd.read_csv(file_info, sep='\t')
    # sequencing platform
    df_mtx = df_mtx.loc[
             :, df_info.loc[df_info['platform'] == '10X Genomics', 'folder']]
    df_na = df_mtx.applymap(lambda x: isnan(x))
    cell_no_na = df_mtx.shape[1] - np.sum(df_na, axis=1)
    cell_no_na = cell_no_na.sort_values(ascending=False)

    # condition1
    select_cell = cell_no_na.loc[cell_no_na > 1].index
    df_select = df_mtx.loc[select_cell, :]
    df_select_na = df_select.applymap(lambda x: isnan(x))
    col_no_na = df_select.shape[0] - np.sum(df_select_na, axis=0)
    col_no_na = col_no_na.sort_values(ascending=False)

    # condition2
    df_select = df_select.loc[:, col_no_na.loc[col_no_na > 1].index]
    num_cell = int(np.min(np.min(df_select)))
    num_cell = 50
    folders = df_select.columns
    cells = df_select.index
    list_input = []
    for folder in folders:
        for cell in cells:
            if isnan(df_select.loc[cell, folder]):
                continue
            else:
                list_input.append(
                    {'folder': folder,
                     'cell': cell.replace(' ', '_').replace('/', '.')}
                )

    pool = Pool(processes=20)
    func_extract = partial(extract_data, path_data, num_cell)
    list_mat = pool.map(func_extract, list_input)
    pool.close()

    list_meta = []
    for idx, df_mat in enumerate(list_mat):
        df_meta = pd.DataFrame(columns=['column', 'folder', 'cell'])
        df_meta['column'] = df_mat.columns
        df_meta['folder'] = list_input[idx]['folder']
        df_meta['cell'] = \
            list_input[idx]['cell'].replace('_', ' ').replace('.', '/')
        list_meta.append(df_meta)

    out_mat = pd.concat(list_mat, axis=1, join='inner')
    out_meta = pd.concat(list_meta, axis=0)
    out_mat.to_csv(os.path.join(path_out, 'matrix.txt'), sep='\t')
    out_meta.to_csv(os.path.join(path_out, 'meta.txt'), sep='\t', index=None)

    return


def data_set2(path_data, path_out):
    # same cell numbers;
    # select several batches using different sequencing platform;
    # same cell types for each batches
    # path_data = '/lustre/tianlab/zhangyu/BEE/collect_data'
    # path_out = ''
    # path_out = ''
    file_mtx = os.path.join(path_data, 'mtx_lab_cell.txt')
    file_info = os.path.join(path_data, 'lab_info.txt')

    df_mtx = pd.read_csv(file_mtx, sep='\t', index_col=0)
    df_info = pd.read_csv(file_info, sep='\t')
    # select dataset
    select_folder = \
        ['MouseAtlas_SingleCell_Han2018', 'MouseEmbryo_SingleCell_Cao2019',
         'MouseBrainSpinalcordAtlas_SingleCell_Rosenberg2018',
         'MouseBrainDentateGyrus_SingleCell_Hochgerner2018',
         'MouseCerebralcortex_SingleCell_Loo2018',
         'MouseNeocortex_SingleCell_Tasic2018']
    df_mtx = df_mtx.loc[:, select_folder]
    df_na = df_mtx.applymap(lambda x: isnan(x))
    cell_no_na = df_mtx.shape[1] - np.sum(df_na, axis=1)
    cell_no_na = cell_no_na.sort_values(ascending=False)

    # condition1
    select_cell = cell_no_na.loc[cell_no_na > 2].index
    df_select = df_mtx.loc[select_cell, :]
    df_select_na = df_select.applymap(lambda x: isnan(x))
    col_no_na = df_select.shape[0] - np.sum(df_select_na, axis=0)
    col_no_na = col_no_na.sort_values(ascending=False)

    # condition2
    df_select = df_select.loc[:, col_no_na.loc[col_no_na > 1].index]
    # sel_50 = np.sum(df_select < 50, axis=1)
    # cell_drop = sel_50.loc[sel_50 > 0].index.tolist()
    # df_select = df_select.drop(labels=cell_drop, axis=0)
    num_cell = int(np.min(np.min(df_select)))
    num_cell = 50
    folders = df_select.columns
    cells = df_select.index
    list_input = []
    for folder in folders:
        for cell in cells:
            if isnan(df_select.loc[cell, folder]):
                continue
            else:
                list_input.append(
                    {'folder': folder,
                     'cell': cell.replace(' ', '_').replace('/', '.')}
                )

    pool = Pool(processes=20)
    func_extract = partial(extract_data, path_data, num_cell)
    list_mat = pool.map(func_extract, list_input)
    pool.close()

    list_meta = []
    for idx, df_mat in enumerate(list_mat):
        df_meta = pd.DataFrame(columns=['column', 'folder', 'cell'])
        df_meta['column'] = df_mat.columns
        df_meta['folder'] = list_input[idx]['folder']
        df_meta['cell'] = \
            list_input[idx]['cell'].replace('_', ' ').replace('.', '/')
        list_meta.append(df_meta)

    path_out = '/lustre/tianlab/zhangyu/BEE/real_data/multi_platform_balance'
    out_mat = pd.concat(list_mat, axis=1, join='inner')
    out_meta = pd.concat(list_meta, axis=0)
    out_mat.to_csv(os.path.join(path_out, 'matrix.txt'), sep='\t')
    out_meta.to_csv(os.path.join(path_out, 'meta.txt'), sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    path_input = '/lustre/tianlab/zhangyu/BEE/collect_data'
    path_output = '/lustre/tianlab/zhangyu/BEE/real_data/data1'

    time_end = time()
    print(time_end - time_start)
