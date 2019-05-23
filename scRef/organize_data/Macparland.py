#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: Macparland.py
# @time: 4/9/19 8:56 PM

from time import time
import pandas as pd
import os
import re
from collections import defaultdict
import numpy as np
from multiprocessing import Pool
from functools import partial


def exp_cell_sub_folder(path_count, selected_barcode, folder_nm):
    if folder_nm[-6:] == 'select':
        return

    barcode_pattern = re.compile(r'barcodes.tsv')
    gene_pattern = re.compile(r'genes.tsv')
    mat_pattern = re.compile(r'matrix.mtx')

    folder = os.path.join(path_count, folder_nm)
    if os.path.isdir(folder):
        for file in os.listdir(folder):
            if barcode_pattern.search(file):
                os.rename(os.path.join(folder, file),
                          os.path.join(folder, 'barcodes.tsv'))
            elif gene_pattern.search(file):
                os.rename(os.path.join(folder, file),
                          os.path.join(folder, 'genes.tsv'))
            elif mat_pattern.search(file):
                os.rename(os.path.join(folder, file),
                          os.path.join(folder, 'matrix.mtx'))

        folder_select = os.path.join(path_count, f"{folder_nm}_select")
        if not os.path.exists(folder_select):
            os.mkdir(folder_select)
        dict_barcode_colnum = {}
        with open(os.path.join(folder, 'barcodes.tsv'), 'r') as r_bar:
            with open(os.path.join(folder_select, 'barcodes.tsv'),
                      'w') as w_bar:
                for idx, line in enumerate(r_bar):
                    barcode = f"{folder_nm}_{line.strip()}"
                    barcode = barcode.replace('-', '_')
                    if barcode in selected_barcode:
                        w_bar.write(f"{barcode}\n")
                        dict_barcode_colnum[barcode] = idx + 1
        set_select_colnum = set(dict_barcode_colnum.values())
        dict_colnum_newcolnum = {}
        for idx, colnum in enumerate(dict_barcode_colnum.values()):
            dict_colnum_newcolnum[colnum] = idx + 1

        with open(os.path.join(folder, 'matrix.mtx'), 'r') as r_f:
            with open(os.path.join(folder_select, 'matrix_tmp.mtx'),
                      'w') as w_f:
                for i, line in enumerate(r_f):
                    if i > 2:
                        tmp_line = line.strip().split(' ')
                        col = int(tmp_line[1])
                        if col in set_select_colnum:
                            w_f.write(f"{tmp_line[0]} "
                                      f"{dict_colnum_newcolnum[col]} "
                                      f"{tmp_line[2]}\n")

        sub_mtx = pd.read_csv(
            os.path.join(folder_select, 'matrix_tmp.mtx'),
            sep=' ', header=None)
        max_col = np.max(sub_mtx.iloc[:, 1])
        sum_counts = np.sum(sub_mtx.iloc[:, 2])
        max_row = 33694
        with open(os.path.join(folder_select, 'matrix_header.mtx'),
                  'w') as w_f:
            w_f.write('%%MatrixMarket matrix coordinate integer general\n'
                      '%\n')
            w_f.write(f"{max_row} {max_col} {sum_counts}\n")

        os.system(f"cat "
                  f"{os.path.join(folder_select, 'matrix_header.mtx')} "
                  f"{os.path.join(folder_select, 'matrix_tmp.mtx')} > "
                  f"{os.path.join(folder_select, 'matrix.mtx')}")
        os.remove(os.path.join(folder_select, 'matrix_header.mtx'))
        os.remove(os.path.join(folder_select, 'matrix_tmp.mtx'))
        os.system(f"cp {os.path.join(folder, 'genes.tsv')} "
                  f"{folder_select}")

        os.system(f'Rscript tran_mat.R {folder_select}')
        df_folder = pd.read_csv(os.path.join(folder_select, 'res_mnt.txt'),
                                sep='\t')
        return df_folder
    else:
        return


def exp_cell(path_count, selected_file):
    selected_mat = pd.read_csv(selected_file, sep=',', index_col=0)
    selected_barcode = set(selected_mat.columns)

    pool = Pool(processes=12)
    func = partial(exp_cell_sub_folder,
                   path_count, selected_barcode)
    list_mats = pool.map(func, os.listdir(path_count))
    pool.close()

    df_count = pd.concat(list_mats, axis=1, sort=False)

    idx = df_count.index
    idx.name = "Gene"
    df_count.index = idx
    df_count.to_csv(os.path.join(path_count, 'cell_exp.txt'), sep='\t')

    return


def exp_raw(file_counts, ref_file, file_cluster, file_celltype, path_output):
    df_ref = pd.read_csv(ref_file, sep='\t', dtype='object')
    df_counts = pd.read_csv(file_counts, sep='\t', index_col=0)
    df_cell_cluster = pd.read_csv(file_cluster, sep='\t', dtype='object')
    df_cluster_celltype = pd.read_csv(file_celltype, sep='\t', dtype='object')

    cell_types = []
    for num in df_cell_cluster['Cluster#']:
        cell_types.append(
            df_cluster_celltype.loc[
                df_cluster_celltype['Cluster'] == num, 'Celltype'].values[0])
    df_cell_cluster['Celltype'] = cell_types
    df_meta = df_cell_cluster.loc[:, ['CellName', 'Celltype']]
    df_meta.columns = ['ID', 'Cluster']
    df_meta.to_csv(os.path.join(path_output, 'cell_meta.txt'),
                   sep='\t', na_rep='NA', index=False)

    list_clusters = []
    for cluster in df_ref['Cell_type']:
        sub_cells = np.array(df_cell_cluster.loc[
            df_cell_cluster['Celltype'] == cluster, 'CellName']).tolist()
        df_cluster = df_counts.loc[:, sub_cells]
        se_cluster = np.sum(df_cluster, axis=1)
        folder = np.array(
            df_ref.loc[df_ref['Cell_type'] == cluster, 'Dir']).tolist()
        cl_id = np.array(
            df_ref.loc[df_ref['Cell_type'] == cluster, 'CL_ID']).tolist()
        se_cluster.name = \
            f"{folder[0]}|cell_exp.txt|{cluster}|{','.join(cl_id)}"
        list_clusters.append(se_cluster)
    df_clusters = pd.concat(list_clusters, axis=1, sort=False)

    df_clusters.to_csv(os.path.join(path_output, 'exp_raw.txt'),
                       sep='\t', na_rep='NA')

    return


def sum_by_cl(path):
    df_raw = pd.read_csv(os.path.join(path, 'exp_raw.txt'), sep='\t',
                         index_col=0)
    path_pattern = re.compile(r'/Human.+/?')
    pattern_2cl = re.compile(r',')
    path_mouse = path_pattern.search(path).group()[1:-1]
    cols = df_raw.columns
    dict_cl = defaultdict(list)
    for col in cols:
        list_col = col.strip().split('|')
        if list_col[-1][0:2] == "CL":
            if pattern_2cl.search(col):
                cls = col.strip().split(',')
                dict_cl[cls[0]].append(col)
                dict_cl[cls[1]].append(col)
            else:
                dict_cl[list_col[-1]].append(col)

    list_cl = []
    for cl_id in dict_cl.keys():
        colnms = dict_cl[cl_id]
        df_cl = df_raw.loc[:, colnms]
        cl_sum = np.sum(df_cl, axis=1)
        cl_sum.name = f"{path_mouse}|{cl_id}"
        cl_sum.index = df_raw.index
        list_cl.append(cl_sum)

    df_cls = pd.concat(list_cl, axis=1)
    df_cls.to_csv(os.path.join(path, 'exp_sum.txt'), sep='\t', na_rep='NA')

    return


if __name__ == '__main__':
    time_start = time()
    """
    exp_cell(
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/Counts',
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/GSE115469_Data.csv')
    """
    exp_raw(
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/Counts/cell_exp.txt',
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/cell_type.txt',
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/Cell_clusterID_cycle.txt',
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/cluster_celltype.txt',
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/')
    sum_by_cl(
        '/home/zy/single_cell/data_human/L2/'
        'HumanLiver_SingleCell_MacParland2018/')
    time_end = time()
    print(time_end - time_start)
