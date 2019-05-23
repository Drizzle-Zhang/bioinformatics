#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: PijuanSala.py
# @time: 4/22/19 8:55 PM

from time import time
import pandas as pd
import os
import re
from collections import defaultdict
import numpy as np
from multiprocessing import Pool
from functools import partial
import subprocess


def select_cells(sub_cells, list_subfile):
    df_subfile = pd.read_csv(list_subfile['filename'], sep='\t', index_col=0)
    sel_cells = set(sub_cells).intersection(set(df_subfile.columns))
    df_sel = df_subfile.loc[:, sel_cells]

    return df_sel


def exp_raw(file_counts, ref_file, file_cluster, path_output):
    df_ref = pd.read_csv(ref_file, sep='\t', dtype='object')
    df_ref = df_ref.dropna()
    df_ref.index = range(df_ref.shape[0])
    for i in range(df_ref.shape[0]):
        if type(df_ref.loc[i, 'CL_ID']) == str:
            df_ref.loc[i, 'CL_ID'] = df_ref.loc[i, 'CL_ID'].strip()
        df_ref.loc[i, 'Cell_type'] = df_ref.loc[i, 'Cell_type'].strip()

    df_counts = pd.read_csv(file_counts, sep='\t')
    idx_counts = df_counts.index
    idx_counts.name = 'Gene'
    df_counts.index = idx_counts
    df_counts.to_csv(os.path.join(path_output, 'cell_exp.txt'), sep='\t')

    df_cell_cluster = pd.read_csv(file_cluster, sep=',', dtype='object')
    uniq_barcodes = ['-'.join([df_cell_cluster.loc[idx, 'barcode'],
                               df_cell_cluster.loc[idx, 'sample']])
                     for idx in df_cell_cluster.index]
    df_cell_cluster['uniq_barcode'] = uniq_barcodes
    df_meta = df_cell_cluster.loc[:, ['uniq_barcode', 'celltype']]
    df_meta.columns = ['ID', 'Cluster']
    df_meta = df_meta.dropna()
    df_meta.to_csv(os.path.join(path_output, 'cell_meta.txt'),
                   sep='\t', na_rep='NA', index=False)

    num_cols = df_meta.shape[0]
    size_subfile = 1000
    path_subfiles = os.path.join(path_output, 'subfiles')
    if not os.path.exists(path_subfiles):
        os.mkdir(path_subfiles)
    num_file = num_cols // size_subfile + 1
    list_subfiles = []
    for i in range(num_file):
        if i == num_file - 1:
            list_subfiles.append(
                {'filename': os.path.join(path_subfiles, 'subfile' + str(i)),
                 'start': str(2 + i * size_subfile), 'end': str(num_cols)})
        else:
            list_subfiles.append(
                {'filename': os.path.join(path_subfiles, 'subfile' + str(i)),
                 'start': str(2 + i * size_subfile),
                 'end': str(1 + (i + 1) * size_subfile)})
    subprocesses = []
    process = 40
    for i in range(len(list_subfiles)):
        if i % process == 0:
            for sub in subprocesses:
                sub.wait()
            subprocesses = []
        subprocesses.append(
            subprocess.Popen(
                f"cut -f 1,{list_subfiles[i]['start']}-"
                f"{list_subfiles[i]['end']} {file_counts} > "
                f"{list_subfiles[i]['filename']}",
                shell=True))

    list_clusters = []
    for cluster in df_ref['Cell_type']:
        sub_cells = np.array(df_cell_cluster.loc[
            df_cell_cluster['Celltype'] == cluster, 'CellName']).tolist()
        pool = Pool(processes=40)
        func_sel = partial(select_cells, sub_cells)
        list_df_sel = pool.map(func_sel, list_subfiles)
        pool.close()
        df_cluster = pd.concat(list_df_sel, axis=1, sort=False)
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
    path_pattern = re.compile(r'/Mouse.+/?')
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
    exp_raw(
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/L2/L2_wy/'
        'MouseEmbryo_SingleCell_PijuanSala2019/raw_count.mtx',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/L2/L2_wy/'
        'MouseEmbryo_SingleCell_PijuanSala2019/cell_type.txt',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/L2/L2_wy/'
        'MouseEmbryo_SingleCell_PijuanSala2019/atlas/meta.csv',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/L2/L2_wy/'
        'MouseEmbryo_SingleCell_PijuanSala2019/')
    sum_by_cl(
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/L2/L2_wy/'
        'MouseEmbryo_SingleCell_PijuanSala2019/')
    time_end = time()
    print(time_end - time_start)
