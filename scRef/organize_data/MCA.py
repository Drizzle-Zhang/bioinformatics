#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: MCA.py
# @time: 4/16/19 10:33 PM

from time import time
import pandas as pd
import os
import re
from collections import defaultdict
import numpy as np
from multiprocessing import Pool
from functools import partial


def merge_count(path_count, output_file):
    for idx, file in enumerate(os.listdir(path_count)):
        sub_count = pd.read_csv(os.path.join(path_count, file), sep=' ')
        if idx == 0:
            df_merge = sub_count
        else:
            df_merge = pd.concat([df_merge, sub_count], axis=1, sort=False)
    df_merge.to_csv(output_file, sep='\t', na_rep='NA')

    return


def exp_raw_one_countfile(clusters, df_cell_cluster, df_ref,
                          path_counts, file_counts):
    df_counts = pd.read_csv(os.path.join(path_counts, file_counts),
                            sep=' ', index_col=0)
    count_cols = df_counts.columns
    list_clusters = []
    for cluster in clusters:
        sub_cells = np.array(df_cell_cluster.loc[
                                 df_cell_cluster[
                                     'Cluster'] == cluster, 'ID']).tolist()
        inter_cols = set(sub_cells).intersection(set(count_cols))
        if len(inter_cols) != 0:
            df_cluster = df_counts.loc[:, sub_cells]
            se_cluster = np.sum(df_cluster, axis=1)
            path_folder = np.array(
                df_ref.loc[df_ref['Cell_type'] == cluster, 'Dir']).tolist()
            cl_id = np.array(
                df_ref.loc[df_ref['Cell_type'] == cluster, 'CL_ID']).tolist()
            se_cluster.name = \
                f"{path_folder[0]}|{file_counts}|{cluster}|{','.join(cl_id)}"
            list_clusters.append(se_cluster)
    if list_clusters:
        df_clusters = pd.concat(list_clusters, axis=1, sort=False)
        return df_clusters
    else:
        return


def exp_raw(path_counts, ref_file, file_cluster, path_output):
    df_ref = pd.read_csv(ref_file, sep='\t', dtype='object')
    df_ref = df_ref.dropna()
    df_ref.index = range(df_ref.shape[0])
    for i in range(df_ref.shape[0]):
        if type(df_ref.loc[i, 'CL_ID']) == str:
            df_ref.loc[i, 'CL_ID'] = df_ref.loc[i, 'CL_ID'].strip()
        df_ref.loc[i, 'Cell_type'] = df_ref.loc[i, 'Cell_type'].strip()
    df_cell_cluster = pd.read_csv(file_cluster, sep=',',
                                  dtype='object', index_col=0)
    df_cell_cluster = df_cell_cluster.loc[:, ['Cell.name', 'Annotation']]
    df_cell_cluster.columns = ['ID', 'Cluster']
    df_cell_cluster.to_csv(os.path.join(path_output, 'cell_meta.txt'),
                           sep='\t', index=None)
    count_files = os.listdir(path_counts)
    clusters = set(np.array(df_ref['Cell_type']).tolist())

    pool = Pool(processes=20)
    func_one_countfile = partial(
        exp_raw_one_countfile, clusters, df_cell_cluster, df_ref, path_counts)
    list_countfiles = pool.map(func_one_countfile, count_files)
    pool.close()
    df_merge = pd.concat(list_countfiles, axis=1, sort=False)
    df_merge_idx = df_merge.index
    df_merge_idx.name = 'Gene'
    df_merge.index = df_merge_idx

    df_merge.to_csv(os.path.join(path_output, 'exp_raw.txt'),
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
    df_cls_idx = df_cls.index
    df_cls_idx.name = 'Gene'
    df_cls.index = df_cls_idx
    df_cls.to_csv(os.path.join(path, 'exp_sum.txt'), sep='\t', na_rep='NA')

    return


if __name__ == '__main__':
    time_start = time()
    """
    merge_count(
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/rmbatch_dge',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/cell_exp.txt')
    """
    exp_raw(
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/rmbatch_dge',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/cell_type.txt',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/MCA_CellAssignments.csv',
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018')

    sum_by_cl(
        '/home/disk/yangjing/tianlab/tools/scRef-master/Reference/'
        'MouseAtlas_SingleCell_Han2018/')

    time_end = time()
    print(time_end - time_start)
