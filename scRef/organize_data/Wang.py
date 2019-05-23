#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: merge_files_Wang.py
# @time: 4/8/19 8:53 PM

from time import time
import pandas as pd
import os
import re
from collections import defaultdict
import numpy as np


def exp_raw(file_counts, ref_file, file_cluster, path_output):
    df_ref = pd.read_csv(ref_file, sep='\t', dtype='object')
    df_counts = pd.read_csv(file_counts, sep='\t', index_col=0)
    df_cell_cluster = pd.read_csv(file_cluster, sep='\t', dtype='object')

    col_pattern = re.compile(r'_\d+w')
    pre_pattern = re.compile(r'.+?_')
    suf_pattern = re.compile(r'w_.+')
    new_cols = []
    for col in df_counts.columns:
        if col_pattern.search(col):
            week = int(col_pattern.search(col).group()[1:-1]) - 2
            prefix = pre_pattern.match(col).group()
            suffix = suf_pattern.search(col).group()
            new_cols.append(f"{prefix}{week}{suffix}")
        else:
            new_cols.append(col)
    df_counts.columns = new_cols
    df_counts.to_csv(os.path.join(path_output, 'cell_exp.txt'),
                     sep='\t', na_rep='NA')
    df_meta = df_cell_cluster.loc[:, ['Cell', 'Cluster']]
    df_meta.columns = ['ID', 'Cluster']
    df_meta.to_csv(os.path.join(path_output, 'cell_meta.txt'),
                   sep='\t', na_rep='NA', index=False)

    clusters = set(np.array(df_ref['Cell_type']).tolist())
    list_clusters = []
    for cluster in clusters:
        sub_cells = np.array(df_cell_cluster.loc[
            df_cell_cluster['Cluster'] == cluster, 'Cell']).tolist()
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
    exp_raw(
        '/home/zy/single_cell/data_human/L2/HumanKidney_SingleCell_Wang2018/'
        'counts.txt',
        '/home/zy/single_cell/data_human/L2/HumanKidney_SingleCell_Wang2018/'
        'cell_type.txt',
        '/home/zy/single_cell/data_human/L2/HumanKidney_SingleCell_Wang2018/'
        'GSE2944332_meta.txt',
        '/home/zy/single_cell/data_human/L2/HumanKidney_SingleCell_Wang2018/')
    sum_by_cl(
        '/home/zy/single_cell/data_human/L2/HumanKidney_SingleCell_Wang2018/')

    time_end = time()
    print(time_end - time_start)
