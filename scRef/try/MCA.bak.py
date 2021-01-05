#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: MCA.bak.py
# @time: 8/5/20 10:35 PM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
import random
from glob import glob


path_MCA = '/home/disk/scRef/MouseAtlas_SingleCell_Han2018'
path_count = os.path.join(path_MCA, 'rmbatch_dge')
path_combine = os.path.join(path_MCA, 'combinedMCA')


def sample_count(file):
    file_count = os.path.join(path_count, file)
    mtx_count = pd.read_csv(file_count, sep=' ', index_col=0)
    all_cell_ids = list(mtx_count.columns)
    num_sample = int(len(all_cell_ids) / 50)
    print(file, '          ', num_sample)
    sample_cell_ids = random.sample(all_cell_ids, num_sample)
    sample_mtx_count = mtx_count.loc[:, sample_cell_ids]

    return sample_mtx_count


file_counts = os.listdir(path_count)
pool = Pool(30)
list_sample = pool.map(sample_count, file_counts)
pool.close()

df_sample = pd.concat(list_sample, axis=1, join='inner')
file_sample = os.path.join(path_combine, 'sample.txt')
df_sample.to_csv(file_sample, sep='\t')


# uniform cell names
file_mat = os.path.join(path_combine, 'MCA_combined_mouse.txt')
df_mat = pd.read_csv(file_mat, sep='\t', index_col=0)
file_name = os.path.join(path_combine, 'MCA_combined_mouse_uniformCellName.txt')
df_name = pd.read_csv(file_name, sep='\t', header=None)
uniform_names = df_name[1].unique().tolist()
uniform_names.sort()

list_series = []
for cell_name in uniform_names:
    ori_names = df_name.loc[df_name[1] == cell_name, 0].tolist()
    if len(ori_names) == 1:
        sub_series = df_mat.loc[:, ori_names]
        sub_series.name = cell_name
        list_series.append(sub_series)
    else:
        sub_series = np.sum(df_mat.loc[:, ori_names], axis=1)
        sub_series.name = cell_name
        list_series.append(sub_series)

df_mat_uniform = pd.concat(list_series, axis=1)
file_uniform = os.path.join(path_combine, 'MCA_combined_mouse_uniform.txt')
df_mat_uniform.to_csv(file_uniform, sep='\t')

# combine cell ref from all tissue
path_MCA_ref = '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/MCA_ref'
path_MCA = '/home/disk/scRef/MouseAtlas_SingleCell_Han2018'
path_combine = os.path.join(path_MCA, 'combinedMCA')
file_concat_inner = os.path.join(path_combine, 'MCA_concat_inner.txt')
file_concat_outer = os.path.join(path_combine, 'MCA_concat_outer.txt')

files = glob(f"{path_MCA_ref}/*_mouse.txt")

list_df = []
for file in files:
    df_sub = pd.read_csv(file, sep='\t', index_col=0)
    list_df.append(df_sub)

df_concat_inner = pd.concat(list_df, axis=1, join='inner')
df_concat_inner.to_csv(file_concat_inner, sep='\t')
df_concat_outer = pd.concat(list_df, axis=1, join='outer')
df_concat_outer = df_concat_outer.fillna(0)
df_concat_outer.to_csv(file_concat_outer, sep='\t')

if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
