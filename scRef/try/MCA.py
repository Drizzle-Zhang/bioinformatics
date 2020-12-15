#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: MCA.py
# @time: 8/26/20 4:21 PM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from glob import glob


# combine expression profile by folder
def combine_by_folder(folder):
    tissues = (df_info.loc[df_info['folder'] == folder, 'tissue']).tolist()
    list_count = []
    list_anno = []
    for tissue in tissues:
        count_file = os.path.join(path_MCA_by_tissue, f"{folder}/{tissue}.txt")
        annotation_file = \
            os.path.join(path_MCA_by_tissue,
                         f"{folder}/CellAssignments_{tissue}.txt")
        if not os.path.exists(count_file):
            continue
        if not os.path.exists(annotation_file):
            continue
        df_count = pd.read_csv(count_file, sep='\t', index_col=0)
        df_annotation = pd.read_csv(annotation_file, sep='\t', index_col=0)
        overlap_cells = \
            list(set(df_count.columns).intersection(set(df_annotation.index)))
        df_count_overlap = df_count.loc[:, overlap_cells]
        df_annotation_overlap = df_annotation.loc[overlap_cells, :]
        list_count.append(df_count_overlap)
        list_anno.append(df_annotation_overlap)

    if len(list_count) < 1:
        print(folder)
        return
    sub_count_combine = pd.concat(list_count, axis=1, join='outer', sort=False)
    sub_count_combine = sub_count_combine.fillna(0)
    sub_anno_combine = pd.concat(list_anno, axis=0, join='inner')
    sub_anno_combine['CellType'] = ['NA'] * sub_anno_combine.shape[0]
    sub_anno_combine['TissueCellType'] = ['NA'] * sub_anno_combine.shape[0]
    for cellid in sub_anno_combine.index:
        tissue_name = sub_anno_combine.loc[cellid, 'Tissue']
        annotation = sub_anno_combine.loc[cellid, 'Annotation']
        cell_type = annotation.split('(')[0].split('_')[0]
        sub_anno_combine.loc[cellid, 'CellType'] = cell_type
        sub_anno_combine.loc[cellid, 'TissueCellType'] = \
            f"{tissue_name}|{cell_type}"

    cell_types = sub_anno_combine['TissueCellType'].unique().tolist()
    list_cell_num = []
    list_series = []
    for cell_type in cell_types:
        cell_ids = sub_anno_combine.loc[
            sub_anno_combine['TissueCellType'] == cell_type,
            'Cell.name'].tolist()
        list_cell_num.append({'TissueCellType': cell_type,
                              'CellNum': len(cell_ids)})
        if type_combine == 'mean':
            sub_mtx = np.mean(sub_count_combine.loc[:, cell_ids], axis=1)
        elif type_combine == 'sum':
            sub_mtx = np.sum(sub_count_combine.loc[:, cell_ids], axis=1)
        else:
            print('Please set type of combine!')
            # return
        sub_mtx.name = cell_type
        list_series.append(sub_mtx)
    mtx_out = pd.concat(list_series, sort=False, axis=1)
    mtx_out = mtx_out.astype(np.float16)

    annotation_out = sub_anno_combine.loc[:,
                     ['Tissue', 'Annotation', 'CellType', 'TissueCellType']]
    annotation_out = annotation_out.drop_duplicates()
    df_num = pd.DataFrame(list_cell_num)
    annotation_out = pd.merge(
        annotation_out, df_num, on='TissueCellType', how='left')

    # write files
    file_combine = os.path.join(
        path_MCA_by_tissue, f"{folder}/Combine_{type_combine}.txt")
    mtx_out.to_csv(file_combine, sep='\t')
    file_count = os.path.join(
        path_MCA_by_tissue, f"{folder}/Count_all_batch.txt")
    file_annotation = os.path.join(
        path_MCA_by_tissue, f"{folder}/CellAssignments_all_batch.txt")
    sub_count_combine.to_csv(file_count, sep='\t')
    sub_anno_combine.to_csv(file_annotation, sep='\t')

    return mtx_out, annotation_out


if __name__ == '__main__':
    time_start = time()
    # generate cell expression values by tissue
    path_MCA = '/home/disk/scRef/MouseAtlas_SingleCell_Han2018'
    path_MCA_by_tissue = os.path.join(path_MCA, 'MCA_by_tissue')
    path_combine = os.path.join(path_MCA, 'combinedMCA')
    type_combine = 'mean'

    file_info = os.path.join(path_MCA_by_tissue, 'folder.info')
    df_info = pd.read_csv(file_info, sep='\t')
    folders = df_info['folder'].unique().tolist()

    pool = Pool(10)
    list_tuple = pool.map(combine_by_folder, folders)
    pool.close()

    list_1 = [sub_tuple[0] for sub_tuple in list_tuple if sub_tuple]
    list_2 = [sub_tuple[1] for sub_tuple in list_tuple if sub_tuple]
    file_combine_all = \
        os.path.join(path_combine, f'Count_all_tissue_{type_combine}.txt')
    count_combine_all = pd.concat(list_1, axis=1, join='outer', sort=False)
    count_combine_all = count_combine_all.fillna(0)
    count_combine_all.to_csv(file_combine_all, sep='\t')

    file_comine_meta = \
        os.path.join(path_combine, f'metadata_all_tissue.txt')
    df_meta = pd.concat(list_2, axis=0, join='outer', sort=False)
    df_meta.to_csv(file_comine_meta, sep='\t', index=None)

    time_end = time()
    print(time_end - time_start)
