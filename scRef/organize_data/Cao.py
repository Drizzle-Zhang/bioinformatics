#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: merge_files_Cao.py
# @time: 3/27/19 2:32 PM

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
import subprocess
from collections import defaultdict
import re
# import rpy2.robjects as ro


def make_seurat_file(cell_annotation, file_gene, file_matrix, path_seurat,
                     cell_type):
    sub_cell_annotation = cell_annotation.loc[
        cell_annotation['Main_cell_type'] == cell_type, 'sample']

    path_cell_type = os.path.join(
        path_seurat, cell_type.replace(' ', '_').replace('&', 'and'))
    print(cell_type)
    os.makedirs(path_cell_type)
    sub_cell_annotation.to_csv(
        os.path.join(path_cell_type, 'barcodes.tsv'),
        sep='\t', index=None, header=False)
    os.system(f"cp {os.path.join(path_seurat, file_gene)} "
              f"{path_cell_type}")
    sub_cols = np.array(sub_cell_annotation.index + 1).tolist()

    dict_col_idx = {}
    for idx, col in enumerate(sub_cols):
        dict_col_idx[col] = idx + 1
    set_sub_cols = set(sub_cols)

    with open(file_matrix, 'r') as r_f:
        with open(os.path.join(path_cell_type, 'matrix_tmp.mtx'),
                  'w') as w_f:
            for i, line in enumerate(r_f):
                if i > 1:
                    tmp_line = line.strip().split(' ')
                    col = int(tmp_line[1])
                    if col in set_sub_cols:
                        w_f.write(f"{tmp_line[0]}\t{dict_col_idx[col]}"
                                  f"\t{tmp_line[2]}\n")

    sub_mtx = pd.read_csv(
        os.path.join(path_cell_type, 'matrix_tmp.mtx'),
        sep='\t', header=None)
    max_col = np.max(sub_mtx.iloc[:, 1])
    sum_counts = np.sum(sub_mtx.iloc[:, 2])
    max_row = 26183
    with open(os.path.join(path_cell_type, 'matrix_header.mtx'),
              'w') as w_f:
        w_f.write('%%MatrixMarket matrix coordinate integer general\n')
        w_f.write(f"{max_row}\t{max_col}\t{sum_counts}\n")

    os.system(f"cat "
              f"{os.path.join(path_cell_type, 'matrix_header.mtx')} "
              f"{os.path.join(path_cell_type, 'matrix_tmp.mtx')} > "
              f"{os.path.join(path_cell_type, 'matrix.mtx')}")
    os.remove(os.path.join(path_cell_type, 'matrix_header.mtx'))
    os.remove(os.path.join(path_cell_type, 'matrix_tmp.mtx'))

    return


def generate_seurat_files(file_gene, file_cell, file_cleaned,
                          file_matrix, path_seurat):
    cell_annotation = pd.read_csv(file_cell, sep=',')
    cell_annotation['index'] = cell_annotation.index
    cell_annotation.index = cell_annotation['sample']
    with open(file_cleaned, 'r') as r_clean:
        cleaded_cells = [cell.strip() for cell in r_clean]
    cell_annotation_cleaded = cell_annotation.loc[cleaded_cells, :]
    cell_annotation_cleaded.index = cell_annotation_cleaded['index']
    cell_types = np.array(cell_annotation_cleaded['Main_cell_type']).tolist()
    cell_types = set([cell for cell in cell_types if cell is not np.nan])
    print(cell_types)

    pool = Pool(processes=20)
    func = partial(make_seurat_file, cell_annotation_cleaded, file_gene,
                   file_matrix, path_seurat)
    pool.map(func, cell_types)
    pool.close()

    return


def sum_count(list_subfile):
    df_sub = pd.read_csv(list_subfile['filename'], sep='\t')
    sub_merge = np.sum(df_sub, axis=1)

    return sub_merge


def split_sum(path, file_count):
    df_meta = pd.read_csv(os.path.join(path, 'barcodes.tsv'), sep='\t')
    num_cols = df_meta.shape[0]
    size_subfile = 1000
    path_subfiles = os.path.join(path, 'subfiles')
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
                f"{list_subfiles[i]['end']} "
                f"{os.path.join(path, file_count)} > "
                f"{list_subfiles[i]['filename']}", shell=True))

    for sub in subprocesses:
        sub.wait()

    pool = Pool(processes=40)
    list_df_sel = pool.map(sum_count, list_subfiles)
    pool.close()
    df_merge = pd.concat(list_df_sel, axis=1, sort=False)

    os.system(f"rm -rf {path_subfiles}")

    return df_merge


def trans_and_sum(path):
    if not os.path.exists(
            os.path.join(path, f"{path.strip().split('/')[-1]}_res_mnt.txt")):
        os.system(f'Rscript tran_mat.R {path}')

        os.system(f"head -1 {os.path.join(path, 'res_mnt_pre.txt')} "
                  f"> {os.path.join(path, 'res_header_pre.txt')}")
        os.system(f"sed -i '1d' {os.path.join(path, 'res_mnt_pre.txt')}")
        with open(os.path.join(path, 'res_header.txt'), 'w') as w_res:
            with open(os.path.join(path, 'res_header_pre.txt'), 'r') as r_mtx:
                w_res.write('Gene\t' + r_mtx.readline())
        os.system(f"cat {os.path.join(path, 'res_header.txt')} "
                  f"{os.path.join(path, 'res_mnt_pre.txt')} > "
                  f"{os.path.join(path, 'res_mnt.txt')}")
        os.system(f"rm {os.path.join(path, 'res_mnt_pre.txt')} "
                  f"{os.path.join(path, 'res_header_pre.txt')} "
                  f"{os.path.join(path, 'res_header.txt')}")

        cell = \
            path.strip().split('/')[-1].replace('_', ' ').replace('and', '&')

        cell_mat = split_sum(path, 'res_mnt.txt')
        cell_merge = np.sum(cell_mat, axis=1)
        cell_merge.name = cell
        os.rename(os.path.join(path, 'res_mnt.txt'),
                  os.path.join(path,
                               f"{path.strip().split('/')[-1]}_res_mnt.txt"))

    else:
        cell = \
            path.strip().split('/')[-1].replace('_', ' ').replace('and', '&')
        cell_mat = split_sum(
            path, f"{path.strip().split('/')[-1]}_res_mnt.txt")
        cell_merge = np.sum(cell_mat, axis=1)
        cell_merge.name = cell

    return cell_merge


def get_exp_raw_pre(path_seurat, path_output):
    folders = os.listdir(path_seurat)
    sub_paths = [os.path.join(path_seurat, folder) for folder in folders]

    for idx, path in enumerate(sub_paths):
        print(path)
        if idx == 0:
            df_merge = trans_and_sum(path)
        else:
            df_sub = trans_and_sum(path)
            df_merge = pd.concat([df_merge, df_sub], axis=1, sort=False)

    print(df_merge.shape)
    df_merge.to_csv(os.path.join(path_output, 'exp_raw_pre.txt'), sep='\t')

    return


def exp_raw(file_meta, file_cleaned, path_output):
    df_exp_pre = pd.read_csv(os.path.join(path_output, 'exp_raw_pre.txt'),
                             sep='\t', index_col=0)
    df_genes = pd.read_csv(os.path.join(path_output, 'genes.tsv'), sep='\t',
                           header=None)
    genes = df_genes.iloc[:, 1]
    genes.name = 'Gene'
    df_exp_pre.index = genes

    cell_annotation = pd.read_csv(file_meta, sep=',')
    cell_annotation.index = cell_annotation['sample']
    with open(file_cleaned, 'r') as r_clean:
        cleaded_cells = [cell.strip() for cell in r_clean]
    df_meta = cell_annotation.loc[cleaded_cells, ['sample', 'Main_cell_type']]
    df_meta.columns = ['ID', 'Cluster']
    df_meta = df_meta.dropna()
    df_meta.to_csv(os.path.join(path_output, 'cell_meta.txt'),
                   sep='\t', index=False)

    df_ref = pd.read_csv(os.path.join(path_output, 'cell_type.txt'), sep='\t')
    new_cols = []
    for col in df_exp_pre.columns:
        new_cols.append(f"MouseEmbryo_SingleCell_Cao2019|cell_exp|{col}|"
                        f"{df_ref.loc[df_ref['Cell_type'] == col, 'CL_ID'].values[0]}")
    df_exp_pre.columns = new_cols

    df_exp_pre.to_csv(os.path.join(path_output, 'exp_raw.txt'), sep='\t')

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
    """
    generate_seurat_files(
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019/genes.tsv',
        '/home/disk/scRef/L2/L2_wy/
        MouseEmbryo_SingleCell_Cao2019/cell_annotate.csv',
        '/home/disk/scRef/L2/L2_wy/
        MouseEmbryo_SingleCell_Cao2019/cleaned_cells.txt',
        '/home/disk/scRef/L2/L2_wy/
        MouseEmbryo_SingleCell_Cao2019/GSE119945_gene_count.txt',
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019/data_Seurat')
    
    get_exp_raw_pre(
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019/data_Seurat',
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019')
    """
    exp_raw(
        '/home/disk/scRef/L2/L2_wy/'
        'MouseEmbryo_SingleCell_Cao2019/cell_annotate.csv',
        '/home/disk/scRef/L2/L2_wy/'
        'MouseEmbryo_SingleCell_Cao2019/cleaned_cells.txt',
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019/')
    sum_by_cl(
        '/home/disk/scRef/L2/L2_wy/MouseEmbryo_SingleCell_Cao2019/')

    time_end = time()
    print(time_end - time_start)
