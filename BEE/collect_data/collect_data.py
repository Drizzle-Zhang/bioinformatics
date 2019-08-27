#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: collect_data.py
# @time: 7/23/19 8:10 PM

from time import time
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial


def get_onecelldata_from_mat(select_cols, path_exp, sep_symbol,
                             path_subfiles, subfile):
    os.system(f"cut -d '{sep_symbol}' "
              f"-f 1,{subfile['start']}-{subfile['end']} "
              f"{path_exp} > "
              f"{os.path.join(path_subfiles, subfile['filename'])}")
    df_exp = pd.read_csv(
        os.path.join(path_subfiles, subfile['filename']),
        sep=sep_symbol, index_col=0)
    cols = set(select_cols).intersection(set(df_exp.columns))
    select_exp = df_exp.loc[:, cols]

    return select_exp


def split_merge_data(df_meta, path_exp, raw_colnms, path_folder_cl_name):
    if len(raw_colnms) == 1:
        select_cols = np.array(df_meta.loc[
                                   df_meta['Cluster'] == raw_colnms[0], 'ID'],
                               dtype='str').tolist()
    else:
        select_cols = []
        for col in raw_colnms:
            select_col = np.array(
                df_meta.loc[df_meta['Cluster'] == col, 'ID']).tolist()
            select_cols.extend(select_col)

    df_meta.index = df_meta['ID']
    select_meta = df_meta.loc[select_cols, :]

    list_path_exp = os.path.split(path_exp)
    if list_path_exp[1] != 'cell_exp.txt':
        pre_mtx = os.path.join(
            path_folder_cl_name, 'premtx_' + list_path_exp[1])
        header = os.path.join(
            path_folder_cl_name, 'header_' + list_path_exp[1])
        pre_header = os.path.join(
            path_folder_cl_name, 'preheader_' + list_path_exp[1])
        os.system(f"head -1 {path_exp} > {pre_header}")

        with open(pre_header, 'r') as r_mtx:
            line = r_mtx.readline()
            list_tab = line.strip().split('\t')
            list_comma = line.strip().split(',')
            list_blank = line.strip().split(' ')
            if len(list_tab) > 1 and len(list_comma) == 1 and len(list_blank) == 1:
                sep_symbol = '\t'
                num_cols = len(list_tab)
            elif len(list_tab) == 1 and len(list_comma) > 1 and len(list_blank) == 1:
                sep_symbol = ','
                num_cols = len(list_comma)
            elif len(list_tab) == 1 and len(list_comma) == 1 and len(list_blank) > 1:
                sep_symbol = ' '
                num_cols = len(list_blank)
            else:
                print(f"sep error {path_folder_cl_name}")
                return

            if line[0:4] != 'Gene':
                with open(header, 'w') as w_res:
                    w_res.write('Gene' + sep_symbol + line)
                os.system(f"cp {path_exp} {pre_mtx}")
                os.system(f"sed -i '1d' {pre_mtx}")
                os.system(f"cat {header} {pre_mtx} > {path_exp}")
                os.system(f"rm {pre_header} {header} {pre_mtx}")
            else:
                num_cols -= 1
                os.system(f"rm {pre_header}")
    else:
        num_cols = df_meta.shape[0]
        sep_symbol = '\t'

    size_subfile = 3000
    if num_cols < size_subfile:
        df_exp = pd.read_csv(path_exp, sep=sep_symbol, index_col=0)
        cols = set(select_cols).intersection(set(df_exp.columns))
        select_exp = df_exp.loc[:, cols]
    else:
        path_subfiles = \
                os.path.join(path_folder_cl_name, 'subfiles')
        if not os.path.exists(path_subfiles):
            os.mkdir(path_subfiles)
        num_file = num_cols // size_subfile + 1
        list_subfiles = []
        for i in range(num_file):
            if i == num_file - 1:
                list_subfiles.append({'filename': 'subfile' + str(i),
                                      'start': str(2 + i * size_subfile),
                                      'end': str(num_cols)})
            else:
                list_subfiles.append({'filename': 'subfile' + str(i),
                                      'start': str(2 + i * size_subfile),
                                      'end': str(1 + (i + 1) * size_subfile)})
        list_select_exp = []
        for subfile in list_subfiles:
            expre = get_onecelldata_from_mat(
                select_cols, path_exp, sep_symbol, path_subfiles, subfile
            )
            list_select_exp.append(expre)
        select_exp = pd.concat(list_select_exp, axis=1, sort=False)
        os.system(f"rm -rf {path_subfiles}")
        # select_exp = select_exp.dropna()

    return select_exp, select_meta


def get_one_mat(path_folder, out_folder, df_ref, df_meta, cl_name):
    if cl_name is np.nan:
        print(f'{cl_name}         NA')
        return
    else:
        print(cl_name)
        cl_file = cl_name.replace(' ', '_').replace('/', '.')
        path_folder_cl_name = os.path.join(out_folder, cl_file)
        if not os.path.exists(path_folder_cl_name):
            os.mkdir(path_folder_cl_name)
        pre_mtx_name = f"{cl_file}_pre_mtx.txt"
        pre_meta_name = f"{cl_file}_pre_meta.txt"
        pre_path_mtx = os.path.join(path_folder_cl_name, pre_mtx_name)
        pre_path_meta = os.path.join(path_folder_cl_name, pre_meta_name)
        mtx_name = f"{cl_file}_mtx.txt"
        meta_name = f"{cl_file}_meta.txt"
        path_mtx = os.path.join(path_folder_cl_name, mtx_name)
        path_meta = os.path.join(path_folder_cl_name, meta_name)
        if os.path.exists(path_mtx) and os.path.exists(path_meta):
            print(f"{cl_name}  existed")
            return

    df_ref_cl = df_ref.loc[df_ref['CL_name'] == cl_name, :]
    try:
        raw_colnms = np.array(df_ref_cl['Cell_type']).tolist()
    except KeyError:
        print(f'{cl_name}         Cell_type  error')
        return
    if len(raw_colnms) < 1:
        print(f'{cl_name}         none')
        return

    if os.path.exists(f"{path_folder}/raw_count/cell_exp"):
        folder_exp = f"{path_folder}/raw_count/cell_exp"
        sub_exps = []
        sub_metas = []
        for file_exp in os.listdir(folder_exp):
            path_exp = os.path.join(folder_exp, file_exp)
            sub_exp, sub_meta = split_merge_data(
                df_meta, path_exp, raw_colnms, path_folder_cl_name
            )
            sub_exps.append(sub_exp)
            sub_metas.append(sub_meta)
        select_exp = pd.concat(sub_exps, axis=1, sort=False, join='outer')
        select_meta = pd.concat(sub_metas, axis=0, sort=False, join='outer')
        # select_exp = select_exp.fillna(0)
    elif os.path.exists(f"{path_folder}/raw_count/cell_exp.txt"):
        path_exp = f"{path_folder}/raw_count/cell_exp.txt"
        select_exp, select_meta = split_merge_data(
            df_meta, path_exp, raw_colnms, path_folder_cl_name
        )
    else:
        print(f"{cl_name}         error: no count matrix file")
        return

    select_exp.to_csv(pre_path_mtx, sep='\t', na_rep='NA')
    os.rename(pre_path_mtx, path_mtx)
    select_meta.to_csv(pre_path_meta, sep='\t', index=None)
    os.rename(pre_path_meta, path_meta)

    return


def count_cell(out_folder, cl_name):
    if cl_name is np.nan:
        print(f'{cl_name}         NA')
        return
    else:
        print(cl_name)
        cl_file = cl_name.replace(' ', '_').replace('/', '.')
        path_folder_cl_name = os.path.join(out_folder, cl_file)
        mtx_name = f"{cl_file}_mtx.txt"
        path_mtx = os.path.join(path_folder_cl_name, mtx_name)

        df_meta = pd.read_csv(path_mtx, sep='\t', index_col=0)

        return {'cell_num': df_meta.shape[1], 'cell_name': cl_name}


def generate_mat_by_lab(path_data, path_output):

    list_folder = os.listdir(path_data)
    for folder in list_folder:
        path_folder = os.path.join(path_data, folder)
        try:
            df_ref = pd.read_csv(
                os.path.join(path_folder, 'cell_type.txt'),
                sep='\t', dtype='object', index_col=False)
        except FileNotFoundError:
            print(folder, 'cell_type.txt  error')
            continue

        try:
            df_ref = df_ref.loc[df_ref['Relationship'] == 'equal', :]
        except KeyError:
            print(folder, 'cell_type.txt  error  Relationship')

        path_meta = f"{path_folder}/raw_count/cell_meta.txt"
        try:
            df_meta = pd.read_csv(path_meta, sep='\t', dtype='object')
        except FileNotFoundError:
            print(folder, 'cell_meta  error')
            continue
        # if df_meta.shape[0] < 8000:
        #     continue

        list_lab = folder.strip().split('_')
        year = list_lab[2][-4:]
        print(folder)
        # if int(year) >= 2017:
        #     print(folder)
        # else:
        #     continue

        out_folder = os.path.join(path_output, folder)
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        cl_names = set(np.array(df_ref['CL_name']).tolist())
        if df_meta.shape[0] < 500000:
            pool = Pool(processes=40)
        else:
            pool = Pool(processes=5)
        func_get_one_mat = partial(
            get_one_mat, path_folder, out_folder, df_ref, df_meta
        )
        pool.map(func_get_one_mat, cl_names)
        pool.close()

    list_lab_info = []
    list_lab_cell = []
    for folder in list_folder:
        path_folder = os.path.join(path_data, folder)

        try:
            df_ref = pd.read_csv(
                os.path.join(path_folder, 'cell_type.txt'),
                sep='\t', dtype='object', index_col=False)
        except FileNotFoundError:
            print(folder, 'cell_type.txt  error')
            continue

        try:
            df_ref = df_ref.loc[df_ref['Relationship'] == 'equal', :]
        except KeyError:
            print(folder, 'cell_type.txt  error  Relationship')

        path_meta = f"{path_folder}/raw_count/cell_meta.txt"
        try:
            df_meta = pd.read_csv(path_meta, sep='\t', dtype='object')
        except FileNotFoundError:
            print(folder, 'cell_meta  error')
            continue
        # if df_meta.shape[0] < 8000:
        #     continue

        list_lab = folder.strip().split('_')
        year = list_lab[2][-4:]
        list_lab_info.append({'folder': folder,
                              'lab_lastname': list_lab[2][:-4],
                              'year': year,
                              'organ': list_lab[0]})
        #
        # if int(year) >= 2017:
        #     list_lab_info.append({'folder': folder,
        #                           'lab_lastname': list_lab[2][:-4],
        #                           'year': year,
        #                           'organ': list_lab[0]})
        # else:
        #     continue

        out_folder = os.path.join(path_output, folder)
        cl_names = set(np.array(df_ref['CL_name']).tolist())
        pool = Pool(processes=40)
        func_count_cell = partial(count_cell, out_folder)
        list_cells = pool.map(func_count_cell, cl_names)
        pool.close()

        df_cells = pd.DataFrame(list_cells)
        df_cells.index = df_cells['cell_name']
        df_cells = df_cells['cell_num']
        list_lab_cell.append(df_cells)

    df_lab_info = pd.DataFrame(list_lab_info)
    df_lab_cell = pd.concat(list_lab_cell, axis=1, sort=False, join='outer')
    df_lab_cell.columns = df_lab_info['folder']

    df_lab_info.to_csv(
        os.path.join(path_output, 'lab_info.txt'), sep='\t', index=None
    )
    df_lab_cell.to_csv(
        os.path.join(path_output, 'mtx_lab_cell.txt'), sep='\t'
    )

    return


if __name__ == '__main__':
    time_start = time()
    generate_mat_by_lab(
        '/lustre/tianlab/scRef/MouseReference_v1',
        '/lustre/tianlab/zhangyu/BEE/collect_data')

    time_end = time()
    print(time_end - time_start)
