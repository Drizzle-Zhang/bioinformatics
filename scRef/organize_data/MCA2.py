#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: MCA2.py
# @time: 9/23/19 11:01 AM

from time import time
import pandas as pd
import os
import re
from multiprocessing import Pool
from functools import partial


def trans_sep(info_folder, path_data, path_out, df_assign, file):
    df_file = pd.read_csv(os.path.join(path_data, file), sep=' ')
    cells = df_file.columns
    try:
        sub_assign = df_assign.loc[cells, :]
    except KeyError:
        print(file, "   not exist in cellassignment file")
        return
    sub_assign = sub_assign.dropna()

    tissue = info_folder.loc[info_folder['file'] == file, 'tissue'].tolist()[0]
    folder = info_folder.loc[info_folder['file'] == file, 'folder'].tolist()[0]
    out_folder = os.path.join(path_out, folder)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    df_file.to_csv(os.path.join(out_folder, tissue + '.txt'), sep='\t')
    sub_assign.to_csv(
        os.path.join(out_folder, 'CellAssignments_' + tissue + '.txt'),
        sep='\t')

    return


def main(path_data, file_assign, path_out):
    files = os.listdir(path_data)
    tissue = []
    folder = []
    tissue_pattern = re.compile(r'.+_rm')
    for file in files:
        post_tissue = tissue_pattern.search(file).group()[:-3]
        tissue.append(post_tissue)
        organ = post_tissue.strip().split('.')[0]
        if organ[-1].isdigit():
            organ = organ[:-1]
        folder.append(organ)

    info_folder = pd.DataFrame({'file': files, 'tissue': tissue,
                                'folder': folder})
    info_folder.to_csv(os.path.join(path_out, 'folder.info'), sep='\t',
                       index=None)

    df_assign = pd.read_csv(file_assign, sep=',', index_col=0)
    df_assign.index = df_assign['Cell.name']

    pool = Pool(processes=10)
    func = partial(trans_sep, info_folder, path_data, path_out, df_assign)
    pool.map(func, files)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    main(
        '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/rmbatch_dge',
        '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/'
        'MCA_CellAssignments.csv',
        '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/MCA_by_tissue')
    time_end = time()
    print(time_end - time_start)