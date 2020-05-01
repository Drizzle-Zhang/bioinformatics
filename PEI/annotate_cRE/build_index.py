#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: build_index.py
# @time: 2020/4/28 21:45

from time import time
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from multiprocessing import Pool
from subprocess import check_output
import sys
sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
from node10_preparation import merge_standard_bed


def merge_cell_tissue():
    dict_merge = dict(
        path=path_dhs_merge,
        term_name='all_cellline_tissue',
        accession_ids=
        ['cell_line/DHS/GRCh38tohg19_standard/all_celllines',
         'tissue/DHS/GRCh38tohg19_cluster/all_organs'],
        flank_percent=1.0)
    merge_standard_bed('/local/zy/PEI/mid_data', dict_merge)

    file_merge = os.path.join(path_dhs_merge, 'all_cellline_tissue.bed')
    file_merge_index = os.path.join(path_dhs_merge, 'all_index.txt')
    with open(file_merge_index, 'w') as w_index:
        fmt = "{chrom}\t{start}\t{end}\t{dhs_id}\t{index}\n"
        with open(file_merge, 'r') as r_merge:
            for index, line in enumerate(r_merge):
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                w_index.write(fmt.format(**locals()))

    return


def drop_dup(x):
    if x.shape[0] == 1:
        return x
    else:
        max_overlap = np.max(x.iloc[:, -1])
        row_out = x.loc[x.iloc[:, -1] == max_overlap, :]
        return row_out


def sub_generate_index(dict_in):
    str_term = dict_in['str_term']
    path_term = dict_in['path_term']
    file_ref = dict_in['file_ref']

    file_term = os.path.join(path_term, str_term + '.bed')
    file_intersect = os.path.join(path_term, 'intersect.txt')
    file_index = os.path.join(path_term, 'index.txt')

    os.system(f"bedtools intersect -a {file_term} -b {file_ref} -wao | "
              f"cut -f 4,12,13,14 > {file_intersect}")
    len_term = int(str(check_output(f"wc -l {file_term}",
                                    shell=True).strip()).split(' ')[0][2:])
    len_intersect = int(str(check_output(
        f"wc -l {file_intersect}", shell=True).strip()).split(' ')[0][2:])
    if len_term == len_intersect:
        os.system(f"cut -f 1,2,3 {file_intersect} > {file_index}")
    else:
        df_intersect = pd.read_csv(file_intersect, sep='\t', header=None)
        df_0 = df_intersect.loc[df_intersect[3] == 0, :]
        if df_0.shape[0] > 0:
            print('Error: ' + str_term)
        df_pn = (df_intersect.loc[df_intersect[3] > 0, :]).copy()
        df_pn['key'] = df_pn[0]
        df_pn_uniq = df_pn.groupby('key').apply(drop_dup)
        df_pn_uniq = df_pn_uniq.drop(['key', 3], axis=1)
        df_pn_uniq = df_pn_uniq.drop_duplicates(subset=0)
        assert df_pn_uniq.shape[0] == len_term
        df_pn_uniq.to_csv(file_index, sep='\t', header=None, index=None)

    os.remove(file_intersect)

    return


def generate_index_file():
    file_ref = os.path.join(path_dhs_merge, 'all_index.txt')
    list_input = []
    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_dhs_cell, 'meta.reference.tsv'), sep='\t')
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_dhs_cell, str_term)
        list_input.append({'str_term': str_term, 'path_term': path_term,
                           'file_ref': file_ref})

    # tissue
    df_meta_tissue = pd.read_csv(
        os.path.join(path_dhs_tissue_stan, 'meta.reference.tsv'), sep='\t')
    for i in range(df_meta_tissue.shape[0]):
        str_life_organ = \
            (df_meta_tissue.loc[i, 'Biosample life_organ']).replace(' ', '_')
        str_term = (df_meta_tissue.loc[i, 'Biosample term name']).replace(
            ' ', '_').replace('/', '+').replace("'", "--")
        path_term = os.path.join(
            path_dhs_tissue_stan, f"{str_life_organ}/{str_term}")
        list_input.append(
            {'str_term': str_term, 'path_term': path_term,
             'file_ref': file_ref})

    # tissue cluster
    df_meta_cluster = pd.read_csv(meta_suborgan, sep='\t')
    df_meta_cluster = df_meta_cluster.dropna()
    df_merge = pd.merge(df_meta_cluster, df_meta_tissue,
                        on=['Biosample life stage', 'Biosample organ',
                            'Biosample term name'])
    df_merge = df_merge.loc[:, ['Biosample life_organ', 'Biosample suborgan']]
    df_merge = df_merge.drop_duplicates()
    for i in df_merge.index:
        str_life_organ = \
            (df_merge.loc[i, 'Biosample life_organ']).replace(' ', '_')
        str_term = (df_merge.loc[i, 'Biosample suborgan']).replace(
            ' ', '_').replace('/', '+').replace("'", "--")
        path_term = os.path.join(
            path_dhs_tissue_cluster, f"{str_life_organ}/{str_term}")
        list_input.append(
            {'str_term': str_term, 'path_term': path_term,
             'file_ref': file_ref})

    pool = Pool(processes=40)
    pool.map(sub_generate_index, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    path_dhs_cell = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_cluster = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster'
    path_dhs_tissue_stan = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_standard'
    path_dhs_merge = '/local/zy/PEI/mid_data/database_feature/DHS_index'
    meta_suborgan = '/local/zy/PEI/origin_data/meta_file/meta.reference.tsv'

    merge_cell_tissue()

    generate_index_file()

    time_end = time()
    print(time_end - time_start)