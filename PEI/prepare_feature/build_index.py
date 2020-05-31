#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: build_index.py
# @time: 2020/4/28 21:45

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from subprocess import check_output
import sys
# sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
sys.path.append(
    '/lustre/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE')
from preparation import merge_standard_bed


def merge_cell_tissue():
    dict_merge = [dict(
        path=path_dhs_merge,
        term_name='all_cellline_tissue',
        accession_ids=
        ['cell_line/DHS/GRCh38tohg19_standard/all_celllines',
         'tissue/DHS/GRCh38tohg19_cluster/all_organs'],
        flank_percent=1.0)]
    merge_standard_bed(path_mid, dict_merge, num_cpu)

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


def calculate_score(df_in):
    if df_in.shape[0] == 1:
        return df_in
    else:
        max_score = np.max(df_in.loc[:, 1])
        row_out = df_in.loc[df_in.loc[:, 1] == max_score, :]
        return row_out


def sub_generate_index(dict_in):
    str_term = dict_in['str_term']
    path_term = dict_in['path_term']
    file_ref = dict_in['file_ref']

    file_term = os.path.join(path_term, str_term + '.bed')
    file_intersect = os.path.join(path_term, 'intersect.txt')
    file_index = os.path.join(path_term, 'index.txt')
    file_ref_index = os.path.join(path_term, 'ref_index.txt')
    # if os.path.isfile(file_index):
    #     return

    os.system(f"bedtools intersect -a {file_term} -b {file_ref} -wao "
              f"-sorted | cut -f 4,5,12,13,14 > {file_intersect}")
    len_term = int(str(check_output(f"wc -l {file_term}",
                                    shell=True).strip()).split(' ')[0][2:])
    len_intersect = int(str(check_output(
        f"wc -l {file_intersect}", shell=True).strip()).split(' ')[0][2:])
    if len_term == len_intersect:
        os.system(f"cut -f 1,2,3,4,5 {file_intersect} > {file_index}")
    else:
        df_intersect = pd.read_csv(file_intersect, sep='\t', header=None)
        df_0 = df_intersect.loc[df_intersect[4] == 0, :]
        if df_0.shape[0] > 0:
            print('Error: ' + str_term)
        df_pn = (df_intersect.loc[df_intersect[4] > 0, :]).copy()
        df_pn['key'] = df_pn[0]
        df_pn_uniq = df_pn.groupby('key').apply(drop_dup)
        df_pn_uniq = df_pn_uniq.drop(['key'], axis=1)
        df_pn_uniq = df_pn_uniq.drop_duplicates(subset=0)
        # assert df_pn_uniq.shape[0] == len_term
        if df_pn_uniq.shape[0] != len_term:
            print(path_term)
            # df_bed = pd.read_csv(file_term, sep='\t', header=None)
            # list_bed = df_bed[3].tolist()
            # set_bed=set(df_bed[3].tolist())
            # for dhs_id in set_bed:
            #     count_id = list_bed.count(dhs_id)
            #     if count_id > 1:
            #         print(dhs_id)
        df_pn_uniq.to_csv(file_index, sep='\t', header=None, index=None)

        df_pn_uniq = df_pn_uniq.iloc[:, 1:].copy()
        df_pn_uniq['key1'] = df_pn_uniq[2]
        df_ref_uniq = df_pn_uniq.groupby('key1').apply(calculate_score)
        df_ref_uniq = df_ref_uniq.drop(['key1', 4], axis=1)
        df_ref_uniq = df_ref_uniq.drop_duplicates()
        df_ref_uniq.to_csv(file_ref_index, sep='\t', header=None, index=None)

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

    life_organs = list(set(df_meta_tissue['Biosample life_organ'].tolist()))
    for life_organ in life_organs:
        str_life_organ = life_organ.replace(' ', '_')
        path_term = os.path.join(path_dhs_tissue_cluster, str_life_organ)
        list_input.append(
            {'str_term': str_life_organ, 'path_term': path_term,
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

    pool = Pool(processes=num_cpu)
    pool.map(sub_generate_index, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    path_dhs_cell = \
        path_mid + '/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_cluster = \
        path_mid + '/tissue/DHS/GRCh38tohg19_cluster'
    path_dhs_tissue_stan = \
        path_mid + '/tissue/DHS/GRCh38tohg19_standard'
    path_dhs_merge = path_mid + '/database_feature/DHS_index'
    meta_suborgan = path_origin + '/meta_file/meta.reference.tsv'

    merge_cell_tissue()
    print(
        'Bed file incorperating cell line and tissue data has been generated')

    generate_index_file()

    time_end = time()
    print(time_end - time_start)
