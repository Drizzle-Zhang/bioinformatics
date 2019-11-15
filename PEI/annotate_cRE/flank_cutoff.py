#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: flank_cutoff.py
# @time: 11/13/19 8:22 PM

from time import time
from prepare_bed_file_mu02 import merge_bed
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from functools import partial
from subprocess import check_output


def count_rows(dict_in):
    flank_percent = dict_in['flank_percent']
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    split_out = os.path.join(path_out, f"{term_name}.bed")
    output = check_output(f"wc -l {split_out}", shell=True)
    num_rows = int(str(output).split(' ')[0][2:])
    label = path_out.split('/')[-1]
    df_bed = pd.read_csv(split_out, sep='\t', usecols=[1, 2])
    length = df_bed.iloc[:, 1] - df_bed.iloc[:, 0]
    medium = length.quantile(0.5)
    percent_75 = length.quantile(0.75)
    percent_90 = length.quantile(0.9)
    percent_95 = length.quantile(0.95)
    percent_99 = length.quantile(0.99)
    maximum = np.max(length)

    return {'flank_percent': flank_percent, 'num_rows': num_rows,
            'label': label, 'medium': medium, 'percent_75': percent_75,
            'percent_90': percent_90, 'percent_95': percent_95,
            'percent_99': percent_99, 'maximum': maximum}


def generate_flank_plot_file(path_in, path_flank, list_select):
    if os.path.exists(path_flank):
        os.system(f"rm -rf {path_flank}")
    os.mkdir(path_flank)

    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')
    list_input = []
    array_flank = np.linspace(0, 1, 51)
    for line in list_select:
        sub_meta = df_meta.loc[
                   (df_meta['Biosample term name'] == line[0]) &
                   (df_meta['Biosample life stage'] == line[1]) &
                   (df_meta['Biosample organ'].apply(
                         lambda x: line[2] in x.strip().split(','))), :]
        path_out = os.path.join(
            path_flank, '_'.join([val.replace(' ', '-') for val in line]))
        os.mkdir(path_out)
        accession_ids = sub_meta['File accession'].tolist()
        for flank in array_flank:
            term_name = f'flank_{str(flank)}'
            list_input.append(
                dict(path=path_out,
                     term_name=term_name,
                     accession_ids=accession_ids,
                     flank_percent=flank))

    pool = Pool(processes=80)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

    pool = Pool(processes=20)
    list_out = pool.map(count_rows, list_input)
    pool.close()
    df_out = pd.DataFrame(list_out)
    df_out.to_csv(
        os.path.join(path_flank, 'flank_count.txt'), sep='\t', index=None
    )

    return


if __name__ == '__main__':
    time_start = time()
    list_test = [['brain', 'embryonic', 'brain'],
                 ['frontal cortex', 'adult', 'brain'],
                 ['sigmoid colon', 'adult', 'intestine'],
                 ['B cell', 'adult', 'blood'],
                 ['lung', 'embryonic', 'lung'],
                 ['heart', 'child', 'heart']]
    # DHS
    path_dhs = '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
               'DNase-seq/GRCh38tohg19/'
    path_dhs_flank = '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
                     'DNase-seq/GRCh38tohg19/flank'
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_test)

    list_test = [['stomach', 'adult', 'stomach'],
                 ['layer of hippocampus', 'adult', 'brain'],
                 ['sigmoid colon', 'adult', 'intestine'],
                 ['B cell', 'adult', 'blood'],
                 ['placenta', 'embryonic', 'extraembryonic component'],
                 ['heart left ventricle', 'adult', 'heart']]
    # H3K4me3
    path_h3k4me3 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3'
    path_h3k4me3_flank = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3/flank'
    generate_flank_plot_file(path_h3k4me3, path_h3k4me3_flank, list_test)

    # H3K27ac
    path_h3k27ac = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac'
    path_h3k27ac_flank = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac/flank'
    generate_flank_plot_file(path_h3k27ac, path_h3k27ac_flank, list_test)

    time_end = time()
    print(time_end - time_start)
