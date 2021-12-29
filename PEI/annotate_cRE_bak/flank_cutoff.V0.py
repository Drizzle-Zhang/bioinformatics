#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: flank_cutoff.py
# @time: 11/13/19 8:22 PM

from time import time
from prepare_bed_file import merge_bed
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
    df_bed = pd.read_csv(split_out, sep='\t', usecols=[1, 2], header=None)
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


def generate_flank_plot_file(path_in, path_flank, list_dict):
    if os.path.exists(path_flank):
        os.system(f"rm -rf {path_flank}")
    os.mkdir(path_flank)

    list_input = []
    array_flank = np.linspace(0, 1, 51)
    for line in list_dict:
        path_out = os.path.join(path_flank, line['name'])
        os.mkdir(path_out)
        accession_ids = line['files']
        for flank in array_flank:
            term_name = f'flank_{str(flank)}'
            list_input.append(
                dict(path=path_out,
                     term_name=term_name,
                     accession_ids=accession_ids,
                     flank_percent=flank))

    pool = Pool(processes=40)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

    pool = Pool(processes=40)
    list_out = pool.map(count_rows, list_input)
    pool.close()
    df_out = pd.DataFrame(list_out)
    df_out.to_csv(
        os.path.join(path_flank, 'flank_count.txt'), sep='\t', index=None
    )

    return


def generate_experiment_input(path_in, list_select):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    list_out = []
    for line in list_select:
        sub_meta = df_meta.loc[
                   df_meta['Experiment accession'] == line, :]
        sub_name = line
        accession_ids = sub_meta['File accession'].tolist()
        list_out.append({'name': sub_name, 'files': accession_ids})

    return list_out


def generate_term_input(path_in, list_select):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    list_out = []
    for line in list_select:
        sub_meta = df_meta.loc[
                   (df_meta['Biosample term name'] == line[0]) &
                   (df_meta['Biosample life stage'] == line[1]) &
                   (df_meta['Biosample organ'].apply(
                         lambda x: line[2] in x.strip().split(','))), :]
        sub_name = '_'.join([val.replace(' ', '-') for val in line])
        accession_ids = sub_meta['Experiment accession'].tolist()
        list_out.append({'name': sub_name, 'files': accession_ids})

    return list_out


if __name__ == '__main__':
    time_start = time()
    # experiments
    # DHS
    list_test = ['ENCSR000EIJ', 'ENCSR503HIB', 'ENCSR792ZXA', 'ENCSR946DXB',
                 'ENCSR095GWE', 'ENCSR499IFY']
    path_dhs = '/home/zy/driver_mutation/data/ENCODE/DNase-seq/GRCh38tohg19/'
    path_dhs_flank = '/home/zy/driver_mutation/data/ENCODE/' \
                     'DNase-seq/GRCh38tohg19/flank'
    list_dict_dhs = generate_experiment_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs)

    # # H3K4me3
    # list_test = ['ENCSR432GOP', 'ENCSR477BHF', 'ENCSR377ILM', 'ENCSR813ZEY',
    #              'ENCSR780FXX', 'ENCSR107RDP']
    # path_h3k4me3 = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K4me3'
    # path_h3k4me3_flank = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K4me3/flank'
    # list_dict_h3k4me3 = generate_experiment_input(path_h3k4me3, list_test)
    # generate_flank_plot_file(
    #     path_h3k4me3, path_h3k4me3_flank, list_dict_h3k4me3)
    #
    # # H3K27ac
    # list_test = ['ENCSR743DDX', 'ENCSR668GBL', 'ENCSR203KCB', 'ENCSR726WVB',
    #              'ENCSR150QXE', 'ENCSR792VLP']
    # path_h3k27ac = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K27ac'
    # path_h3k27ac_flank = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K27ac/flank'
    # list_dict_h3k27ac = generate_experiment_input(path_h3k27ac, list_test)
    # generate_flank_plot_file(
    #     path_h3k27ac, path_h3k27ac_flank, list_dict_h3k27ac)

    # DHS
    list_test = [['brain', 'embryonic', 'brain'],
                 ['frontal cortex', 'adult', 'brain'],
                 ['sigmoid colon', 'adult', 'intestine'],
                 ['limb', 'embryonic', 'limb'],
                 ['lung', 'embryonic', 'lung'],
                 ['cerebellar cortex', 'adult', 'brain']]
    path_dhs = '/home/zy/driver_mutation/data/ENCODE/DNase-seq/' \
               'GRCh38tohg19_experiment'
    path_dhs_flank = '/home/zy/driver_mutation/data/ENCODE/' \
                     'DNase-seq/GRCh38tohg19_experiment/flank'
    list_dict_dhs = generate_term_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs)

    # list_test = [['stomach', 'adult', 'stomach'],
    #              ['layer of hippocampus', 'adult', 'brain'],
    #              ['sigmoid colon', 'adult', 'intestine'],
    #              ['adrenal gland', 'adult', 'adrenal gland'],
    #              ['placenta', 'embryonic', 'extraembryonic component'],
    #              ['heart left ventricle', 'adult', 'heart']]
    # # H3K4me3
    # path_h3k4me3 = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K4me3_experiment'
    # path_h3k4me3_flank = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K4me3_experiment/flank'
    # list_dict_h3k4me3 = generate_term_input(path_h3k4me3, list_test)
    # generate_flank_plot_file(
    #     path_h3k4me3, path_h3k4me3_flank, list_dict_h3k4me3)
    #
    # # H3K27ac
    # path_h3k27ac = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K27ac_experiment'
    # path_h3k27ac_flank = \
    #     '/home/zy/driver_mutation/data/ENCODE/' \
    #     'histone_ChIP-seq/GRCh38tohg19/H3K27ac_experiment/flank'
    # list_dict_h3k27ac = generate_term_input(path_h3k27ac, list_test)
    # generate_flank_plot_file(
    #     path_h3k27ac, path_h3k27ac_flank, list_dict_h3k27ac)

    time_end = time()
    print(time_end - time_start)
