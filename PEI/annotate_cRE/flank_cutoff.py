#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: flank_cutoff.py
# @time: 11/13/19 8:22 PM

from time import time
from node10_preparation import merge_bed, merge_standard_bed
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


def generate_flank_plot_file(path_in, path_flank, list_dict, merge_type=1):
    if os.path.exists(path_flank):
        os.system(f"rm -rf {path_flank}")
    os.mkdir(path_flank)

    list_input = []
    array_flank = np.linspace(0, 1.5, 76)
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

    if merge_type == 1:
        pool = Pool(processes=40)
        func_merge = partial(merge_bed, path_in)
        pool.map(func_merge, list_input)
        pool.close()
    elif merge_type == 2:
        pool = Pool(processes=40)
        func_merge = partial(merge_standard_bed, path_in)
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


def generate_organ_input(path_in, list_select):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'meta.reference.tsv'), sep='\t')

    list_out = []
    for line in list_select:
        sub_meta = df_meta.loc[df_meta['Biosample organ'] == line, :]
        sub_name = line.replace(' ', '_')
        accession_ids = []
        for sub_dict in sub_meta.to_dict('records'):
            organ = sub_dict['Biosample organ']
            life_stage = sub_dict['Biosample life stage']
            term = sub_dict['Biosample term name']
            term_name = \
                term.replace(' ', '_').replace('/', '+').replace("'", "--")
            accession_ids.append(
                f"{organ.replace(' ', '_')}/{life_stage.replace(' ', '_')}/"
                f"{term_name}/{term_name}"
            )
        list_out.append({'name': sub_name, 'files': accession_ids})

    return list_out


def generate_all_input(path_in):
    folders = os.listdir(path_in)

    accession_ids = []
    for folder in folders:
        if os.path.isdir(os.path.join(path_in, folder)) & (folder != 'flank'):
            accession_ids.append(f"{folder}/{folder}")
    list_out = {'name': 'all_organs', 'files': accession_ids}

    return list_out


if __name__ == '__main__':
    time_start = time()
    # experiments
    list_test = ['ENCSR000EIJ', 'ENCSR503HIB', 'ENCSR792ZXA', 'ENCSR946DXB',
                 'ENCSR095GWE', 'ENCSR499IFY']
    path_dhs = '/local/zy/PEI/data/ENCODE/DNase-seq/GRCh38tohg19/'
    path_dhs_flank = '/local/zy/PEI/data/ENCODE/' \
                     'DNase-seq/GRCh38tohg19/flank'
    list_dict_dhs = generate_experiment_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs)

    # term
    list_test = [['brain', 'embryonic', 'brain'],
                 ['frontal cortex', 'adult', 'brain'],
                 ['sigmoid colon', 'adult', 'intestine'],
                 ['limb', 'embryonic', 'limb'],
                 ['lung', 'embryonic', 'lung'],
                 ['cerebellar cortex', 'adult', 'brain']]
    path_dhs = '/local/zy/PEI/data/ENCODE/DNase-seq/' \
               'GRCh38tohg19_experiment'
    path_dhs_flank = '/local/zy/PEI/data/ENCODE/' \
                     'DNase-seq/GRCh38tohg19_experiment/flank'
    list_dict_dhs = generate_term_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs)

    # organ
    list_test = ['adrenal gland', 'brain', 'intestine', 'musculature of body',
                 'testis', 'heart']
    path_dhs = '/local/zy/PEI/data/DHS/GRCh38tohg19_standard'
    path_dhs_flank = '/local/zy/PEI/data/DHS/GRCh38tohg19_standard/flank'
    list_dict_dhs = generate_organ_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs, 2)

    # all
    path_dhs = '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster'
    path_dhs_flank = '/local/zy/PEI/data/DHS/GRCh38tohg19_standard/flank'
    list_dict_dhs = generate_organ_input(path_dhs, list_test)
    generate_flank_plot_file(path_dhs, path_dhs_flank, list_dict_dhs, 2)

    time_end = time()
    print(time_end - time_start)
