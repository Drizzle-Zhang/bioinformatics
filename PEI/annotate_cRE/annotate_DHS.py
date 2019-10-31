#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: annotate_DHS.py
# @time: 10/28/19 4:08 PM

from time import time
import os
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np


def generate_file_list(path_in, path_out):
    folder_1 = os.listdir(path_in)
    list_input = []
    for element_1 in folder_1:
        path_1 = os.path.join(path_in, element_1)
        if os.path.isdir(path_type):
            path_type_out = os.path.join(path_out, biotype)
            if not os.path.exists(path_type_out):
                os.makedirs(path_type_out)
            life_stages = os.listdir(path_type)
            for life_stage in life_stages:
                path_life = os.path.join(path_type, life_stage)
                path_life_out = os.path.join(path_type_out, life_stage)
                if not os.path.exists(path_life_out):
                    os.makedirs(path_life_out)
                files = os.listdir(path_life)
                for file in files:
                    list_input.append(dict(path_out=path_life_out,
                                           file=file,
                                           path_in=path_life,
                                           biotype=biotype,
                                           life_stage=life_stage,
                                           label=f'{biotype}|{life_stage}|'
                                                 f'{file}'))

    return list_input


def sub_stan(type_bed, dict_in, col_score=6):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    with open(file_in, 'r') as r_f:
        with open(os.path.join(path_out, file), 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{label}\t{score}\t.\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                score = max(list_line[col_score].strip().split(','))
                w_f.write(fmt.format(**locals()))

    return


def standardize_bed(path_in, path_out, type_bed):
    list_input = generate_file_list(path_in, path_out)
    pool = Pool(processes=40)
    func_stan = partial(sub_stan, type_bed)
    pool.map(func_stan, list_input)
    pool.close()

    return


def annotate_dhs_promoter(path_promoter_in, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    os.system(f"bedtools intersect -a {file_in} -b {path_promoter_in} -loj "
              f"> {bedtools_out}")
    os.system(f"cut -f 1,2,3,4,5,6,10,11 {bedtools_out} > "
              f"{os.path.join(path_out, file)}")
    # os.remove(bedtools_out)

    return


def annotate_dhs_histone(dict_in):
    path_out = dict_in['path_out']
    path_in_1 = dict_in['path_in_1']
    path_in_2 = dict_in['path_in_2']
    file = dict_in['file']
    file_in_1 = os.path.join(path_in_1, file)
    file_in_2 = os.path.join(path_in_2, file)
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    os.system(f"bedtools intersect -a {file_in_1} -b {file_in_2} -loj "
              f"> {bedtools_out}")
    os.system(f"cut -f 1,2,3,4,5,6,7,8,12,13 {bedtools_out} > "
              f"{os.path.join(path_out, file)}")
    # os.remove(bedtools_out)

    return


def annotate_dhs(path_dhs_in, path_promoter_in, path_h3k27ac_in, path_dhs_out):
    tmp_path_out = path_dhs_out + '_tmp'
    if not os.path.exists(tmp_path_out):
        os.mkdir(tmp_path_out)

    list_dhs = generate_file_list(path_dhs_in, tmp_path_out)
    pool = Pool(processes=40)
    func_pro = partial(annotate_dhs_promoter, path_promoter_in)
    pool.map(func_pro, list_dhs)
    pool.close()

    list_dhs_pro = generate_file_list(tmp_path_out, path_dhs_out)
    list_h3k27ac = generate_file_list(path_h3k27ac_in, path_dhs_out)
    df_dhs_pro = pd.DataFrame(list_dhs_pro,
                              columns=['path_out', 'file', 'path_in',
                                       'biotype', 'life_stage', 'label'])
    df_dhs_pro.columns = ['path_out', 'file', 'path_in_1',
                          'biotype', 'life_stage', 'label']
    df_h3k27ac = pd.DataFrame(list_h3k27ac,
                              columns=['path_out', 'file', 'path_in',
                                       'biotype', 'life_stage', 'label'])
    df_h3k27ac.columns = ['path_out', 'file', 'path_in_2',
                          'biotype', 'life_stage', 'label']
    df_merge = pd.merge(df_dhs_pro, df_h3k27ac,
                        on=['path_out', 'file',
                            'biotype', 'life_stage', 'label'], how='inner')
    list_dict = df_merge.to_dict('records')

    pool = Pool(processes=40)
    pool.map(annotate_dhs_histone, list_dict)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # build DHS reference by organ

    # standardization
    # DHS
    path_dhs = '/home/zy/driver_mutation/data/DHS/hg19'
    path_dhs_stan = '/home/zy/driver_mutation/data/DHS/hg19_standard'
    standardize_bed(path_dhs, path_dhs_stan, 'DHS')

    # H3K27ac
    path_h3k27ac = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/hg19/H3K27ac'
    path_h3k27ac_stan = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/' \
        'hg19/H3K27ac_standard'
    standardize_bed(path_h3k27ac, path_h3k27ac_stan, 'H3K27ac')

    # annotate DHS
    path_promoter = '/home/zy/driver_mutation/data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'
    path_anno = '/home/zy/driver_mutation/data/DHS/hg19_annotation'
    annotate_dhs(path_dhs_stan, path_promoter, path_h3k27ac_stan, path_anno)

    time_end = time()
    print(time_end - time_start)
