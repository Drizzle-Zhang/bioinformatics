#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: annotate_DHS.py
# @time: 10/28/19 4:08 PM

from time import time
import os
from multiprocessing import Pool
from functools import partial


def sub_stan(type_bed, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    with open(file_in, 'r') as r_f:
        with open(os.path.join(path_out, file)) as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{label}\t{score}\t.\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                score = list_line[6]
                w_f.write(fmt.format(**locals()))

    return


def standardize_bed(path_in, path_out, type_bed):
    types = os.listdir(path_in)
    list_input = []
    for biotype in types:
        path_type = os.path.join(path_in, biotype)
        path_type_out = os.path.join(path_out, biotype)
        if not os.path.exists(path_type_out):
            os.makedirs(path_type_out)
        if os.path.isdir(path_type):
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
                                           path_in=path_life))

    pool = Pool(processes=40)
    func_stan = partial(sub_stan, type_bed)
    pool.map(func_stan, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
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

    time_end = time()
    print(time_end - time_start)
