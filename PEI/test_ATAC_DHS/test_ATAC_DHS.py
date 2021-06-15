#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: test_ATAC_DHS.py
# @time: 2020/4/22 15:17

from time import time
import os
import sys
from multiprocessing import Pool
from functools import partial
from subprocess import check_output
import pandas as pd
import numpy as np
from collections import defaultdict
sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
from node10_preparation_cellline import sub_hg38tohg19


def calculate_overlap(file_ref, path_hg19, dict_in):
    sub_file = os.path.join(path_hg19, dict_in['File accession'] + '.bed')
    file_intersect = os.path.join(
        path_hg19, dict_in['File accession'] + '.intersect')

    os.system(f"bedtools intersect -wa -wb -a {sub_file} -b {file_ref} > "
              f"{file_intersect}")

    sub_num = int(str(check_output(f"wc -l {sub_file}",
                                   shell=True).strip()).split(' ')[0][2:])
    intersect_num = int(str(check_output(
        f"cut -f 1,2,3 {file_intersect} | sort | uniq | wc -l",
        shell=True).strip()).split(' ')[0][2:-1])

    return {'File accession': dict_in['File accession'],
            'overlap_ratio': intersect_num/sub_num}


def test_atac_dhs(file_ref, path_hg19, path_hg38):
    list_input = []
    files_hg38 = os.listdir(path_hg38)
    for file in files_hg38:
        list_input.append({'File accession': file[:-4], 'Assembly': 'GRCh38'})

    pool = Pool(processes=10)
    func_hg38tohg19 = partial(sub_hg38tohg19, path_hg38, path_hg19)
    pool.map(func_hg38tohg19, list_input)
    pool.close()

    pool = Pool(processes=10)
    func_calc = partial(calculate_overlap, file_ref, path_hg19)
    result = pool.map(func_calc, list_input)
    pool.close()

    df_res = pd.DataFrame(result)
    df_res.to_csv(os.path.join(path_hg19, 'results.tsv'), sep='\t', index=None)

    return


def test_brain(sample):
    term = sample.split('_')[0]
    sample_id = sample.split('_')[1]
    file_sample = os.path.join(path_root, f"{term}/{sample_id}.bed")
    file_sample_tmp = file_sample + '.tmp'

    os.system(f"grep -w {sample} {file_all_tar} > {file_sample_tmp}")
    with open(file_sample, 'w') as w_bed:
        fmt = "{chrom}\t{start}\t{end}\t{region_id}\t{score}\t.\n"
        with open(file_sample_tmp, 'r') as r_sample:
            for line_sample in r_sample:
                list_line_sample = line_sample.strip().split(' ')
                region_id = list_line_sample[0]
                score = list_line_sample[1]
                if float(score) < 1.64:
                    continue
                list_id_region = dict_id_region[region_id]
                if list_id_region:
                    list_id_region = list_id_region[0]
                    chrom = list_id_region[0]
                    start = list_id_region[1]
                    end = list_id_region[2]
                    w_bed.write(fmt.format(**locals()))

    if term == 'brainCortex':
        file_ref = file_ref_cortex
    elif term == 'brainCerebellum':
        file_ref = file_ref_cerebellum
    else:
        return
    file_intersect = file_sample + '.intersect'
    os.system(f"bedtools intersect -wa -wb -a {file_sample} -b {file_ref} > "
              f"{file_intersect}")

    sub_num = int(str(check_output(f"wc -l {file_sample}",
                                   shell=True).strip()).split(' ')[0][2:])
    intersect_num = int(str(check_output(
        f"cut -f 1,2,3 {file_intersect} | sort | uniq | wc -l",
        shell=True).strip()).split(' ')[0][2:-1])

    return {'File accession': sample,
            'overlap_ratio': intersect_num/sub_num}


if __name__ == '__main__':
    time_start = time()
    # A549
    file_ref_a549 = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19/A549/A549.bed'
    path_hg19_a549 = '/local/zy/PEI/test_ATAC_DHS/A549/ATAC/hg19'
    path_hg38_a549 = '/local/zy/PEI/test_ATAC_DHS/A549/ATAC/GRCh38'
    test_atac_dhs(file_ref_a549, path_hg19_a549, path_hg38_a549)

    # brain
    file_ref_cortex = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster/' \
        'adult_brain/cerebral_cortex/cerebral_cortex.bed'
    file_ref_cerebellum = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster/' \
        'adult_brain/cerebellum/cerebellum.bed'
    path_root = '/local/zy/PEI/test_ATAC_DHS/brain'
    file_samples = os.path.join(path_root, 'brain_samples.txt')
    file_ref_regions = \
        os.path.join(path_root, 'PIP-04_all_TARs.70pc.active.bed')
    file_all_tar = os.path.join(path_root, 'brain_TARs.txt')
    # build dict
    dict_id_region = defaultdict(list)
    with open(file_ref_regions, 'r') as r_regions:
        for line in r_regions:
            list_line = line.strip().split('\t')
            dict_id_region[list_line[3]].append(list_line[0:3])

    samples = (np.loadtxt(file_samples, dtype='str')).tolist()
    pool = Pool(processes=40)
    result = pool.map(test_brain, samples)
    pool.close()

    df_res = pd.DataFrame(result)
    df_res.to_csv(os.path.join(path_root, 'results.tsv'), sep='\t', index=None)

    time_end = time()
    print(time_end - time_start)
