#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: atac_map.py
# @time: 8/13/19 9:11 PM

from time import time
import os
import re
from multiprocessing import Pool
from functools import partial
from collections import defaultdict


def func_align(path_out, file_fq):
    pre_pattern = re.compile(r'/[^/]+\.')
    pre_nm = pre_pattern.search(file_fq).group()[1:-1]
    file_sam = os.path.join(path_out, f'{pre_nm}.sam')
    file_bam = os.path.join(path_out, f'{pre_nm}.bam')
    file_sort = os.path.join(path_out, f'{pre_nm}.sorted')
    file_sortbam = os.path.join(path_out, f'{pre_nm}.sorted.bam')
    file_sort_dedup = os.path.join(path_out, f'{pre_nm}.sorted.dedup.bam')
    os.system(f"bowtie2 -p 6 -x /home/genomewide/refgenome/hg38/hg38 "
              f"{file_fq} -S {file_sam}")
    os.system(f"/home/disk/liupengfei/software/samtools-1.5/samtools "
              f"view -bS {file_sam} > {file_bam}")
    os.system(f"/home/disk/liupengfei/software/samtools-1.5/samtools "
              f"sort {file_bam} -T {file_sort} -o {file_sortbam}")
    os.system(f"/home/disk/liupengfei/software/samtools-1.5/samtools "
              f"rmdup -s {file_sortbam} {file_sort_dedup}")
    # os.remove(file_sam)
    # os.remove(file_bam)
    # os.remove(file_sortbam)

    return


def align(dir_fastq, path_out):
    files = os.listdir(dir_fastq)
    files_fq = []
    for file in files:
        if file.endswith('.fastq'):
            files_fq.append(os.path.join(dir_fastq, file))

    pool = Pool(processes=5)
    func = partial(func_align, path_out)
    pool.map(func, files_fq)
    pool.close()

    return


def func_map(dir_bam, dict_bam, ref_peaks, dir_count, sample):
    file_merge = os.path.join(dir_bam, f"{sample}.merge.bam")
    os.system(f"/home/disk/liupengfei/software/samtools-1.5/samtools "
              f"merge {file_merge} {' '.join(dict_bam[sample])}")

    file_count = os.path.join(dir_count, f"{sample}_counts.bed")
    os.system(f"bedtools coverage -abam {file_merge} -b {ref_peaks} "
              f"> {file_count}")

    return


def map_refpeaks(dir_bam, dir_count, ref_peaks):
    files = os.listdir(dir_bam)
    dict_bam = defaultdict(list)
    sample_pattern = re.compile(r'.+_L')
    for file in files:
        if file.endswith('.sorted.dedup.bam'):
            sample = sample_pattern.search(file).group()[:-2]
            dict_bam[sample].append(os.path.join(dir_bam, file))

    pool = Pool(processes=5)
    func = partial(func_map, dir_bam, dict_bam, ref_peaks, dir_count)
    pool.map(func, dict_bam.keys())
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    dir_fastq = '/home/yzj/ZF/0804'
    dir_bam = '/home/zy/ATAC_RNAseq_Survival/ATAC_map/bams'
    dir_count = '/home/zy/ATAC_RNAseq_Survival/ATAC_map/count_files'
    ref_peaks = '/home/zy/ATAC_RNAseq_Survival/map_bam_peaks/ATAC_peaks.txt'
    # align(dir_fastq, dir_bam)
    map_refpeaks(dir_bam, dir_count, ref_peaks)

    time_end = time()
    print(time_end - time_start)
