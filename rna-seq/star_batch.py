#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: star_batch.py
# @time: 2018/10/23 21:51

import sys
from time import time
import os
import re
import subprocess


def star(path_trim_fq, path_index, path_output, process=2):
    # read fastq files of paired-end sequencing
    fastq_files = os.listdir(path_trim_fq)
    pair1 = []
    pair2 = []
    for file in fastq_files:
        if file[-5:] == '_1.fq':
            pair1.append(file)
        elif file[-5:] == '_2.fq':
            pair2.append(file)

    # verify matching correctly and join path and filename
    pair1.sort()
    pair2.sort()
    pattern = re.compile(r'SRR(\d+?)_')
    path_pair1 = []
    path_pair2 = []
    for i in range(len(pair2)):
        pair_1 = pattern.search(pair1[i]).group()[3:-1]
        pair_2 = pattern.search(pair2[i]).group()[3:-1]
        assert pair_1 == pair_2
        path_pair1.append(os.path.join(path_trim_fq, pair1[i]))
        path_pair2.append(os.path.join(path_trim_fq, pair2[i]))

    # run tophat
    os.system('STAR --runThreadN 20 --twopassMode Basic --outSAMstrandField intronMotif --outFilterIntronMotifs '
              'RemoveNoncanonical --genomeDir ' + path_index + ' --readFilesIn ' + ','.join(path_pair1) + ' ' +
              ','.join(path_pair2) + ' --outFileNamePrefix ' + path_output + ' --outSAMtype BAM SortedByCoordinate '
              '--quantMode GeneCounts TranscriptomeSAM')


if __name__ == '__main__':
    start = time()
    star(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time()
    print(end - start)
