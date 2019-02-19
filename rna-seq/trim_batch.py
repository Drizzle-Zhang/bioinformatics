#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: trim_batch.py
# @time: 2018/10/23 10:48

import sys
from time import time
import os
import re
import subprocess


def trim_batch(path_fastq, path_output):
    # read fastq files of paired-end sequencing
    fastq_files = os.listdir(path_fastq)
    pair1 = []
    pair2 = []
    for file in fastq_files:
        if file[-8:] == '_1.fastq':
            pair1.append(file)
        elif file[-8:] == '_2.fastq':
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
        path_pair1.append(os.path.join(path_fastq, pair1[i]))
        path_pair2.append(os.path.join(path_fastq, pair2[i]))

    # run trim_galore
    subprocesses = []
    for i in range(len(pair1)):
        subprocesses.append(subprocess.Popen('trim_galore --illumina --stringency 3 --paired ' + path_pair1[i] + ' '
                                             + path_pair2[i] + ' -o ' + path_output, shell=True))

    for sub in subprocesses:
        sub.wait()

    return


if __name__ == '__main__':
    start = time()
    trim_batch(sys.argv[1], sys.argv[2])
    end = time()
    print(end - start)
