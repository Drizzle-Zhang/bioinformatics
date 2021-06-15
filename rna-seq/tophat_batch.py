#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: tophat_batch.py
# @time: 2018/10/22 16:08

import sys
from time import time
import os
import re
import subprocess


def tophat(path_trim_fq, path_index, path_output, process=2, num_sample=6):
    # read fastq files of paired-end sequencing`
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

    # run tophat2
    subprocesses = []
    for i in range(len(pair1)):
        if i % process == 0:
            for sub in subprocesses:
                sub.wait()
            subprocesses = []
        if i >= num_sample:
            break
        subprocesses.append(subprocess.Popen('tophat2 -p 10 -o ' +
                                             os.path.join(path_output, pattern.search(pair1[i]).group()[:-1])
                                             + ' ' + path_index + ' ' + path_pair1[i]
                                             + ' ' + path_pair2[i], shell=True))

    """
    # run tophat
    os.system('tophat2 -p 20 -o ' + path_output + ' ' + path_index + ' ' + ','.join(path_pair1) + ' '
              + ','.join(path_pair2))"""

    return


if __name__ == '__main__':
    start = time()
    tophat(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time()
    print(end - start)
