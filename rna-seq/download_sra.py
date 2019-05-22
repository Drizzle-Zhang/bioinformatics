#!/usr/bin/env python
# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: download_sra.py
# @time: 2018/10/22 11:37

import subprocess
from time import time
import sys


def download(ls_sra, path_out, process):
    with open(ls_sra, 'r') as fi:
        downlist = []
        for line in fi:
            tmp_line = line.strip().split(',')
            if tmp_line[0] != 'Run':
                downlist.append(tmp_line[0])

    subprocesses = []
    for i, one in enumerate(downlist):
        if i % process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        subprocesses.append(
            subprocess.Popen(
                'fastq-dump --split-3 -O ' + path_out + ' ' + one, shell=True))

    return


if __name__ == '__main__':
    start = time()
    download(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    end = time()
    print(end - start)
    # python36 download_sra.py SraRunInfo.csv ./ > output
