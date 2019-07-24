#!/usr/bin/env python
# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: download_EMAT.py
# @time: 2019/5/23 21:40

import subprocess
from time import time
import sys
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial


def mv_file(path_tmp, path_out, one_file):
    tmp_file = os.path.join(path_tmp, one_file)
    out_file = os.path.join(path_out, one_file)
    if os.path.exists(tmp_file):
        if not os.path.exists(out_file):
            os.system(f"mv {tmp_file} {out_file}")

    return


def download(list_ebi, path_out, path_tmp, process):
    df_ebi = pd.read_csv(list_ebi, sep='\t')
    list_ebi = np.array(df_ebi['Comment[FASTQ_URI]']).tolist()

    subprocesses = []
    file_fastq = []
    for i, one in enumerate(list_ebi):
        one_file = one.strip().split('/')[-1]
        if i % process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
            pool = Pool(processes=process)
            func_mv = partial(mv_file, path_tmp, path_out)
            pool.map(func_mv, file_fastq)
            pool.close()
        if os.path.exists(os.path.join(path_out, one_file)):
            continue
        else:
            subprocesses.append(
                subprocess.Popen(
                    f"wget -O {os.path.join(path_tmp, one_file)} {one} "
                    f"> {os.path.join(path_tmp, one + '.log')}",
                    shell=True))
            print(one_file)
            file_fastq.append(one_file)

    for sub_process in subprocesses:
        sub_process.wait()
    pool = Pool(processes=process)
    func_mv = partial(mv_file, path_tmp, path_out)
    pool.map(func_mv, file_fastq)
    pool.close()

    return


if __name__ == '__main__':
    start = time()
    download(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    # [1] a tsv-file generated by EBI, the file contains ftp links
    # [2] path saving downloaded fastq.gz files
    # [3] path saving tmp files
    # [4] number of processes
    end = time()
    print(end - start)
    # python36 download_sra.py [1] [2] [3] [4]