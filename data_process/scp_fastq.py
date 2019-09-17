#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: cp_data.py
# @time: 9/11/19 12:32 AM

from time import time
import os
from multiprocessing import Pool
from functools import partial


# def scp_gzip(path_output, input_file):
#     file_late = os.path.join(path_output, os.path.split(input_file)[1])
#     os.system(
#         f"scp -P 12306 ssh root@10.157.72.101:{input_file} {path_output}")
#     os.system(f"gzip --best {file_late}")
#
#     return


def scp(path_output, input_file):
    file_late = os.path.join(path_output, os.path.split(input_file)[1])
    if not os.path.exists(file_late):
        if not os.path.exists(file_late + '.gz'):
            if os.path.split(input_file)[1] != "cut_read1.fastq":
                os.system(
                    f"scp -P 12306 ssh root@10.157.72.101:"
                    f"{input_file} {path_output}")

    return


def gzip(path_output, input_file):
    file_late = os.path.join(path_output, os.path.split(input_file)[1])
    if os.path.exists(file_late):
        if not os.path.exists(file_late + '.gz'):
            os.system(f"gzip {file_late}")

    return


def main_func(file_list, path_output, num_proc=40):
    with open(file_list, 'r') as r_f:
        list_files = [file.strip() for file in r_f]

    for i in range(len(list_files)//num_proc + 1):
        sub_list_files = list_files[num_proc * i: num_proc * (i + 1)]
        if i == len(list_files)//num_proc:
            sub_list_files = list_files[num_proc * i:]
        pool = Pool(processes=num_proc)
        func_scp = partial(scp, path_output)
        func_gzip = partial(gzip, path_output)
        pool.map(func_scp, sub_list_files)
        pool.map(func_gzip, sub_list_files)
        pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    main_func('/home/disk/public_data/fastq.txt',
              '/home/disk/public_data/yangjingjing',
              num_proc=40)
    time_end = time()
    print(time_end - time_start)
