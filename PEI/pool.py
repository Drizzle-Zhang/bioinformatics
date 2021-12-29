#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: pool.py
# @time: 2020/3/24 10:43

from time import time
from multiprocessing import Pool


def func(x):
    return x*x


if __name__ == '__main__':
    time_start = time()
    pool = Pool(3)
    list_input = [1, 2, 3, 4, 5, 6]
    # res = [pool.apply_async(func, (i, )) for i in list_input]
    res = pool.map_async(func, list_input)
    pool.close()
    pool.join()
    list_output = res.get()
    time_end = time()
    print(time_end - time_start)
