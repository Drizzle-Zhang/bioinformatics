#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_feature.py
# @time: 2020/5/6 18:29

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool, Process
from functools import partial



if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)