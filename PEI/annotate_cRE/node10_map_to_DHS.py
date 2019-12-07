#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: node10_map_to_DHS.py
# @time: 12/6/19 3:50 PM

from time import time
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from subprocess import check_output
from itertools import combinations

if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)